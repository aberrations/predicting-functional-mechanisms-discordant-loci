
# load libraries ----------------------------------------------------------

library(vroom)
library(caroline)
library(stringr)
library(coloc)
library(locuscomparer)
library(ggplot2)


# load data ----------------------------------------------------------
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/raw_data/STARNET')
maf <- read.delim('MAF.txt')
aor <- read.delim('AORlist.txt')
vaf <- read.delim('VAFlist.txt')
liv <- read.delim('LIVlist.txt')
mam <- read.delim('MAMlist.txt')
skm <- read.delim('SKMlist.txt')
suf <- read.delim('SUFlist.txt')

# reformat STARNET files  ----------------------------------------------------------

se_from_t <- function(starnet_file) {
  starnet_file$se <- starnet_file$beta / starnet_file$t.stat
  return(starnet_file)
}

starnet_eQTL_reformatter <- function(starnet_file, MAF_file, Sample_Size) {
  starnet_file$MAF <- NA
  genes <- unique(starnet_file$gene)
  new_starnet_file <- starnet_file[1,]
  for (ensembl_gene in genes) {
    gene_starnet_file <- subset(starnet_file, gene == ensembl_gene)
    gene_starnet_file <- subset(gene_starnet_file, SNP %in% MAF_file$SNP)
    MAF_file_gene <- subset(MAF_file, SNP %in% gene_starnet_file$SNP)
    gene_starnet_file <- gene_starnet_file[!duplicated(gene_starnet_file$SNP),]
    MAF_file_gene <- MAF_file_gene[!duplicated(MAF_file_gene$SNP),]
    MAF_file_gene <- MAF_file_gene[(order(MAF_file_gene$SNP)),]
    gene_starnet_file <- gene_starnet_file[(order(gene_starnet_file$SNP)),]
    gene_starnet_file$MAF <- MAF_file_gene$MAF
    new_starnet_file <- rbind(new_starnet_file,gene_starnet_file)
  }
  new_starnet_file <- na.omit(new_starnet_file)
  
  
  new_starnet_file <- se_from_t(new_starnet_file)
  new_starnet_file <- data.frame(SNP = new_starnet_file$snpid,
                             CHR = new_starnet_file$chr,
                             POS = new_starnet_file$pos,
                             BETA = new_starnet_file$beta,
                             SE = new_starnet_file$se,
                             P = new_starnet_file$p.value,
                             MAF = new_starnet_file$MAF,
                             N = Sample_Size,
                             P_adj_FDR = new_starnet_file$padj_fdr,
                             GENE = new_starnet_file$gene)
  return(new_starnet_file)
}

vaf = starnet_eQTL_reformatter(starnet_file = vaf, MAF_file = maf, Sample_Size = 509)
aor = starnet_eQTL_reformatter(starnet_file = aor, MAF_file = maf, Sample_Size = 515)
mam = starnet_eQTL_reformatter(starnet_file = mam, MAF_file = maf, Sample_Size = 520)
liv = starnet_eQTL_reformatter(starnet_file = liv, MAF_file = maf, Sample_Size = 523)
skm = starnet_eQTL_reformatter(starnet_file = skm, MAF_file = maf, Sample_Size = 514)
suf = starnet_eQTL_reformatter(starnet_file = suf, MAF_file = maf, Sample_Size = 550)

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data/STARNET')
write.delim(vaf, 'VAF.txt')
write.delim(aor, 'AOR.txt')
write.delim(mam, 'MAM.txt')
write.delim(liv, 'LIV.txt')
write.delim(skm, 'SKM.txt')
write.delim(suf, 'SUF.txt')


# declare colocalization functions -------------------------------------------------------

colocalize_discordant_loci <- function(QTL_directory,GWAS_directory,results_directory,prior){
  setwd(QTL_directory)
  QTL_files <- list.files(pattern = "*.txt")
  setwd(GWAS_directory)
  GWAS_files <- list.files(path = GWAS_directory, pattern = '*_gwas_in_loci.txt')
  print(paste('Colocalization of QTL data from', length(QTL_files), 'Tissues with',length(GWAS_files),'GWAS.'))
  for (GWAS in GWAS_files) {
    setwd(GWAS_directory)
    GWAS.df <- vroom(GWAS)
    GWAS.df <- subset(GWAS.df, EAF < 1 & EAF > 0)
    GWAS_name <- str_remove(GWAS, pattern = '_gwas_in_loci.txt')
    colocalization_summary_file <- data.frame(Tissue_Name = NA,
                                              Gene_ID = NA,
                                              lead_snp = NA,
                                              lead_snp_effect = NA,
                                              nsnps = NA,
                                              PPH0 = NA,
                                              PPH1 = NA,
                                              PPH2 = NA,
                                              PPH3 = NA,
                                              PPH4 = NA)
    setwd(QTL_directory)
    for (tissue_QTL in QTL_files) {
      tissue_QTL.df <- vroom(tissue_QTL)
      tissue_name <- str_remove(tissue_QTL, pattern = '.txt')
      tissue_QTL.df <- subset(tissue_QTL.df, !is.na(SE))
      gene_ids <- unique(tissue_QTL.df$GENE)
      for (gene_id in gene_ids) {
        QTL <- subset(tissue_QTL.df, GENE == gene_id)
        start_position <- min(QTL$POS)
        stop_position <- max(QTL$POS)
        chromosome <- unique(QTL$CHR)
        GWAS_in_Locus <- subset(GWAS.df, POS > start_position & POS < stop_position & CHR == chromosome)
        QTL <- QTL[(QTL$SNP %in% GWAS_in_Locus$RSID),]
        GWAS_in_Locus <- GWAS_in_Locus[(GWAS_in_Locus$RSID %in% QTL$SNP),]
        if (nrow(GWAS_in_Locus) > 0 & nrow(QTL) > 0){
          if(GWAS_name == 'whradjbmi'){
            locus_compare_file_name <- "WHRadjBMI_GWAS"
            GWAS_name_locus <- "WHRadjBMI GWAS"
            datasetA = list(beta = GWAS_in_Locus$BETA, MAF = GWAS_in_Locus$EAF, N = GWAS_in_Locus$N, p = GWAS_in_Locus$P, SNP = GWAS_in_Locus$RSID, varbeta = (GWAS_in_Locus$SE)^2, type = 'quant')
          }
          if(GWAS_name == 't2dadjbmi'){
            locus_compare_file_name <- "T2DadjBMI_GWAS"
            GWAS_name_locus <- "T2DadjBMI GWAS"
            datasetA = list(beta = GWAS_in_Locus$BETA, MAF = GWAS_in_Locus$EAF, N = GWAS_in_Locus$N, p = GWAS_in_Locus$P, SNP = GWAS_in_Locus$RSID, varbeta = (GWAS_in_Locus$SE)^2, s = 0.09, type = 'cc')
          }
          datasetB = list(beta = QTL$BETA, MAF = QTL$MAF, N = QTL$N, p = QTL$P, SNP = QTL$SNP, varbeta = (QTL$SE)^2, type = 'quant')
          colocalization_file <- coloc.abf(dataset1 = datasetA, dataset2 = datasetB, p12 = prior)
          colocalization_file_summary <- colocalization_file$summary
          colocalization_results <- colocalization_file$results
          variant_index = which.max(colocalization_results$SNP.PP.H4)
          lead_variant = QTL$SNP[variant_index]
          lead_variant_effect = subset(QTL, SNP == lead_variant)$BETA[1]
          locus_information <- c(tissue_name,
                                 gene_id,
                                 lead_variant,
                                 lead_variant_effect,
                                 colocalization_file_summary)
          PPH4 = colocalization_file_summary[6]
          colocalization_summary_file <- rbind(colocalization_summary_file, locus_information)
          if(PPH4 > 0.50) {
            locus_compare_gwas <- data.frame(rsid = GWAS_in_Locus$RSID, pval = GWAS_in_Locus$P)
            locus_compare_eQTL <- data.frame(rsid = QTL$SNP, pval = QTL$P)
            locuscompare_plot <- locuscompare(in_fn2 = locus_compare_gwas, in_fn1 = locus_compare_eQTL, title1 = paste(tissue_name,gene_id,"eQTL"), title2 = GWAS_name_locus, genome = "hg38")
            setwd(results_directory)
            dir.create('results', showWarnings = F)
            write.delim(colocalization_results, file = paste('results/',tissue_name,'_',gene_id,'_',locus_compare_file_name,"_coloc_summary.txt",sep = ""))
            setwd(paste(results_directory,"/locuscompareplots",sep = ""))
            ggsave(filename = paste(tissue_name,'_',gene_id,'_',locus_compare_file_name,".jpg",sep = ""),locuscompare_plot)
            setwd(QTL_directory)
          }
          print(paste('Colocalization of',gene_id,'locus is complete.'))
        }
      }
    }
    setwd(results_directory)
    colocalization_summary_file_name <- paste(GWAS_name,'_coloc_summary_file.txt',sep = '')
    colocalization_summary_file <- subset(colocalization_summary_file, Tissue_Name != 'character')
    write.delim(colocalization_summary_file,colocalization_summary_file_name)
  }
  print('Colocalization analysis is complete.')
}

colocalize_discordant_loci(QTL_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data/STARNET',
                           GWAS_directory = '/Users/ya8eb/Documents/Research/Functional_Annotation_T2D_BFD/results/discordant_loci_annotation/formatted_data/GWAS',
                           results_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/results/STARNET',
                           prior = 1e-05)

setwd('/Users/ya8eb/Documents/Research/Functional_Annotation_T2D_BFD/results/discordant_loci_annotation/results/STARNET/1e-05')
T2D_results <- read.delim('t2dadjbmi_coloc_summary_file.txt')
WHR_results <- read.delim('whradjbmi_coloc_summary_file.txt')

#significant colocs
sig_T2D_results <- subset(T2D_results, PPH4 > 0.40)
sig_WHR_results <- subset(WHR_results, PPH4 > 0.40)
shared <- subset(sig_T2D_results, Gene_ID %in% sig_WHR_results$Gene_ID)
shared <- unique(shared$Gene_ID)
