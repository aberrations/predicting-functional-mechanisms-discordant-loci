# load libraries ----------------------------------------------------------
library(vroom)
library(caroline)
library(stringr)
library(coloc)
library(locuscomparer)
library(ggplot2)

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data')
gene_list <- read.delim('discordant_loci_genes.txt')

# load data ---------------------------------------------------------------
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/raw_data/METSIM')
adipose_eQTL <- vroom('geneQTL_METSIM_n426_adipose_summaryStats_250kb_all.txt')
adipose_eQTL <- subset(adipose_eQTL, ENSG %in% gene_list$Gene.stable.ID.version)

adipose_sQTL <- vroom('spliceQTL_METSIM_n426_adipose_summaryStats_250kb_all.txt')
adipose_sQTL <- subset(adipose_sQTL, Gene %in% unique(adipose_eQTL$Gene))


metsim_eqtl_formatter <- function(eqtl) {
  eqtl$CHR <- sub(":.*","",eqtl$Variant)
  eqtl$POS <- str_extract(eqtl$Variant,":.*_")
  eqtl$POS <- str_remove(eqtl$POS,":")
  eqtl$POS <- str_remove(eqtl$POS,"_")
  N = 426
  degfreedom = N - 2
  tvalue = qt(eqtl$pval/2, df = degfreedom)
  eqtl$SE = abs(eqtl$beta/tvalue)
  
  eqtl_reformatted <- data.frame(
    SNP = eqtl$rsID,
    CHR = eqtl$CHR,
    POS = eqtl$POS,
    BETA = eqtl$beta,
    SE = eqtl$SE,
    P = eqtl$pval,
    EAF = eqtl$EAF,
    N = 426,
    GENE = eqtl$Gene,
    VARIANT = eqtl$Variant)
  return(eqtl_reformatted)
}
adipose_eQTL_reformatted <- metsim_eqtl_formatter(adipose_eQTL)

metsim_sqtl_formatter <- function(sqtl) {
  sqtl$CHR <- sub(":.*", "", sqtl$Variant)
  sqtl$POS <- str_extract(sqtl$Variant,":.*_")
  sqtl$POS <- str_remove(sqtl$POS,":")
  sqtl$POS <- str_remove(sqtl$POS,"_")
  N = 426
  degfreedom = N - 2
  tvalue = qt(sqtl$pval/2, df = degfreedom)
  sqtl$SE = abs(sqtl$beta/tvalue)
  
  sqtl_reformatted <- data.frame(
    SNP = sqtl$rsID,
    CHR = sqtl$CHR,
    POS = sqtl$POS,
    BETA = sqtl$beta,
    SE = sqtl$SE,
    P = sqtl$pval,
    EAF = sqtl$EAF,
    N = 426,
    GENE = sqtl$Gene,
    SPLICE = sqtl$Splice,
    VARIANT = sqtl$Variant)
  return(sqtl_reformatted)
}
adipose_sQTL_reformatted <- metsim_sqtl_formatter(adipose_sQTL)

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data/METSIM')
write.delim(adipose_sQTL_reformatted,'METSIM_adipose_sQTL.txt')
write.delim(adipose_eQTL_reformatted,'METSIM_adipose_eQTL.txt')
gc()


# colocalization ----------------------------------------------------------
colocalize_discordant_loci_eQTL <- function(QTL_directory,GWAS_directory,results_directory,prior){
  setwd(QTL_directory)
  QTL_files <- list.files(pattern = "*eQTL.txt")
  setwd(GWAS_directory)
  GWAS_files <- list.files(path = GWAS_directory, pattern = '*_gwas_in_loci.txt')
  print(paste('Colocalization of QTL data from', length(QTL_files), 'Tissues with',length(GWAS_files),'GWAS.'))
  for (GWAS in GWAS_files) {
    setwd(GWAS_directory)
    GWAS.df <- vroom(GWAS)
    GWAS.df <- subset(GWAS.df, EAF < 1 & EAF > 0 & !duplicated(RSID))
    GWAS_name <- str_remove(GWAS, pattern = '_gwas_in_loci.txt')
    colocalization_summary_file <- data.frame(GENE = NA,
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
      tissue_QTL.df <- subset(tissue_QTL.df, EAF < 1 & EAF > 0)
      tissue_QTL.df <- subset(tissue_QTL.df, !is.na(SE))
      gene_ids <- unique(tissue_QTL.df$GENE)
      for (gene_id in gene_ids) {
        QTL <- subset(tissue_QTL.df, GENE == gene_id)
        QTL <- subset(QTL, !duplicated(SNP))
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
            datasetA = list(beta = GWAS_in_Locus$BETA, MAF = GWAS_in_Locus$EAF, N = GWAS_in_Locus$N, p = GWAS_in_Locus$P, snp = GWAS_in_Locus$RSID, varbeta = (GWAS_in_Locus$SE)^2, type = 'quant')
          }
          if(GWAS_name == 't2dadjbmi'){
            locus_compare_file_name <- "T2DadjBMI_GWAS"
            GWAS_name_locus <- "T2DadjBMI GWAS"
            datasetA = list(beta = GWAS_in_Locus$BETA, MAF = GWAS_in_Locus$EAF, N = GWAS_in_Locus$N, p = GWAS_in_Locus$P, snp = GWAS_in_Locus$RSID, varbeta = (GWAS_in_Locus$SE)^2, s = 0.09, type = 'cc')
          }
          datasetB = list(snp = QTL$SNP, beta = QTL$BETA, MAF = QTL$EAF, N = QTL$N, p = QTL$P, snp = QTL$SNP, varbeta = (QTL$SE)^2, type = 'quant')
          colocalization_file <- coloc.abf(dataset1 = datasetA, dataset2 = datasetB, p12 = prior)
          colocalization_summary <- colocalization_file$summary
          colocalization_results <- colocalization_file$results
          variant_index = which.max(colocalization_results$SNP.PP.H4)
          lead_variant = QTL$SNP[variant_index]
          lead_variant_effect = subset(QTL, SNP == lead_variant)$BETA[1]
          gene_id = unique(QTL$GENE)
          locus_information <- c(gene_id,
                                 lead_variant,
                                 lead_variant_effect,
                                 colocalization_summary)
          PPH4 = colocalization_summary[6]
          colocalization_summary_file <- rbind(colocalization_summary_file, locus_information)
          if(PPH4 > 0.50) {
            locus_compare_gwas <- data.frame(rsid = GWAS_in_Locus$RSID, pval = GWAS_in_Locus$P)
            locus_compare_eQTL <- data.frame(rsid = QTL$SNP, pval = QTL$P)
            locuscompare_plot <- locuscompare(in_fn2 = locus_compare_gwas, in_fn1 = locus_compare_eQTL, title1 = paste(gene_id,"Adipose eQTL"), title2 = GWAS_name_locus, genome = "hg38")
            setwd(results_directory)
            dir.create(path = 'results', showWarnings = FALSE)
            write.delim(df = colocalization_results, file = paste('results/',gene_id,'_',locus_compare_file_name,"_coloc_results.txt",sep = ""))
            setwd(paste(results_directory,"/locuscompareplots",sep = ""))
            ggsave(filename = paste(gene_id,'_',locus_compare_file_name,".jpg",sep = ""),locuscompare_plot)
            setwd(QTL_directory)
          }
          print(paste('Colocalization of',gene_id,'locus is complete.'))
        }
      }
    }
    setwd(results_directory)
    colocalization_summary_file_name <- paste(GWAS_name,'_coloc_summary_file.txt',sep = '')
    colocalization_summary_file <- na.omit(colocalization_summary_file)
    write.delim(colocalization_summary_file,colocalization_summary_file_name)
  }
  print('Colocalization analysis is complete.')
}
colocalize_discordant_loci_eQTL(QTL_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data/METSIM',
                           GWAS_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data/GWAS',
                           results_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/results/METSIM/eQTL',
                           prior = 1e-05)
colocalize_discordant_loci_sQTL <- function(QTL_directory,GWAS_directory,results_directory,prior){
  setwd(QTL_directory)
  QTL_files <- list.files(pattern = "*sQTL.txt")
  setwd(GWAS_directory)
  GWAS_files <- list.files(path = GWAS_directory, pattern = '*_gwas_in_loci.txt')
  print(paste('Colocalization of QTL data from', length(QTL_files), 'Tissues with',length(GWAS_files),'GWAS.'))
  for (GWAS in GWAS_files) {
    setwd(GWAS_directory)
    GWAS.df <- vroom(GWAS)
    GWAS.df <- subset(GWAS.df, EAF < 1 & EAF > 0 & !duplicated(RSID))
    GWAS_name <- str_remove(GWAS, pattern = '_gwas_in_loci.txt')
    colocalization_summary_file <- data.frame(GENE = NA,
                                              SPLICE = NA,
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
      tissue_QTL.df <- subset(tissue_QTL.df, EAF < 1 & EAF > 0)
      tissue_QTL.df <- subset(tissue_QTL.df, !is.na(SE))
      splice_ids <- unique(tissue_QTL.df$SPLICE)
      for (splice_id in splice_ids) {
        QTL <- subset(tissue_QTL.df, SPLICE == splice_id)
        QTL <- subset(QTL, !duplicated(SNP))
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
            datasetA = list(beta = GWAS_in_Locus$BETA, MAF = GWAS_in_Locus$EAF, N = GWAS_in_Locus$N, p = GWAS_in_Locus$P, snp = GWAS_in_Locus$RSID, varbeta = (GWAS_in_Locus$SE)^2, type = 'quant')
          }
          if(GWAS_name == 't2dadjbmi'){
            locus_compare_file_name <- "T2DadjBMI_GWAS"
            GWAS_name_locus <- "T2DadjBMI GWAS"
            datasetA = list(beta = GWAS_in_Locus$BETA, MAF = GWAS_in_Locus$EAF, N = GWAS_in_Locus$N, p = GWAS_in_Locus$P, snp = GWAS_in_Locus$RSID, varbeta = (GWAS_in_Locus$SE)^2, s = 0.09, type = 'cc')
          }
          datasetB = list(beta = QTL$BETA, MAF = QTL$EAF, N = QTL$N, p = QTL$P, snp = QTL$SNP, varbeta = (QTL$SE)^2, type = 'quant')
          colocalization_file <- coloc.abf(dataset1 = datasetA, dataset2 = datasetB, p12 = prior)
          colocalization_summary <- colocalization_file$summary
          colocalization_results <- colocalization_file$results
          variant_index = which.max(colocalization_results$SNP.PP.H4)
          lead_variant = QTL$SNP[variant_index]
          lead_variant_effect = subset(QTL, SNP == lead_variant)$BETA[1]
          gene_id = unique(QTL$GENE)
          
          locus_information <- c(gene_id,
                                 splice_id,
                                 lead_variant,
                                 lead_variant_effect,
                                 colocalization_summary)
          PPH4 = colocalization_summary[6]
          colocalization_summary_file <- rbind(colocalization_summary_file, locus_information)
          if(PPH4 > 0.50) {
            locus_compare_gwas <- data.frame(rsid = GWAS_in_Locus$RSID, pval = GWAS_in_Locus$P)
            locus_compare_eQTL <- data.frame(rsid = QTL$SNP, pval = QTL$P)
            locuscompare_plot <- locuscompare(in_fn2 = locus_compare_gwas, in_fn1 = locus_compare_eQTL, title1 = paste(unique(QTL$GENE),splice_id,"sQTL"), title2 = GWAS_name_locus, genome = "hg38")
            setwd(results_directory)
            dir.create(path = 'results', showWarnings = FALSE)
            write.delim(df = colocalization_results, file = paste('results/',unique(QTL$GENE),'_',splice_id,'_',locus_compare_file_name,"_coloc_results.txt",sep = ""))
            setwd(paste(results_directory,"/locuscompareplots",sep = ""))
            ggsave(filename = paste(unique(QTL$GENE),'_',splice_id,'_',locus_compare_file_name,".jpg",sep = ""),locuscompare_plot)
            setwd(QTL_directory)
          }
          print(paste('Colocalization of',splice_id,'locus is complete.'))
        }
      }
    }
    setwd(results_directory)
    colocalization_summary_file_name <- paste(GWAS_name,'_coloc_summary_file.txt',sep = '')
    colocalization_summary_file <- na.omit(colocalization_summary_file)
    write.delim(colocalization_summary_file,colocalization_summary_file_name)
  }
  print('Colocalization analysis is complete.')
}
colocalize_discordant_loci_sQTL(QTL_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data/METSIM',
                                GWAS_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data/GWAS',
                                results_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/results/METSIM/sQTL',
                                prior = 1e-05)

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/results/METSIM/eQTL')
T2D_eQTL <- read.delim('t2dadjbmi_coloc_summary_file.txt')
WHR_eQTL <- read.delim('whradjbmi_coloc_summary_file.txt')

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/results/METSIM/sQTL')
T2D_sQTL <- read.delim('t2dadjbmi_coloc_summary_file.txt')
WHR_sQTL <- read.delim('whradjbmi_coloc_summary_file.txt')

