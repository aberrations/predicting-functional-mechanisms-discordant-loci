
# load libraries ----------------------------------------------------------
library(coloc)
library(caroline)
library(vroom)
library(stringr)
library(ggplot2)
library(locuscomparer)
library(tidyverse)
library(qdapTools)


# clean eQTL files --------------------------------------------------------
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data')
genes_in_discordant_loci <- read.delim('discordant_loci_genes.txt')
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/raw_data/GTEx')
lookup_table <- as.data.frame(vroom('GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz'))

clean_GTEx_eQTL_file <- function(GTEx_file_directory,GTEx_eQTL_file_name, lookup_table, gene_list,formatted_file_directory) {
  gc()
  setwd(GTEx_file_directory)
  GTEx_eQTL <- as.data.frame(vroom(GTEx_eQTL_file_name))
  GTEx_eQTL <- subset(GTEx_eQTL, gene_id %in% gene_list$Gene.stable.ID.version)
  gc()
  GTEx_eQTL <- merge(GTEx_eQTL,lookup_table, by = "variant_id")
  GTEx_eQTL <- data.frame(
    GENE = GTEx_eQTL$gene_id,
    SNP = GTEx_eQTL$rs_id_dbSNP151_GRCh38p7,
    CHR = (str_replace(GTEx_eQTL$chr,'chr','')),
    POS = GTEx_eQTL$variant_pos,
    REF = GTEx_eQTL$ref,
    ALT = GTEx_eQTL$alt,
    BETA = GTEx_eQTL$slope,
    VARBETA = (GTEx_eQTL$slope_se)^2,
    SE = GTEx_eQTL$slope_se,
    EAF = GTEx_eQTL$maf,
    N = GTEx_eQTL$ma_count,
    P = GTEx_eQTL$pval_nominal)
  setwd(formatted_file_directory)
  write.delim(GTEx_eQTL,str_remove(GTEx_eQTL_file_name,'.gz'))
  gc()
}

list_of_eQTL_files <- c('Adipose_Subcutaneous.allpairs.txt.gz',
                        'Adipose_Visceral_Omentum.allpairs.txt.gz',
                        'Liver.allpairs.txt.gz',
                        'Muscle_Skeletal.allpairs.txt.gz',
                        'Pancreas.allpairs.txt.gz')
clean_GTEx_eQTL_file(GTEx_file_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/raw_data/GTEx',
                     GTEx_eQTL_file_name = list_of_eQTL_files[5],
                     lookup_table = lookup_table,
                     gene_list = genes_in_discordant_loci,
                     formatted_file_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data/GTEx')


for (eQTL_file in list_of_eQTL_files) {
  clean_GTEx_eQTL_file(GTEx_file_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/raw_data/GTEx',
                       GTEx_eQTL_file_name = eQTL_file,
                       lookup_table = lookup_table,
                       gene_list = genes_in_discordant_loci,
                       formatted_file_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data/GTEx')
}

# declare functions -------------------------------------------------------

colocalize_discordant_loci <- function(QTL_directory,GWAS_directory,results_directory,prior){
  setwd(QTL_directory)
  QTL_files <- list.files(pattern = "*.txt")
  setwd(GWAS_directory)
  GWAS_files <- list.files(pattern = '*_gwas_in_loci.txt')
  print(paste('Colocalization of QTL data from', length(QTL_files), 'Tissues with',length(GWAS_files),'GWAS.'))
  for (GWAS in GWAS_files) {
    setwd(GWAS_directory)
    GWAS.df <- vroom(GWAS)
    GWAS_name <- str_remove(GWAS, pattern = '_gwas_in_loci.txt')
    GWAS.df <- subset(GWAS.df, EAF < 1 & EAF > 0 & !duplicated(RSID))
    colocalization_summary_file <- data.frame(Tissue_Name = character(0),
                                              Gene_ID = character(0),
                                              lead_snp = character(0),
                                              lead_snp_effect = numeric(0),
                                              lead_snp_info = character(0),
                                              nsnps = numeric(0),
                                              PPH0 = numeric(0),
                                              PPH1 = numeric(0),
                                              PPH2 = numeric(0),
                                              PPH3 = numeric(0),
                                              PPH4 = numeric(0))
    colocalization_summary_file[1,] <- NA
    setwd(QTL_directory)
    for (tissue_QTL in QTL_files) {
      tissue_QTL.df <- vroom(tissue_QTL)
      tissue_QTL.df <- subset(tissue_QTL.df, EAF < 1 & EAF > 0)
      tissue_name <- str_remove(string = tissue_QTL, pattern = '.allpairs.txt')
      tissue_name <- str_replace(string = tissue_name, pattern = '_', replacement = ' ')
      print(paste('Colocalization of eQTLs in',tissue_name,'has begun.'))
      genes <- unique(tissue_QTL.df$GENE)
      for (gene in genes) {
        QTL <- subset(tissue_QTL.df, GENE == gene)
        QTL <- subset(QTL, !duplicated(SNP))
        GWAS_in_Locus <- GWAS.df[(GWAS.df$RSID %in% QTL$SNP),]
        QTL <- QTL[(QTL$SNP %in% GWAS_in_Locus$RSID),]
        start_position <- min(GWAS_in_Locus$POS)
        stop_position <- max(GWAS_in_Locus$POS)
        chromosome <- unique(GWAS_in_Locus$CHR)
        if (nrow(GWAS_in_Locus) > 0 & nrow(QTL) > 0){
          if(GWAS_name == 'whradjbmi') {
            datasetA = list(beta = GWAS_in_Locus$BETA, MAF = GWAS_in_Locus$EAF, N = GWAS_in_Locus$N, p = GWAS_in_Locus$P, snp = GWAS_in_Locus$RSID, varbeta = (GWAS_in_Locus$SE)^2, type = 'quant')
            GWAS_title <- 'WHRadjBMI'
          }
          if(GWAS_name == 't2dadjbmi') {
            datasetA = list(beta = GWAS_in_Locus$BETA, MAF = GWAS_in_Locus$EAF, N = GWAS_in_Locus$N, p = GWAS_in_Locus$P, snp = GWAS_in_Locus$RSID, varbeta = (GWAS_in_Locus$SE)^2, s = 0.09, type = 'cc')
            GWAS_title <- 'T2DadjBMI'
          }
          if(exists('datasetA')){
            datasetB = list(beta = QTL$BETA, MAF = QTL$EAF, N = QTL$N, p = QTL$P, snp = QTL$SNP, varbeta = QTL$VARBETA, type = 'quant')
            colocalization_file <- coloc.abf(dataset1 = datasetA, dataset2 = datasetB, p12 = prior)
            colocalization_results <- colocalization_file$results
            colocalization_summary = colocalization_file$summary
            variant_index = which.max(colocalization_results$SNP.PP.H4)
            lead_variant = colocalization_results$snp[variant_index]
            lead_variant_effect = subset(QTL, SNP == lead_variant)$BETA
            lead_variant_info = subset(QTL, SNP == lead_variant)$SNP
            locus_information <- c(tissue_name, gene, lead_variant, lead_variant_effect, lead_variant_info, colocalization_file$summary)
            colocalization_summary_file <- rbind(colocalization_summary_file, locus_information)
            if(colocalization_summary[6] > 0.50) {
              QTL_locuscompare <- data.frame(rsid = QTL$SNP, pval = QTL$P)
              GWAS_in_Locus_locuscompare <- data.frame(rsid = GWAS_in_Locus$RSID, pval = GWAS_in_Locus$P)
              locuscompareplot <- locuscompare(in_fn1 = QTL_locuscompare, in_fn2 = GWAS_in_Locus_locuscompare, title1 = gene, title2 = GWAS_title,
                                               genome = 'hg19', population = 'EUR')
              locuscompareplot_filename <- paste('locuscompareplots/',gene,'_',GWAS_title,'_',tissue_name,'_locuscompareplots.jpg',sep = '')
              setwd(results_directory)
              dir.create(path = 'results', showWarnings = FALSE)
              write.delim(colocalization_results, paste('results/',gene,'_',GWAS_title,'_',tissue_name,'_coloc_summary_file.txt',sep = ''))
              ggsave(plot = locuscompareplot, filename = locuscompareplot_filename, height = 4.5, width = 6)
              setwd(QTL_directory)
              gc()
            }
            print(paste('Colocalization of',gene,'locus is complete.'))
          }
        }
      }
      print(paste('Colocalization of eQTLs in',tissue_name,' is complete.'))
    }
    setwd(results_directory)
    colocalization_summary_file_name <- paste(GWAS_name,'_coloc_summary_file.txt',sep = '')
    colocalization_summary_file <- na.omit(colocalization_summary_file)
    write.delim(colocalization_summary_file,colocalization_summary_file_name)
  }
}

# execute functions -------------------------------------------------------

colocalize_discordant_loci(QTL_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data/GTEx',
                           GWAS_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/formatted_data/GWAS',
                           results_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/results/GTEx',
                           prior = 1e-05)

# load data -------------------------------------------------------

setwd('/Users/ya8eb/Documents/Research/Functional_Annotation_T2D_BFD/results/discordant_loci_annotation/results/GTEx')
whradjbmi_coloc_summary_file <- vroom('whradjbmi_coloc_summary_file.txt')
t2dadjbmi_coloc_summary_file <- vroom('t2dadjbmi_coloc_summary_file.txt')
print(unique(c(unique(t2dadjbmi_coloc_summary_file$Gene_ID),unique(whradjbmi_coloc_summary_file$Gene_ID))))
significant_colocs_whradjbmi <- subset(whradjbmi_coloc_summary_file, PPH4 > 0.40)
significant_colocs_t2ddadjbmi <- subset(t2dadjbmi_coloc_summary_file, PPH4 > 0.40)
significant_shared <- significant_colocs_t2ddadjbmi[(significant_colocs_t2ddadjbmi$Gene_ID %in% significant_colocs_whradjbmi$Gene_ID),]
sig_genes <- unique(significant_shared$Gene_ID)
print(sig_genes)
print(unique(c(unique(t2dadjbmi_coloc_summary_file$Gene_ID),unique(whradjbmi_coloc_summary_file$Gene_ID))))
