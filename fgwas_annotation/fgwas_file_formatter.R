#### format input files for fGWAS, run tissue genomic annotation analysis, and plot results ------------------------------------------------

#### load libraries ------------------------------------------------
{
  library(vroom)
  library(TACTICAL)
  library(caroline)
}

#### format fGWAS genomic data input text file ------------------------------------------------
## declare functions
# do ld_clumping in GWAS, and then annotate signal number. format credible_input.txt file
fGWAS_base_formatter <- function(GWAS)
{
  fGWAS_input_file <- data.frame(SNPID = GWAS$SNP,
                                 CHR = GWAS$CHR,
                                 POS = GWAS$POS,
                                 F = GWAS$EAF,
                                 Z = GWAS$Z,
                                 N = GWAS$N)
  fGWAS_input_file <- fGWAS_input_file[(order(fGWAS_input_file$POS)),]
  fGWAS_input_file <- fGWAS_input_file[(order(fGWAS_input_file$CHR)),]
  fGWAS_input_file$CHR <- paste('chr',fGWAS_input_file$CHR, sep = '')

  return(fGWAS_input_file)
}

fGWAS_base_formatter <- function(GWAS,bayesian_factor_directory)
{
  setwd(bayesian_factor_directory)
  bayesian_factors <- list.files()
  first_bayesian_factor <- bayesian_factors[1]
  bayesian_factors <- bayesian_factors[2:length(bayesian_factors)]
  bayesian_factor.df <- read.delim(first_bayesian_factor)
  
  for (bayesian_factor in bayesian_factors) {
    temp_bayesian_factor.df <- read.delim(bayesian_factor)
    bayesian_factor.df <- rbind(temp_bayesian_factor.df,bayesian_factor.df)
    bayesian_factor.df <- subset(bayesian_factor.df, snp %in% GWAS$SNP)
  }
  bayesian_factor.df <- bayesian_factor.df[(order(bayesian_factor.df$snp)),]
  GWAS <- GWAS[(order(GWAS$SNP)),]
  fGWAS_input_file <- data.frame(SNPID = GWAS$SNP,
                                 CHR = GWAS$CHR,
                                 POS = GWAS$POS,
                                 F = GWAS$EAF,
                                 Z = GWAS$Z,
                                 N = GWAS$N,
                                 LNBF = log(bayesian_factor.df$internal.sum.lABF))
  fGWAS_input_file <- fGWAS_input_file[(order(fGWAS_input_file$POS)),]
  fGWAS_input_file <- fGWAS_input_file[(order(fGWAS_input_file$CHR)),]
  fGWAS_input_file$CHR <- paste('chr',fGWAS_input_file$CHR, sep = '')
  
  return(fGWAS_input_file)
}

# load data
WHRadjBMI <- vroom('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/formatted_data/WHRadjBMI_GWAS.txt')
T2DadjBMI <- vroom('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/formatted_data/T2DadjBMI_GWAS.txt')
credible_input.df <- vroom('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/formatted_data/credible_input_corrected.txt')
WHRadjBMI <- subset(WHRadjBMI, SNP %in% credible_input.df$SNPID)
T2DadjBMI <- subset(T2DadjBMI, SNP %in% credible_input.df$SNPID)
WHRadjBMI <- WHRadjBMI[(order(WHRadjBMI$POS)),]
WHRadjBMI <- WHRadjBMI[(order(WHRadjBMI$CHR)),]
snp_order <- match(WHRadjBMI$SNP,credible_input.df$SNPID)
snps_ordered <- credible_input.df[snp_order,]
snps_ordered$SEGNUMBER <- substr(snps_ordered$SIGNAL,1,nchar(snps_ordered$SIGNAL)-2)

WHRadjBMI <- WHRadjBMI[(order(WHRadjBMI$SNP)),]
T2DadjBMI <- T2DadjBMI[(order(T2DadjBMI$SNP)),]
credible_input.df <- credible_input.df[(order(credible_input.df$SNPID)),]
WHRadjBMI$Z <- WHRadjBMI$BETA/(WHRadjBMI$VARBETA)^0.5
T2DadjBMI$Z <- T2DadjBMI$BETA/(T2DadjBMI$VARBETA)^0.5




#create fGWAS input file
fGWAS_input_formatted <- fGWAS_base_formatter(GWAS=T2DadjBMI, 
                                               bayesian_factor_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/results/coloc_abf/bayesian_factors')
fGWAS_input_WHRadjBMI.df <- fGWAS_base_formatter(GWAS = WHRadjBMI)
fGWAS_input_T2DadjBMI.df <- fGWAS_base_formatter(GWAS = T2DadjBMI)
T2DadjBMI_GWAS_Case_rate <- 0.09

fGWAS_input_T2DadjBMI.df$NCASE <- round(T2DadjBMI_GWAS_Case_rate*fGWAS_input_T2DadjBMI.df$N)
fGWAS_input_T2DadjBMI.df$NCONTROL <- round((1-T2DadjBMI_GWAS_Case_rate)*fGWAS_input_T2DadjBMI.df$N)
fGWAS_input_T2DadjBMI.df$N <- NULL
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/fgwas_annotation/formatted_data')
write.delim(fGWAS_input_WHRadjBMI.df,'fGWAS_input_WHRadjBMI.txt')
write.delim(fGWAS_input_T2DadjBMI.df,'fGWAS_input_T2DadjBMI.txt')

## preparing each tissue's annotation file
fgwas_file_preparation <- function(credible_input, tissue_annotation,GWAS,bayesian_factor_directory){
  setwd(bayesian_factor_directory)
  bayesian_factors <- list.files()
  first_bayesian_factor <- bayesian_factors[1]
  bayesian_factors <- bayesian_factors[2:length(bayesian_factors)]
  bayesian_factor.df <- read.delim(first_bayesian_factor)
  
  for (bayesian_factor in bayesian_factors) {
    temp_bayesian_factor.df <- read.delim(bayesian_factor)
    bayesian_factor.df <- rbind(temp_bayesian_factor.df,bayesian_factor.df)
    bayesian_factor.df <- subset(bayesian_factor.df, snp %in% GWAS$SNP)
  }
  bayesian_factor.df <- bayesian_factor.df[(order(bayesian_factor.df$snp)),]
  
  fgwas_file <- data.frame(SNPID = credible_input$SNPID,
                           CHR = credible_input$CHR,
                           POS = credible_input$POS,
                           Z = GWAS$Z,
                           F = GWAS$EAF,
                           N = GWAS$N,
                           SEGNUMBER = credible_input$SIGNAL,
                           LNBF = log(bayesian_factor.df$internal.sum.lABF))
  annotations <- unique(tissue_annotation$V4)
  unique_segments <- unique(fgwas_file$SEGNUMBER)
  SEGMENTNUMBER = 0
  for (unique_segment in unique_segments) {
    SEGMENTNUMBER = SEGMENTNUMBER+1
    index <- fgwas_file$SEGNUMBER == unique_segment
    fgwas_file$SEGNUMBER[index] <- SEGMENTNUMBER
  }
  for (annotation in annotations) {
    variant_annotation <- c()
    annotation_bed = subset(tissue_annotation, V4 == annotation)
    for (variant in fgwas_file$SNPID) {
      snp = subset(fgwas_file, SNPID == variant)
      regions_with_snp <- subset(annotation_bed, V1 == snp$CHR[1] & V2 < snp$POS[1] & V3 > snp$POS[1])
      if(nrow(regions_with_snp) > 0){
        variant_annotation <- append(variant_annotation, 1)
      } else {
        variant_annotation <- append(variant_annotation, 0)
      }
    }
    if(annotation == '14_Bivalent/poised_TSS') {
      annotation = '14_Bivalent_poised_TSS'
    }
    if(annotation == '18_Quiescent/low_signal') {
      annotation = '18_Quiescent_low_signal'
    }
    fgwas_file[[annotation]] <- variant_annotation
  }
  fgwas_file$SEGNUMBER <- as.numeric(fgwas_file$SEGNUMBER)
  fgwas_file <- fgwas_file[(order(fgwas_file$SEGNUMBER)),]
  return(fgwas_file)
}
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/formatted_data')
Adipose.chromatinStates.df <- read.delim(file = 'Adipose.chromatinStates.bed', header = F)
HippocampusMiddle.chromatinStates.df <- read.delim(file = 'HippocampusMiddle.chromatinStates.bed', header = F)
Islets.chromatinStates.df <- read.delim(file = 'Islets.chromatinStates.bed', header = F)
Liver.chromatinStates.df <- read.delim(file = 'Liver.chromatinStates.bed', header = F)
SkeletalMuscle.chromatinStates.df <- read.delim(file = 'SkeletalMuscle.chromatinStates.bed', header = F)
Coding.df <- read.delim(file = 'Coding.bed', header = F)
Coding.df$V4 <- 'Coding'
write.delim(Coding.df, 'Coding.bed', col.names = F)
WHRadjBMI$Z <- WHRadjBMI$BETA/(WHRadjBMI$VARBETA)^0.5
bf_directory <- '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/results/coloc_abf/bayesian_factors'
Adipose_fGWAS_input <- fgwas_file_preparation(credible_input = snps_ordered,
                                              tissue_annotation = Adipose.chromatinStates.df,
                                              GWAS = WHRadjBMI,
                                              bayesian_factor_directory = bf_directory)
Liver_fGWAS_input <- fgwas_file_preparation(credible_input = snps_ordered, 
                                            tissue_annotation = Liver.chromatinStates.df,
                                            GWAS = WHRadjBMI,
                                            bayesian_factor_directory = bf_directory)
SkeletalMuscle_fGWAS_input <- fgwas_file_preparation(credible_input = snps_ordered, 
                                                     tissue_annotation = SkeletalMuscle.chromatinStates.df,
                                                     GWAS = WHRadjBMI,
                                                     bayesian_factor_directory = bf_directory)
Islets_fGWAS_input <- fgwas_file_preparation(credible_input = snps_ordered, 
                                             tissue_annotation = Islets.chromatinStates.df,
                                             GWAS = WHRadjBMI,
                                             bayesian_factor_directory = bf_directory)
Coding_fGWAS_input <- fgwas_file_preparation(credible_input = snps_ordered, 
                                             tissue_annotation = Coding.df,
                                             GWAS = WHRadjBMI,
                                             bayesian_factor_directory = bf_directory)

Adipose_WHRadjBMI_fGWAS_input <- fgwas_file_preparation(snps = snps_ordered, 
                                                        tissue_annotation = Adipose.chromatinStates.df,
                                                        credible_set_signal = credible_input.df$SIGNAL)
HippocampusMiddle_WHRadjBMI_fGWAS_input <- fgwas_file_preparation(snps = snps_ordered, 
                                                        tissue_annotation = HippocampusMiddle.chromatinStates.df,
                                                        credible_set_signal = credible_input.df$SIGNAL)
Islets_WHRadjBMI_fGWAS_input <- fgwas_file_preparation(snps = snps_ordered, 
                                                        tissue_annotation = Islets.chromatinStates.df,                                                         credible_set_signal = credible_input.df$SIGNAL)
Liver_WHRadjBMI_fGWAS_input <- fgwas_file_preparation(snps = snps_ordered, 
                                                        tissue_annotation = Liver.chromatinStates.df,
                                                      credible_set_signal = credible_input.df$SIGNAL)
SkeletalMuscle_WHRadjBMI_fGWAS_input <- fgwas_file_preparation(snps = fGWAS_input_WHRadjBMI.df, 
                                                        tissue_annotation = SkeletalMuscle.chromatinStates.df,
                                                        credible_set_signal = credible_input.df$SIGNAL)
Islets_T2DadjBMI_fGWAS_input <- fgwas_file_preparation(snps = fGWAS_input_T2DadjBMI.df, 
                                                        tissue_annotation = Islets.chromatinStates.df,
                                                        credible_set_signal = credible_input.df$SIGNAL)
HippocampusMiddle_T2DadjBMI_fGWAS_input <- fgwas_file_preparation(snps = fGWAS_input_T2DadjBMI.df, 
                                                                  tissue_annotation = HippocampusMiddle.chromatinStates.df,
                                                                  credible_set_signal = credible_input.df$SIGNAL)
Islets_T2DadjBMI_fGWAS_input <- fgwas_file_preparation(snps = fGWAS_input_T2DadjBMI.df, 
                                                       tissue_annotation = Islets.chromatinStates.df,
                                                       credible_set_signal = credible_input.df$SIGNAL)
Liver_T2DadjBMI_fGWAS_input <- fgwas_file_preparation(snps = fGWAS_input_T2DadjBMI.df, 
                                                      tissue_annotation = Liver.chromatinStates.df,
                                                      credible_set_signal = credible_input.df$SIGNAL)
SkeletalMuscle_T2DadjBMI_fGWAS_input <- fgwas_file_preparation(snps = fGWAS_input_T2DadjBMI.df, 
                                                               tissue_annotation = SkeletalMuscle.chromatinStates.df,
                                                               credible_set_signal = credible_input.df$SIGNAL)

Coding_T2DadjBMI_fGWAS_input <- fgwas_file_preparation(snps = fGWAS_input_T2DadjBMI.df, 
                                                               tissue_annotation = Coding.df,
                                                               credible_set_signal = credible_input.df$SIGNAL)
Coding_WHRadjBMI_fGWAS_input <- fgwas_file_preparation(snps = fGWAS_input_WHRadjBMI.df, 
                                                       tissue_annotation = Coding.df,
                                                       credible_set_signal = credible_input.df$SIGNAL)
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/fgwas_annotation/formatted_data')
write.delim(Adipose_fGWAS_input, 'Adipose_fGWAS_input.txt')
write.delim(Islets_fGWAS_input, 'Islets_fGWAS_input.txt')
write.delim(Liver_fGWAS_input, 'Liver_fGWAS_input.txt')
write.delim(SkeletalMuscle_fGWAS_input, 'SkeletalMuscle_fGWAS_input.txt')
write.delim(Coding_fGWAS_input, 'Coding_fGWAS_input.txt')

write.delim(Adipose_WHRadjBMI_fGWAS_input, 'Adipose_WHRadjBMI_fGWAS_input.txt')
write.delim(HippocampusMiddle_WHRadjBMI_fGWAS_input, 'HippocampusMiddle_WHRadjBMI_fGWAS_input.txt')
write.delim(Islets_WHRadjBMI_fGWAS_input, 'Islets_WHRadjBMI_fGWAS_input.txt')
write.delim(Liver_WHRadjBMI_fGWAS_input, 'Liver_WHRadjBMI_fGWAS_input.txt')
write.delim(SkeletalMuscle_WHRadjBMI_fGWAS_input, 'SkeletalMuscle_WHRadjBMI_fGWAS_input.txt')
write.delim(Adipose_T2DadjBMI_fGWAS_input, 'Adipose_T2DadjBMI_fGWAS_input.txt')
write.delim(HippocampusMiddle_T2DadjBMI_fGWAS_input, 'HippocampusMiddle_T2DadjBMI_fGWAS_input.txt')
write.delim(Islets_T2DadjBMI_fGWAS_input, 'Islets_T2DadjBMI_fGWAS_input.txt')
write.delim(Liver_T2DadjBMI_fGWAS_input, 'Liver_T2DadjBMI_fGWAS_input.txt')
write.delim(SkeletalMuscle_T2DadjBMI_fGWAS_input, 'SkeletalMuscle_T2DadjBMI_fGWAS_input.txt')
write.delim(Coding_WHRadjBMI_fGWAS_input, 'Coding_WHRadjBMI_fGWAS_input.txt')
write.delim(Coding_T2DadjBMI_fGWAS_input, 'Coding_T2DadjBMI_fGWAS_input.txt')




#### reformat fGWAS annotation results ------------------------------------------------
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/fgwas_annotation/results')
fGWAS_enrichment_results <- function(path, GWAS_name){
  setwd(path)
  bad_annotations <- c('ln(lk):',
                       'X-validation',
                       'ridgeparam:')
  tissues <- list.dirs(path, recursive = F)
  tissue_annotation_file <- data.frame(tissue = NA, annotation = NA, enrichment = NA)
  for (tissue in tissues) {
    annotation_path <- paste(tissue,'/',GWAS_name,'/',sep = '')
    tissue_name <- sub('.+/(.+)', '\\1', tissue)
    setwd(annotation_path)
    tissue_GWAS_annotation <- read.delim('complete_model.ridgeparams', sep = ' ', header = F)
    tissue_GWAS_annotation$tissue <- tissue_name
    tissue_GWAS_annotation <- subset(tissue_GWAS_annotation, V1 != 'ridgeparam:' & V1 != 'X-validation' )
    tissue_GWAS_annotation <- data.frame(tissue = tissue_GWAS_annotation$tissue,
                                         annotation = tissue_GWAS_annotation$V1,
                                         enrichment = tissue_GWAS_annotation$V2)
    tissue_annotation_file <- rbind.data.frame(tissue_annotation_file,tissue_GWAS_annotation)
    setwd(path)
  }
  tissue_annotation_file <- subset(tissue_annotation_file, !(annotation %in% bad_annotations) & !is.na(annotation))
  tissue_annotation_file$enrichment <- 2^(as.numeric(tissue_annotation_file$enrichment))
  return(tissue_annotation_file)
}
fGWAS_enrichment_Combined <- fGWAS_enrichment_results('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/fgwas_annotation/results',
                                                       GWAS_name = 'Combined')

fGWAS_enrichment_WHRadjBMI <- fGWAS_enrichment_results('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/fgwas_annotation/results',
                                  GWAS_name = 'WHRadjBMI')
fGWAS_enrichment_T2DadjBMI <- fGWAS_enrichment_results('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/fgwas_annotation/results',
                                  GWAS_name = 'T2DadjBMI')

write.delim(fGWAS_enrichment_Combined,'fGWAS_enrichment_Combined.txt')
write.delim(fGWAS_enrichment_WHRadjBMI,'fGWAS_enrichment_WHRadjBMI.txt')
write.delim(fGWAS_enrichment_T2DadjBMI,'fGWAS_enrichment_T2DadjBMI.txt')
