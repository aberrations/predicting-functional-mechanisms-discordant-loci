#### format input files for TACTICAL, run tissue classification analysis, and plot results ------------------------------------------------

#### load libraries ------------------------------------------------

library(vroom)
library(TACTICAL)
library(caroline)
library(ggplot2)
library(gridExtra)
library(reshape2)


#### format credible-input text file ------------------------------------------------
## declare functions
{
  # do ld_clumping in GWAS, and then annotate signal number. format credible_input.txt file
  credible_input_formatter <- function(credible_set_directory, GWAS, output_directory)
  {
    linkage_disequilibrium_threshold = 0.5
    credible_input <- data.frame(SIGNAL = NA,
                                 SNP = NA,
                                 CHR = NA,
                                 POS = NA,
                                 VALUE = NA)
    setwd(credible_set_directory)
    credible_set_file_list <- list.files()
    
    for (credible_set_file in credible_set_file_list) {
      credible_set_locus <- read.delim(credible_set_file)
      snp_list <- credible_set_locus$SNP
      toxic_snps <- c('rs60244988','rs61115731')
      snp_list <- subset(snp_list, !(snp_list %in% toxic_snps))
      credible_set_locus$SIGNAL <- NA
      signal_number = 0
      while(length(snp_list) > 0 ) {
        lead_snp = snp_list[1]
        signal_number = signal_number + 1
        snps_in_ld <- haploR::queryHaploreg(query = lead_snp, ldThresh = linkage_disequilibrium_threshold)
        snp_list <- subset(snp_list, !(snp_list %in% snps_in_ld$rsID))
        credible_set_signal <- subset(credible_set_locus, SNP %in% snps_in_ld$rsID)
        locus_coordinates <- unique(credible_set_signal$COORDINATES)[1]
        credible_set_signal$SIGNAL <- paste(locus_coordinates,'_',signal_number,sep='')
        credible_set_signal <- na.omit(credible_set_signal)
        credible_input_GWAS <- subset(GWAS, SNP %in% credible_set_signal$SNP)
        credible_set_signal <- credible_set_signal[(order(credible_set_signal$SNP)),]
        credible_input_GWAS <- credible_input_GWAS[(order(credible_input_GWAS$SNP)),]
        credible_set_signal <- credible_set_signal[(!duplicated(credible_set_signal$SNP)),]
        credible_input_GWAS <- credible_input_GWAS[(!duplicated(credible_input_GWAS$SNP)),]
        SIGNALs <- credible_set_signal$SIGNAL
        SNPs <- credible_set_signal$SNP
        CHRs <- paste('chr',credible_input_GWAS$CHR, sep = '')
        POSs <- credible_input_GWAS$POS
        VALUEs <- credible_set_signal$BF
        credible_input_signal <- data.frame(SIGNAL = SIGNALs,
                                            SNP = SNPs,
                                            CHR = CHRs,
                                            POS = POSs,
                                            VALUE = VALUEs)
        credible_input <- rbind(credible_input,credible_input_signal)
        write.delim(credible_input,file = paste(output_directory,'/credible_input.txt',sep = ''))
      }
    }
    credible_input <- subset(credible_input, !duplicated(SNP) & !is.na(SNP))
    return(credible_input)
  }
  
  # reformat tissue bed files and write tissue path file for TACTICAL
  tissue_path_formatter <- function(tissue_bed_directory,output_directory){
    setwd(tissue_bed_directory)
    tissues <- list.files()
    for(tissue in tissues) {
      setwd(tissue_bed_directory)
      tissue_bed <- read.delim(file = tissue, header = F)
      tissue_bed$V5 <- NULL
      setwd(output_directory)
      write.table(tissue_bed,file = tissue, row.names = F, col.names = F, quote = F)
    }
    setwd(output_directory)
    
    tissue_path_file <- data.frame(tissue_name = gsub("\\..*","",tissues),
                                   file_name = tissues)
    tissue_path_file
    write.table(tissue_path_file, file = 'tissue_path.txt', row.names = F, col.names = F, quote = F)
  }
  tissue_path_formatter(tissue_bed_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/raw_data/tissue_annotation',
                        output_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/formatted_data/tissue_annotation')
  
  # downloaded coding region file and expression specificity scores from TACTICAL
  # did fGWAS tissue annotation analysis in terminal
}

## call functions
{
  #notes on manual changes.
  #SNPs rs60644988 and rs61115731 were removed and manually formatted because of HaploReg issues.
  
  WHRadjBMI <- vroom('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/formatted_data/WHRadjBMI_GWAS.txt')
  T2DadjBMI <- vroom('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/formatted_data/T2DadjBMI_GWAS.txt')
  credible_99 <- vroom('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/credible_set/results/credible_sets/credible_99.txt')
  credible_99 <- credible_99[(order(credible_99$SNP)),]
  WHRadjBMI <- subset(WHRadjBMI, SNP %in% credible_99$SNP)
  WHRadjBMI <- WHRadjBMI[(order(WHRadjBMI$SNP)),]
  #T2DadjBMI <- subset(T2DadjBMI, SNP %in% credible_99$SNP)
  
  credible_set_input <- credible_input_formatter(credible_set_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/credible_set/results/credible_sets/99',
                                                 GWAS = WHRadjBMI,
                                                 output_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/formatted_data')
  credible_input.df <- vroom('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/formatted_data/credible_input.txt')
  credible_input.df <- subset(credible_input.df, !duplicated(SNP) & !is.na(SNP))
  write.delim(credible_input.df,'/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/formatted_data/credible_input_corrected.txt')
}

#### execute TACTICAL ------------------------------------------------
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/formatted_data')
Combined_snp.df <- annotate_snps(snp_file = "credible_input_corrected.txt",
                                 tissue_path_file = "tissue_path.txt",
                                 tissue_annotation_file = "fGWAS_enrichment_Combined_Tissues.txt",
                                 genomic_path_file = "genomic_file_paths.txt",
                                 genomic_annotation_file = "fGWAS_enrichment_Combined_Coding.txt")

Combined_tvec.df <- calculate_tissue_vectors(snp.annotated.df = Combined_snp.df,
                                             tissue_annotation_file = "fGWAS_enrichment_Combined_Tissues.txt",
                                             genomic_annotation_file = "fGWAS_enrichment_Combined_Coding.txt",
                                             ess.annot = "Coding",
                                             ess.file = "gene-expression-specificity-scores.txt")

Combined_tscores.df <- calculate_toa_scores(snp.tissvec.df = Combined_tvec.df)

Combined_classified.df <- tissue_classifier(toa.df=Combined_tscores.df)
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/results')
write.delim(Combined_classified.df, 'Combined_classified.txt')
write.delim(Combined_tscores.df, 'Combined_tscores.txt')

#### make TACTICAL figures ------------------------------------------------
Combined_tscores.df <- read.delim('Combined_tscores.txt')
Combined_classified.df <- read.delim('Combined_classified.txt')
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/figures')
pr.out <- prcomp(dplyr::select(Combined_tscores.df,-one_of("SIGNAL")),scale=TRUE)
pca.df <- as.data.frame(pr.out$x)
pca.df$classification <- Combined_classified.df$classification
pltA <- ggplot(data=dplyr::filter(pca.df,classification!="unclassified"),
               aes(x=PC1,y=PC2,fill=classification)) +
  geom_point(shape=21,alpha=0.8,color="black",size=2) + 
  scale_fill_brewer(palette = "Set1",name="Assigned Tissue") + 
  theme_classic()
pltB <- ggplot(data=dplyr::filter(pca.df,classification!="unclassified"), 
               aes(x=PC2,y=PC3,fill=classification)) +
  geom_point(shape=21,alpha=0.8,color="black",size=2) + 
  scale_fill_brewer(palette = "Set1",name="Assigned Tissue") +
  theme_classic()
p <- grid.arrange(pltA,pltB,nrow=2)
ggsave(plot = p, filename = 'Combined_classified_PCA.jpg', units = 'in', height = 7.5/2, width = 7.5)

#### map discordant genetic association signals to epigenetic annotations ------------------------------------------------
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/formatted_data')
credible_input_corrected.df <- vroom('credible_input_corrected.txt')
credible_99 <- vroom('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/credible_set/formatted_data/credible_99.txt')
discordant_snps_99 <- subset(credible_99, sign(WHRadjBMI_BETA) != sign(T2DadjBMI_BETA))

discordant_genetic_signals <- subset(credible_input_corrected.df, SNPID %in% discordant_snps_99$SNP)
discordant_genetic_signals_mapped <- subset(Combined_snp.df, SNPID %in% discordant_genetic_signals$SNPID)
discordant_genetic_signals_mapped <- discordant_genetic_signals_mapped[,-c(1:5)]
discordant_genetic_signals_mapped_column_sums <- colSums(discordant_genetic_signals_mapped)
discordant_genetic_signals_mapped <- data.frame(SNPs_in_annotation = colSums(discordant_genetic_signals_mapped))
discordant_genetic_signals_mapped <- subset(discordant_genetic_signals_mapped, SNPs_in_annotation > 0)
discordant_genetic_signals_mapped <- data.frame(Annotation = rownames(discordant_genetic_signals_mapped),
                                                SNPs_in_annotation = discordant_genetic_signals_mapped$SNPs_in_annotation)
discordant_genetic_signals_mapped <- discordant_genetic_signals_mapped[(order(discordant_genetic_signals_mapped$SNPs_in_annotation, decreasing = TRUE)),]
discordant_genetic_signals_mapped$Annotation <- str_replace_all(discordant_genetic_signals_mapped$Annotation,
                                                                pattern = ".6_Weak_transcription",
                                                                replacement = " TxWk")
discordant_genetic_signals_mapped$Annotation <- str_replace_all(discordant_genetic_signals_mapped$Annotation,
                                                                pattern = ".9_Active_enhancer_1",
                                                                replacement = " EnhA1")
discordant_genetic_signals_mapped$Annotation <- str_replace_all(discordant_genetic_signals_mapped$Annotation,
                                                                pattern = ".10_Active_enhancer_2",
                                                                replacement = " EnhA2")
discordant_genetic_signals_mapped$Annotation <- str_replace_all(discordant_genetic_signals_mapped$Annotation,
                                                                pattern = ".5_Strong_transcription",
                                                                replacement = " TxSt")
discordant_genetic_signals_mapped$Annotation <- str_replace_all(discordant_genetic_signals_mapped$Annotation,
                                                                pattern = ".11_Weak_enhancer",
                                                                replacement = " EnhWk")
discordant_genetic_signals_mapped$Annotation <- str_replace_all(discordant_genetic_signals_mapped$Annotation,
                                                                pattern = ".17_Weak_repressed_polycomb",
                                                                replacement = " PrcWk")
discordant_genetic_signals_mapped <- subset(discordant_genetic_signals_mapped, Annotation != "Coding")
write.delim(discordant_genetic_signals_mapped, 'discordant_genetic_signals_mapped.txt')
ggplot(data = discordant_genetic_signals_mapped, aes(x = reorder(Annotation, -SNPs_in_annotation), y = SNPs_in_annotation)) + 
  geom_bar(stat = "identity") +
  xlab('Annotation') + 
  ylab('SNPs in Tissue Annotation') + 
  coord_flip()
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/figures')
ggsave(filename = 'discordant_genetic_signals_mapped.png', units = 'in', height = 5, width = 5)

discordant_genetic_signals_mapped_10 <- discordant_genetic_signals_mapped[1:10,]
ggplot(data = discordant_genetic_signals_mapped_10, aes(x = reorder(Annotation, -SNPs_in_annotation), y = SNPs_in_annotation)) + 
  geom_bar(stat = "identity") +
  xlab('Annotation') + ylab('SNPs in Tissue Annotation') + coord_flip()
ggsave(filename = 'discordant_genetic_signals_mapped_10.png', units = 'in', height = 5, width = 5)

#### discordant genetic association signals plotting ------------------------------------------------
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/credible_set/formatted_data')
T2DadjBMI <- read.delim('T2DadjBMI_GWAS_credible.txt')
WHRadjBMI <- read.delim('WHRadjBMI_GWAS_credible.txt')
credible_99 <- read.delim('credible_99.txt')

discordant_snps_99 <- subset(WHRadjBMI, sign(BETA) != sign(T2DadjBMI$BETA))$SNP
discordant_snps_99 <- subset(credible_99, SNP %in% discordant_snps_99)
index <- sign(discordant_snps_99$WHRadjBMI_BETA) == -1
discordant_snps_99$WHRadjBMI_BETA[index] <- discordant_snps_99$WHRadjBMI_BETA[index]*-1
index <- sign(discordant_snps_99$T2DadjBMI_BETA) == 1
discordant_snps_99$T2DadjBMI_BETA[index] <- discordant_snps_99$T2DadjBMI_BETA[index]*-1
unique_signals <- unique(discordant_snps_99$COORDINATES)
representative_snps <- discordant_snps_99[1,]
for (signal in unique_signals) {
  signal_credible <- subset(discordant_snps_99, COORDINATES == signal)
  representative_snp <- which.max(signal_credible$BF)
  signal_credible <- signal_credible[representative_snp,]
  representative_snps <- rbind(representative_snps, signal_credible)
}
rm(signal_credible,representative_snp,unique_signals,signal,index)
representative_snps <- representative_snps[-c(1),]
Combined_discordant_signal_classification <- subset(Combined_classified.df, SIGNAL %in% discordant_genetic_signals$SIGNAL)
Combined_discordant_signal_t_scores <- subset(Combined_tscores.df, SIGNAL %in% discordant_genetic_signals$SIGNAL)
representative_snps$COORDINATES <- paste(representative_snps$COORDINATES,'_1',sep = '')
Combined_discordant_signal_t_scores <- Combined_discordant_signal_t_scores[order(match(x = Combined_discordant_signal_t_scores$SIGNAL,
                                                                                  table = representative_snps$COORDINATES)),]
Combined_discordant_signal_t_scores$SIGNAL <- representative_snps$SNP
# Plot Combined classifications
melted_Combined_discordant_signal_t_scores <- melt(Combined_discordant_signal_t_scores)
ggplot(melted_Combined_discordant_signal_t_scores, aes(fill=variable, y=value, x=SIGNAL)) + 
  geom_bar(position="fill", stat="identity") +
  ylab("Tissue of Action Score") +
  xlab("Representative SNP")
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/figures')
ggsave(filename = 'Combined_discordant_signal_t_scores.png', units = 'in', height = 5, width = 7.5)
