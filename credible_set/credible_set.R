
# load libraries ----------------------------------------------------------

library(vroom)
library(caroline)
library(pheatmap)
library(ggplot2)
library(forestplot)
library(dplyr)
library(stringr)

# get credible set of SNPS in colocalized loci using colocalization BF --------

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/results')
significant_locus_coordinates <- vroom('significant_loci_overlap.txt')
get_credible_sets <- function(coordinates,credible_threshold, results_directory, input_directory) {
  credible_sets <- data.frame(SNP = NA, BF = NA, COORDINATES = NA)
  for (coordinate in 1:nrow(coordinates)) {
    setwd(input_directory)
    locus_information <- coordinates[coordinate,]
    BF_data <- vroom(paste('chr',locus_information$COORDINATES,'_WHRadjBMI_T2DadjBMI_coloc_abf.txt',sep = ''))
    BF_data$SNP.PP.H4 <- as.numeric(BF_data$SNP.PP.H4)
    BF_data <- BF_data[(order(BF_data$SNP.PP.H4, decreasing = TRUE)),]
    BF_row <- 1
    BF_sum <- BF_data$SNP.PP.H4[BF_row]
    while (credible_threshold > BF_sum) {
      BF_row <- BF_row + 1
      BF_sum <- BF_sum + BF_data$SNP.PP.H4[BF_row]
    }
    credible_set_locus <- BF_data[1:BF_row,]
    credible_set_locus <- data.frame(SNP = credible_set_locus$snp, BF = credible_set_locus$SNP.PP.H4, COORDINATES = locus_information$COORDINATES)
    credible_sets <- rbind(credible_set_locus,credible_sets)
    setwd(results_directory)
    write.delim(credible_set_locus,
                file = paste('chr',locus_information$COORDINATES,'_WHRadjBMI_T2DadjBMI_coloc_abf_credible_set.txt',sep = ''))
  }
  return(na.omit(credible_sets))
}

credible_99 <- get_credible_sets(coordinates = significant_locus_coordinates, 
                                 credible_threshold = 0.99,
                                 results_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/credible_set/results/credible_set/99',
                                 input_directory = '/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/results/coloc_abf/bayesian_factors')

# get clustered heatmap of snps with effect directions --------------------
get_effect_directions <- function(credible_set) {
  ordered_heatmap <- data.frame(T2DadjBMI = credible_set$T2DadjBMI_BETA,
                                WHRadjBMI = credible_set$WHRadjBMI_BETA)
  rownames(ordered_heatmap) <- credible_set$SNP
  
  risk_SNPs <- ordered_heatmap$T2DadjBMI > 0
  ordered_heatmap$T2DadjBMI[risk_SNPs] <- -1*ordered_heatmap$T2DadjBMI[risk_SNPs]
  ordered_heatmap$WHRadjBMI[risk_SNPs] <- -1*ordered_heatmap$WHRadjBMI[risk_SNPs]
  ordered_heatmap <- ordered_heatmap[(order(-ordered_heatmap$WHRadjBMI)),]
  ordered_heatmap <- sign(ordered_heatmap)
  return(ordered_heatmap)
}

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/formatted_data')
WHRadjBMI <- vroom('WHRadjBMI_GWAS.txt')
T2DadjBMI <- vroom('T2DadjBMI_GWAS.txt')
WHRadjBMI <- subset(WHRadjBMI, SNP %in% credible_99$SNP)
T2DadjBMI <- subset(T2DadjBMI, SNP %in% credible_99$SNP)


credible_99 <- credible_99[(order(credible_99$SNP)),]
T2DadjBMI <- T2DadjBMI[(order(T2DadjBMI$SNP)),]
WHRadjBMI <- WHRadjBMI[(order(WHRadjBMI$SNP)),]
credible_99$WHRadjBMI_BETA <- WHRadjBMI$BETA
credible_99$T2DadjBMI_BETA <- T2DadjBMI$BETA
credible_99$WHRadjBMI_SE <- (WHRadjBMI$VARBETA)^0.5
credible_99$T2DadjBMI_SE <- (T2DadjBMI$VARBETA)^0.5
discordant_snps_99 <- subset(credible_99, sign(WHRadjBMI_BETA) != sign(T2DadjBMI_BETA))

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/credible_set/figures')
ordered_heatmap <- get_effect_directions(credible_set = credible_99)
ordered_heatmap <- t(ordered_heatmap)
jpeg(filename = 'credible_set_99_heatmap.jpeg', res = 300, height = 5, width = 7.5, units = 'in')
pheatmap(mat = ordered_heatmap, 
         labels_row = rownames(ordered_heatmap), 
         show_colnames = FALSE, 
         cluster_rows = F,
         cluster_cols = F,
         border_color = 'black')
dev.off()

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/credible_set/figures')
#credible set
ordered_heatmap <- get_effect_directions(credible_set = discordant_snps_99)
ordered_heatmap <- t(ordered_heatmap)
jpeg(filename = 'credible_set_99_discordant_heatmap_no_BMI.jpeg', res = 300, height = 3, width = 5, units = 'in')
pheatmap(mat = ordered_heatmap, labels_row = rownames(ordered_heatmap), show_colnames = FALSE)
dev.off()


# SNP annotations of discordant SNPs --------------------------------------
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/credible_set/results')
discordant_rsids <- discordant_snps_99$SNP
write.delim(discordant_rsids, 'discordant_rsids.txt')
# use the online ensembl variant annotation tool
discordant_variant_annotation <- read.delim('discordant_variant_annotation.txt')
summarize_variant_annotations <- function(variant_annotation_df)
{
  annotations <- unique(variant_annotation_df$Consequence)
  variants <- unique(variant_annotation_df$Uploaded_variation)
  summary <- matrix(ncol = length(annotations), nrow = length(variants))
  colnames(summary) = annotations
  rownames(summary) = variants
  iteration = 0
  for (variant in variants) {
    iteration = iteration + 1
    variant_annotation <- subset(variant_annotation_df, Uploaded_variation == variant)
    variant_annotation <- variant_annotation$Consequence
    variant_annotation <- annotations %in% variant_annotation
    summary[iteration,] <- as.numeric(variant_annotation)
  }
  column_names <- str_replace_all(colnames(summary), pattern = '_', replacement = " ")
  column_names <- str_replace_all(column_names, pattern = 'X3', replacement = "3")
  colnames(summary) <- column_names
  summary <- data.frame(summary)
  return(summary)
}
summary_discordant_variant_annotation <- read.delim('summary_discordant_variant_annotation.txt')

summary_discordant_variant_annotation <- summarize_variant_annotations(discordant_variant_annotation)
colnames(summary_discordant_variant_annotation) <- c('intron variant','non-coding transcript variant','regulatory region variant', 'NMD transcript variant',
                                                     'upstream gene variant','downstream gene variant','TF binding site variant','non-coding transcript exon variant',
                                                     'missense variant','3 prime UTR variant NMD transcript variant','3 prime UTR variant')

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/credible_set/figures')
jpeg(filename = "variant_consequences.jpg",
     units = 'in',
     width = 5,
     height = 7.5,
     res = 300)
pheatmap(summary_discordant_variant_annotation,
         show_rownames = F,
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         #fontsize_row = 1,
         fontsize_col = 7,
         color = c('white','black')
         )
dev.off()
write.delim(summary_discordant_variant_annotation, 'summary_discordant_variant_annotation.txt', row.names = T)

# forest plot ------------------------------------------------------------
WHRadjBMI <- WHRadjBMI[(order(WHRadjBMI$POS)),]
WHRadjBMI <- WHRadjBMI[(order(WHRadjBMI$CHR)),]
T2DadjBMI <- T2DadjBMI[(order(T2DadjBMI$POS)),]
T2DadjBMI <- T2DadjBMI[(order(T2DadjBMI$CHR)),]
credible_99 <- credible_99[order(match(credible_99$SNP,T2DadjBMI$SNP)),]
credible_99$WHRadjBMI_BETA <- WHRadjBMI$BETA
credible_99$WHRadjBMI_SE <- (WHRadjBMI$VARBETA)^0.5
credible_99$T2DadjBMI_BETA <- T2DadjBMI$BETA
credible_99$T2DadjBMI_SE <- (T2DadjBMI$VARBETA)^0.5
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/credible_set/formatted_data')
write.delim(T2DadjBMI,'T2DadjBMI_GWAS_credible.txt')
write.delim(WHRadjBMI,'WHRadjBMI_GWAS_credible.txt')
write.delim(credible_99,'credible_99.txt')
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
rm(signal_credible,representative_snp,unique_signals,signal,index,summary_discordant_variant_annotation,summary_discordant_variant_annotation_plot)
representative_snps <- representative_snps[-c(1),]
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/results')
write.delim(df = representative_snps,'representative_snps.txt')

cochrane_from_rmeta <- structure(list(mean  = c(NA, 
                                                NA, 
                                                representative_snps$WHRadjBMI_BETA[1], 
                                                representative_snps$T2DadjBMI_BETA[1],
                                                representative_snps$WHRadjBMI_BETA[2], 
                                                representative_snps$T2DadjBMI_BETA[2], 
                                                representative_snps$WHRadjBMI_BETA[3], 
                                                representative_snps$T2DadjBMI_BETA[3], 
                                                representative_snps$WHRadjBMI_BETA[4], 
                                                representative_snps$T2DadjBMI_BETA[4], 
                                                representative_snps$WHRadjBMI_BETA[5], 
                                                representative_snps$T2DadjBMI_BETA[5], 
                                                representative_snps$WHRadjBMI_BETA[6], 
                                                representative_snps$T2DadjBMI_BETA[6],
                                                NA, 
                                                NA),
                                      lower = c(NA, 
                                                NA, 
                                                representative_snps$WHRadjBMI_BETA[1]-representative_snps$WHRadjBMI_SE[1], 
                                                representative_snps$T2DadjBMI_BETA[1]-representative_snps$T2DadjBMI_SE[1], 
                                                representative_snps$WHRadjBMI_BETA[2]-representative_snps$WHRadjBMI_SE[2], 
                                                representative_snps$T2DadjBMI_BETA[2]-representative_snps$T2DadjBMI_SE[2],
                                                representative_snps$WHRadjBMI_BETA[3]-representative_snps$WHRadjBMI_SE[3], 
                                                representative_snps$T2DadjBMI_BETA[3]-representative_snps$T2DadjBMI_SE[3],
                                                representative_snps$WHRadjBMI_BETA[4]-representative_snps$WHRadjBMI_SE[4], 
                                                representative_snps$T2DadjBMI_BETA[4]-representative_snps$T2DadjBMI_SE[4],
                                                representative_snps$WHRadjBMI_BETA[5]-representative_snps$WHRadjBMI_SE[5], 
                                                representative_snps$T2DadjBMI_BETA[5]-representative_snps$T2DadjBMI_SE[5],
                                                representative_snps$WHRadjBMI_BETA[6]-representative_snps$WHRadjBMI_SE[6], 
                                                representative_snps$T2DadjBMI_BETA[6]-representative_snps$T2DadjBMI_SE[6],
                                                NA, 
                                                NA),
                                      upper = c(NA, 
                                                NA, 
                                                representative_snps$WHRadjBMI_BETA[1]+representative_snps$WHRadjBMI_SE[1],
                                                representative_snps$T2DadjBMI_BETA[1]+representative_snps$T2DadjBMI_SE[1],
                                                representative_snps$WHRadjBMI_BETA[2]+representative_snps$WHRadjBMI_SE[2],
                                                representative_snps$T2DadjBMI_BETA[2]+representative_snps$T2DadjBMI_SE[2],
                                                representative_snps$WHRadjBMI_BETA[3]+representative_snps$WHRadjBMI_SE[3],
                                                representative_snps$T2DadjBMI_BETA[3]+representative_snps$T2DadjBMI_SE[3],
                                                representative_snps$WHRadjBMI_BETA[4]+representative_snps$WHRadjBMI_SE[4],
                                                representative_snps$T2DadjBMI_BETA[4]+representative_snps$T2DadjBMI_SE[4],
                                                representative_snps$WHRadjBMI_BETA[5]+representative_snps$WHRadjBMI_SE[5], 
                                                representative_snps$T2DadjBMI_BETA[5]+representative_snps$T2DadjBMI_SE[5],
                                                representative_snps$WHRadjBMI_BETA[6]+representative_snps$WHRadjBMI_SE[6], 
                                                representative_snps$T2DadjBMI_BETA[6]+representative_snps$T2DadjBMI_SE[6],
                                                NA, 
                                                NA),
                                      group = c(NA, 
                                                NA, 
                                                "WHRadjBMI", 
                                                "T2DadjBMI",
                                                "WHRadjBMI", 
                                                "T2DadjBMI",
                                                "WHRadjBMI", 
                                                "T2DadjBMI",
                                                "WHRadjBMI", 
                                                "T2DadjBMI",
                                                "WHRadjBMI", 
                                                "T2DadjBMI",
                                                "WHRadjBMI", 
                                                "T2DadjBMI",
                                                NA,
                                                NA)
                                      ),
                                 .Names = c("mean", "lower", "upper","group"),
                                 row.names = c(NA, 16L), 
                                 class = "data.frame")

tabletext <- cbind(c("", 
                     "Variant", 
                     representative_snps$SNP[1], 
                     "",
                     representative_snps$SNP[2],
                     "",
                     representative_snps$SNP[3],
                     "",
                     representative_snps$SNP[4],
                     "",
                     representative_snps$SNP[5],
                     "",
                     representative_snps$SNP[6],
                     "",
                     "",
                     ""),
                   c("", 
                     "GWAS", 
                     "WHRadjBMI", 
                     "T2DadjBMI",
                     "WHRadjBMI", 
                     "T2DadjBMI",
                     "WHRadjBMI", 
                     "T2DadjBMI",
                     "WHRadjBMI", 
                     "T2DadjBMI",
                     "WHRadjBMI", 
                     "T2DadjBMI",
                     "WHRadjBMI", 
                     "T2DadjBMI",
                     "",
                     ""),
                   c("",
                     "Effect Size",
                     representative_snps$WHRadjBMI_BETA[1],
                     representative_snps$T2DadjBMI_BETA[1],
                     representative_snps$WHRadjBMI_BETA[2],
                     representative_snps$T2DadjBMI_BETA[2],
                     representative_snps$WHRadjBMI_BETA[3],
                     representative_snps$T2DadjBMI_BETA[3],
                     representative_snps$WHRadjBMI_BETA[4],
                     representative_snps$T2DadjBMI_BETA[4],
                     representative_snps$WHRadjBMI_BETA[5],
                     representative_snps$T2DadjBMI_BETA[5],
                     representative_snps$WHRadjBMI_BETA[6],
                     representative_snps$T2DadjBMI_BETA[6],
                     "",
                     "")
)
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/credible_set/figures')

jpeg(filename = 'BETA_OR_discordant_snps.jpeg',
     units = 'in',
     res = 300,
     height = 5,
     width = 7)

cochrane_from_rmeta %>% 
  #group_by(group) %>%
  forestplot(clip = c(-0.25,0.05),
             # shapes_gp = forestplot::fpShapesGp()(box = c(box = "royalblue",
             #                                line = "darkblue") %>% lapply(function(x) gpar(fill = x, col = "#555555")),
             #                        default = gpar(vertices = TRUE)),
             ci.vertices = TRUE,
             ci.vertices.height = 0.05,
             labeltext = tabletext,
             boxsize = .1,
             line.margin = 0.1,
             xticks = c(-0.25, -0.15, 0, 0.05),
             col = fpColors(box = "royalblue",
                            line = "darkblue"),
             xlab = "Odds Ratio (T2DadjBMI) and Effect Size (WHRadjBMI)",
)

dev.off()
