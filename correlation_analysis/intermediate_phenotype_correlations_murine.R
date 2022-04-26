## load libraries
{
  library(vroom)
  library(caroline)
  library(pheatmap)
}

## load data
{
  setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/correlation_analysis/raw_data')
  murine_correlation_data <- readxl::read_excel(path = 'murine_data_correlations.xlsx')
  significant_murine_correlation_data <- subset(murine_correlation_data, pvalue < 0.05)
  murine_correlation_data <- murine_correlation_data[(order(murine_correlation_data$trait_name)),]
  murine_correlation_data <- murine_correlation_data[(order(murine_correlation_data$gene_symbol)),]
  genes_with_significant_correlations <- unique(significant_murine_correlation_data$gene_symbol)
  phenotypes_with_significant_correlations <- unique(significant_murine_correlation_data$trait_name)
  Zfp36l2 <- murine_correlation_data[(murine_correlation_data$gene_symbol == 'Zfp36l2'),]
  Zfp36l2 <- Zfp36l2[(order(Zfp36l2$pvalue)),]
  Zfp36l2 <- Zfp36l2[(!duplicated(Zfp36l2$trait_name)),]
  murine_correlation_data_heatmap <- rbind.data.frame(Zfp36l2,subset(murine_correlation_data, gene_symbol == 'Evi2a' | gene_symbol == 'Pam'))
  murine_correlation_data_heatmap <- murine_correlation_data_heatmap[(order(murine_correlation_data_heatmap$trait_name)),]
  murine_correlation_data_heatmap <- data.frame(Evi2a = subset(murine_correlation_data_heatmap, gene_symbol == 'Evi2a')$bicor,
                                               Zfp36l2 = subset(murine_correlation_data_heatmap, gene_symbol == 'Zfp36l2')$bicor,
                                               Pam = subset(murine_correlation_data_heatmap, gene_symbol == 'Pam')$bicor,
                                               row.names = unique(murine_correlation_data_heatmap$trait_name))
  murine_correlation_data_heatmap_transposed <- t(murine_correlation_data_heatmap)
  setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/correlation_analysis/figures')
  jpeg(filename = 'murine_correlations_heatmap_vertical.jpeg',
       height = 5,
       width = 7.5,
       units = 'in',
       res = 300)
  pheatmap(mat = murine_correlation_data_heatmap)
  dev.off()
  jpeg(filename = 'murine_correlations_heatmap_horizontal.jpeg',
       height = 5,
       width = 7.5,
       units = 'in',
       res = 300)
  pheatmap(mat = murine_correlation_data_heatmap_transposed)
  dev.off()
}