
# load libraries ----------------------------------------------------------
library(ggplot2)
library(caroline)
# declare functions ---------------------------------------------------------------


# load data ---------------------------------------------------------------
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/results')
summary_figure_old <- read.delim('summary_figure.txt')
GENES <- summary_figure$GENE
write.delim(GENES, 'GENES.txt')
GENES <- read.delim('GENES.txt')
summary_figure$GENE <- GENES$SYMBOL

# clean data --------------------------------------------------------------


# implement functions -----------------------------------------------------


