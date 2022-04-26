# load libraries ----------------------------------------------------------
  library(vroom)
  library(caroline)
  library(pheatmap)
  library(leafcutter)
  library(WGCNA)
  library(qvalue)
  library(stringr)
  library(corrplot)


# load and clean data -----------------------------------------------------
    setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/correlation_analysis/raw_data')
    METSIM_phenotypes <- vroom('METSIM_pheno.csv')
    METSIM_transcript_rawcounts <- vroom('METSIM_n426_adipose_Salmon_transcript_rawCounts.txt.gz')
    METSIM_gene_rawcounts <- vroom('METSIM_n426_adipose_gene_rawCounts.txt.gz')


    QTL_Genes_discordant_loci <- vroom('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_eqtl_colocalization/results/summary_figure.txt')
    QTL_Genes_discordant_loci <- unique(QTL_Genes_discordant_loci$GENE)
    
    QTL_Genes_discordant_loci <- sub("\\..*","",QTL_Genes_discordant_loci)
    METSIM_transcript_rawcounts$gene <- sub("\\..*","",METSIM_transcript_rawcounts$gene)
    METSIM_transcript_rawcounts <- subset(METSIM_transcript_rawcounts, gene %in% QTL_Genes_discordant_loci)
    METSIM_gene_rawcounts$geneID <- sub("\\..*","",METSIM_gene_rawcounts$geneID)
    METSIM_gene_rawcounts <- subset(METSIM_gene_rawcounts, geneID %in% QTL_Genes_discordant_loci)
    
    {
    # Keep Genes expressed in 50%+ of METSIM donors. Raw counts > 10 in 50%+ of donors.
    METSIM_gene_rawcounts_common <- METSIM_gene_rawcounts[,-c(1)] > 10
    row_sums_filtering <- rowSums(METSIM_gene_rawcounts_common)
    row_sums_filtering <- data.frame(geneID = METSIM_gene_rawcounts$geneID,
                                     METSIM_Donor_Expression = row_sums_filtering)
    row_sums_filtering <- subset(row_sums_filtering, METSIM_Donor_Expression > 426/2)
    METSIM_gene_rawcounts_common <- subset(METSIM_gene_rawcounts, 
                                           geneID %in% row_sums_filtering$geneID)
    rm(row_sums_filtering,METSIM_gene_rawcounts)
    
    # Keep Transcripts expressed in 50%+ of METSIM. Raw counts > 10 in 50%+ of donors.
    METSIM_transcript_rawcounts_common <- METSIM_transcript_rawcounts[,-c(1:4)] > 10
    row_sums_filtering <- rowSums(METSIM_transcript_rawcounts_common)
    row_sums_filtering <- data.frame(transcript = METSIM_transcript_rawcounts$transcript,
                                     METSIM_Donor_Expression = row_sums_filtering)
    row_sums_filtering <- subset(row_sums_filtering, METSIM_Donor_Expression > 426/2)
    METSIM_transcript_rawcounts_common <- subset(METSIM_transcript_rawcounts, 
                                                 transcript %in% row_sums_filtering$transcript)
    rm(row_sums_filtering,METSIM_transcript_rawcounts)
    }
    
    # Sort donors by ID
    {
      expression_donors <- as.numeric(colnames(METSIM_gene_rawcounts_common[-c(1)]))
      transcript_donors <- as.numeric(colnames(METSIM_transcript_rawcounts_common[-c(1:4)]))
      phenotype_donors <- METSIM_phenotypes$METSIM_ID
      expression_donors <- expression_donors[expression_donors %in% phenotype_donors & expression_donors %in% transcript_donors]
      transcript_donors <- transcript_donors[transcript_donors %in% phenotype_donors & transcript_donors %in% expression_donors]
      phenotype_donors <- phenotype_donors[phenotype_donors %in% expression_donors]
      expression_donors <- expression_donors[order(expression_donors)]
      transcript_donors <- transcript_donors[order(transcript_donors)]
      phenotype_donors <- phenotype_donors[order(phenotype_donors)]
      expression_donors <- c('geneID',expression_donors)
      transcript_donors <- c('transcript','gene',transcript_donors)
      METSIM_gene_rawcounts_common <- METSIM_gene_rawcounts_common[expression_donors]
      METSIM_transcript_rawcounts_common <- METSIM_transcript_rawcounts_common[transcript_donors]
      METSIM_phenotypes <- subset(METSIM_phenotypes, METSIM_ID %in% phenotype_donors)
      METSIM_phenotypes <- METSIM_phenotypes[(order(METSIM_phenotypes$METSIM_ID)),]
    }
    
    rm(expression_donors,phenotype_donors,transcript_donors)
    setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/correlation_analysis/formatted_data')
    geneIDs <- METSIM_gene_rawcounts_common$geneID
    METSIM_gene_rawcounts_common <- t(METSIM_gene_rawcounts_common)
    colnames(METSIM_gene_rawcounts_common) = METSIM_gene_rawcounts_common[1,]
    METSIM_gene_rawcounts_common <- METSIM_gene_rawcounts_common[-1,]
    write.delim(METSIM_gene_rawcounts_common, "METSIM_gene_rawcounts_common.txt")
    write.delim(METSIM_transcript_rawcounts_common, "METSIM_transcript_rawcounts_common.txt")
    write.delim(METSIM_phenotypes, "METSIM_phenotypes.txt")
    

    
# correlation analysis ----------------------------------------------------
    
expression_correlation_phenotype <- function(transcriptomic_dataframe, phenotypic_dataframe) {
      correlation_data <- data.frame(phenotype = NA, gene = NA, correlation_coefficient = NA, 
                                     correlation_pvalue = NA, correlation_tstatistic = NA)
      phenotypic_dataframe <- phenotypic_dataframe[-c(1:2)]
      phenotype_list <- colnames(phenotypic_dataframe)
      gene_list <- colnames(transcriptomic_dataframe)
      for (phenotype in phenotype_list) {
        phenotype_data <- as.numeric(phenotypic_dataframe[[phenotype]])
        for (gene in gene_list) {
          gene_transcription_data <- as.numeric(transcriptomic_dataframe[[gene]])
          correlation_test <- bicorAndPvalue(x = gene_transcription_data, y = phenotype_data)
          correlation_coefficient <- correlation_test[["bicor"]][1]
          correlation_p <- correlation_test[["p"]][1]
          correlation_t <- correlation_test[["t"]][1]
          correlation_values <- c(phenotype, gene, correlation_coefficient,
                                  correlation_p,correlation_t)
          correlation_data <- rbind.data.frame(correlation_values,correlation_data)
        }
      }
      correlation_data <- na.omit(correlation_data)
      return(correlation_data)
    }
    
expression_correlation_phenotype_transcripts <- function(transcriptomic_dataframe, phenotypic_dataframe) {
      correlation_data <- data.frame(phenotype = NA, gene = NA, transcript = NA, correlation_coefficient = NA, 
                                     correlation_pvalue = NA, correlation_tstatistic = NA)
      phenotypic_dataframe <- phenotypic_dataframe[-c(1:2)]
      phenotype_list <- colnames(phenotypic_dataframe)
      transcript_list <- transcriptomic_dataframe$transcript
      gene_list <- transcriptomic_dataframe$gene
      transcriptomic_dataframe <- transcriptomic_dataframe[-c(1:2)]
      for (phenotype_number in 1:length(phenotype_list)) {
        phenotype_name <- phenotype_list[phenotype_number]
        phenotype_data <- unlist(x = phenotypic_dataframe[phenotype_number], 
                                 use.names = FALSE)
        phenotype_data <- as.numeric(phenotype_data)
        for (transcript_number in 1:length(transcript_list)) {
          gene_name <- gene_list[transcript_number]
          transcript_name <- transcript_list[transcript_number]
          gene_transcription_data <- unlist(x = transcriptomic_dataframe[transcript_number,], 
                                            use.names = FALSE)
          gene_transcription_data <- as.numeric(gene_transcription_data)
          
          correlation_test <- bicorAndPvalue(x = gene_transcription_data, y = phenotype_data)
          correlation_coefficient <- correlation_test[["bicor"]][1]
          correlation_p <- correlation_test[["p"]][1]
          correlation_t <- correlation_test[["t"]][1]
          correlation_values <- c(phenotype_name, gene_name, transcript_name, correlation_coefficient,
                                  correlation_p,correlation_t)
          correlation_data <- rbind.data.frame(correlation_values,correlation_data)
        }
      }
      correlation_data <- na.omit(correlation_data)
      return(correlation_data)
    }
    
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/correlation_analysis/formatted_data')
METSIM_gene_rawcounts_common <- read.delim('METSIM_gene_rawcounts_common.txt')
METSIM_phenotypes <- read.delim('METSIM_phenotypes.txt')
correlations_expression <- expression_correlation_phenotype(transcriptomic_dataframe = METSIM_gene_rawcounts_common,
                                                            phenotypic_dataframe = METSIM_phenotypes)
correlations_transcript <- expression_correlation_phenotype_transcripts(transcriptomic_dataframe = METSIM_transcript_rawcounts_common,
                                                                            phenotypic_dataframe = METSIM_phenotypes)
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/correlation_analysis/results')
write.delim(correlations_expression, 'human_intermediate_phenotype_correlations_expression.txt')
write.delim(correlations_transcript, 'human_intermediate_phenotype_correlations_trascript.txt')
  
  


# plot results ------------------------------------------------------------

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/correlation_analysis/results')
correlations_expression <- vroom('human_intermediate_phenotype_correlations_expression.txt')
correlations_transcript <- vroom('human_intermediate_phenotype_correlations_trascript.txt')
qvalues_expression <- qvalue(p = correlations_expression$correlation_pvalue, fdr.level = 0.05)
qvalues_transcript <- qvalue(p = correlations_transcript$correlation_pvalue, fdr.level = 0.05)
correlations_expression$qvalues <- qvalues_expression$qvalues
correlations_transcript$qvalues <- qvalues_transcript$qvalues

#expression heatmap
correlations_expression_heatmap <- correlations_expression
correlations_expression_heatmap$phenotype <- str_replace_all(correlations_expression_heatmap$phenotype,"_"," ")
correlations_expression_heatmap_significance <- data.frame('THADA-AS' = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000234936')$qvalues),
                                                           #PEPD = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000124299')$qvalues),
                                                           GIN1 = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000145723')$qvalues),
                                                           PPIP5K2 = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000145725')$qvalues),
                                                           PAM = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000145730')$qvalues),
                                                           #SLC38A9 = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000177058')$qvalues),
                                                           #CHST8 = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000124302')$qvalues),
                                                           row.names = unique(correlations_expression_heatmap$phenotype))
correlations_expression_heatmap <- data.frame('THADA-AS' = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000234936')$correlation_coefficient),
                                              #PEPD = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000124299')$correlation_coefficient),
                                              GIN1 = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000145723')$correlation_coefficient),
                                              PPIP5K2 = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000145725')$correlation_coefficient),
                                              PAM = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000145730')$correlation_coefficient),
                                              #SLC38A9 = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000177058')$correlation_coefficient),
                                              #CHST8 = as.numeric(subset(correlations_expression_heatmap, gene == 'ENSG00000124302')$correlation_coefficient),
                                              row.names = unique(correlations_expression_heatmap$phenotype))
correlations_expression_heatmap <- as.matrix.data.frame(correlations_expression_heatmap)
correlations_expression_heatmap_significance <- as.matrix.data.frame(correlations_expression_heatmap_significance)
colnames(correlations_expression_heatmap) <- c('THADA-AS','GIN1','PPIP5K2','PAM')
colnames(correlations_expression_heatmap_significance) <- c('THADA-AS','GIN1','PPIP5K2','PAM')
correlations_expression_heatmap <- t(correlations_expression_heatmap)
correlations_expression_heatmap_significance <- t(correlations_expression_heatmap_significance)
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/correlation_analysis/figures')
jpeg(filename = 'human_correlations_expression_heatmap.jpeg',
          height = 5,
          width = 10,
          units = 'in',
          res = 300)
corrplot(corr = correlations_expression_heatmap, 
         p.mat = correlations_expression_heatmap_significance,
         method = 'square',
         sig.level = c(0.001, 0.01, 0.05),
         insig = 'label_sig',
         tl.col = 'black',
         col.lim = c(-0.30,0.30),
         is.corr = F,
         pch.cex = 0.9,
         pch.col = 'black')
dev.off()
    
# PAM transcripts heatmap
PAM_ensemblID <- 'ENSG00000145730'
correlations_transcript_heatmap <- subset(correlations_transcript, gene == PAM_ensemblID)
correlations_transcript_heatmap$phenotype <- str_replace_all(correlations_transcript_heatmap$phenotype,"_"," ")
correlations_transcript_heatmap_significance <- data.frame(ENST00000504691.1 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000504691.1')$qvalues),
                                                           ENST00000504456.1 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000504456.1')$qvalues),
                                                           ENST00000510006.1 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000510006.1')$qvalues),
                                                           ENST00000379799.3 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000379799.3')$qvalues),
                                                           ENST00000455264.2 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000455264.2')$qvalues),
                                                           ENST00000348126.2 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000348126.2')$qvalues),
                                                           ENST00000346918.2 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000346918.2')$qvalues),
                                                           ENST00000438793.3 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000438793.3')$qvalues),
                                                           row.names = unique(correlations_transcript_heatmap$phenotype))
correlations_transcript_heatmap <- data.frame(ENST00000504691.1 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000504691.1')$correlation_coefficient),
                                              ENST00000504456.1 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000504456.1')$correlation_coefficient),
                                              ENST00000510006.1 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000510006.1')$correlation_coefficient),
                                              ENST00000379799.3 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000379799.3')$correlation_coefficient),
                                              ENST00000455264.2 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000455264.2')$correlation_coefficient),
                                              ENST00000348126.2 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000348126.2')$correlation_coefficient),
                                              ENST00000346918.2 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000346918.2')$correlation_coefficient),
                                              ENST00000438793.3 = as.numeric(subset(correlations_transcript_heatmap, transcript == 'ENST00000438793.3')$correlation_coefficient),
                                              row.names = unique(correlations_transcript_heatmap$phenotype))
correlations_transcript_heatmap <- as.matrix.data.frame(correlations_transcript_heatmap)
correlations_transcript_heatmap_significance <- as.matrix.data.frame(correlations_transcript_heatmap_significance)
correlations_transcript_heatmap <- t(correlations_transcript_heatmap)
correlations_transcript_heatmap_significance <- t(correlations_transcript_heatmap_significance)
    
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/correlation_analysis/figures')
jpeg(filename = 'human_correlations_heatmap_PAM_transcripts.jpeg',
  height = 5,
  width = 10,
  units = 'in',
  res = 300)
corrplot(corr = correlations_transcript_heatmap, 
         p.mat = correlations_transcript_heatmap_significance,
         method = 'square',
         sig.level = c(0.001, 0.01, 0.05),
         insig = 'label_sig',
         tl.col = 'black',
         col.lim = c(-0.30,0.30),
         is.corr = F,
         pch.cex = 0.9,
         pch.col = 'black')
dev.off()
  

# clean memory ------------------------------------------------------------

gc()
    