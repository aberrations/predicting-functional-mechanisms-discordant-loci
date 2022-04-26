##library loading
library(vroom)
library(biomaRt)
library(caroline)

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical_no_brain/raw_data/expression_specificity')
median_expression_counts <- vroom('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')

#for tissue groups, will find expression score with highest score and save it
Adipose <- data.frame(median_expression_counts$Name,
                             median_expression_counts$Description,
                             combined_counts = pmax(median_expression_counts$`Adipose - Subcutaneous`,median_expression_counts$`Adipose - Visceral (Omentum)`))

cleaned_expression_counts <- data.frame(median_expression_counts$Name,
                                        median_expression_counts$Description,
                                        Adipose$combined_counts,
                                        median_expression_counts$Pancreas,
                                        median_expression_counts$Liver,
                                        median_expression_counts$`Muscle - Skeletal`)
rm(Adipose)
colnames(cleaned_expression_counts) <- c('FEATURE_ID','FEATURE_NAME','Adipose','Islets','Liver','SkeletalMuscle')
cleaned_expression_counts$totalexpression <- cleaned_expression_counts$Adipose+cleaned_expression_counts$Islets+cleaned_expression_counts$Liver+cleaned_expression_counts$SkeletalMuscle
cleaned_expression_counts <- subset(cleaned_expression_counts, totalexpression != 0 )

cleaned_expression_counts$Adipose <- cleaned_expression_counts$Adipose/cleaned_expression_counts$totalexpression
cleaned_expression_counts$Islets <- cleaned_expression_counts$Islets/cleaned_expression_counts$totalexpression
cleaned_expression_counts$Liver <- cleaned_expression_counts$Liver/cleaned_expression_counts$totalexpression
cleaned_expression_counts$SkeletalMuscle <- cleaned_expression_counts$SkeletalMuscle/cleaned_expression_counts$totalexpression
cleaned_expression_counts$totalexpression <- cleaned_expression_counts$totalexpression/cleaned_expression_counts$totalexpression
cleaned_expression_counts <- cleaned_expression_counts[(!duplicated(cleaned_expression_counts$FEATURE_NAME)),]
#TACTICAL_expression_counts <- vroom('/Library/Frameworks/R.framework/Versions/4.1/Resources/library/tactical_no_brain/extdata/gene-expression-specificity-scores.txt')
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh=37)
attributes = c('ensembl_gene_id_version','external_gene_name','chromosome_name','start_position','end_position')
gene_ids <- unique(cleaned_expression_counts$FEATURE_NAME)
filters = c('external_gene_name')
biomart_query <- getBM(attributes = attributes,
                       filters = filters,
                       values = gene_ids,
                       mart = ensembl)
biomart_query$chromosome_name <- as.numeric(biomart_query$chromosome_name)
biomart_query <- na.omit(biomart_query)
biomart_query <-biomart_query[(!duplicated(biomart_query$external_gene_name)),]
expression_specificity_scores <- cleaned_expression_counts
expression_specificity_scores <- subset(expression_specificity_scores, FEATURE_NAME %in% biomart_query$external_gene_name)
biomart_query <- subset(biomart_query, external_gene_name %in% expression_specificity_scores$FEATURE_NAME)
expression_specificity_scores <- expression_specificity_scores[(order(expression_specificity_scores$FEATURE_NAME)),]
biomart_query <- biomart_query[(order(biomart_query$external_gene_name)),]

expression_specificity_scores <- data.frame(
  expression_specificity_scores$FEATURE_ID,
  expression_specificity_scores$FEATURE_NAME,
  biomart_query$chromosome_name,
  biomart_query$start_position,
  biomart_query$end_position,
  expression_specificity_scores$Adipose,
  expression_specificity_scores$Liver,
  expression_specificity_scores$Islets,
  expression_specificity_scores$SkeletalMuscle)
colnames(expression_specificity_scores) <- c('FEATURE_ID','FEATURE_NAME','CHR','START','END','Adipose','Liver','Islets','SkeletalMuscle')
expression_specificity_scores$CHR <- as.numeric(expression_specificity_scores$CHR)
expression_specificity_scores <- na.omit(expression_specificity_scores)
expression_specificity_scores <- expression_specificity_scores[(order(expression_specificity_scores$START)),]
expression_specificity_scores <- expression_specificity_scores[(order(expression_specificity_scores$CHR)),]
expression_specificity_scores$CHR <- paste('chr',expression_specificity_scores$CHR,sep = '')
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical_no_brain/formatted_data')
write.delim(expression_specificity_scores, 'gene-expression-specificity-scores.txt')



