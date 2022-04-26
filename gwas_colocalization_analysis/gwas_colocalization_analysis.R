## load libraries
{
  library(coloc)
  library(caroline)
  library(dplyr)
  library(locuscomparer)
  library(ggplot2)
  library(vroom)
  library(hyprcoloc)
}

## load data
{
  setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/raw_data')
  WHRadjBMI <- vroom('whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt')
  BMI <- vroom('bmi.giant-ukbb.meta-analysis.combined.23May2018.txt')
  T2DadjBMI <- vroom('Mahajan.NatGenet2018b.T2Dbmiadj.European.txt', delim = '\t')
}

## declare functions
{
  clean_pulit_gwas <- function(pulit_gwas) {
    reformatted_pulit_gwas <- data.frame(SNP = sub(":.*","",pulit_gwas$SNP),
                                         CHR = pulit_gwas$CHR,
                                         POS = pulit_gwas$POS,
                                         EA = pulit_gwas$Tested_Allele,
                                         NEA = pulit_gwas$Other_Allele,
                                         EAF = pulit_gwas$Freq_Tested_Allele,
                                         BETA = pulit_gwas$BETA,
                                         VARBETA = (pulit_gwas$SE)^2,
                                         P = pulit_gwas$P,
                                         N = pulit_gwas$N,
                                         CHR_POS = paste(pulit_gwas$CHR,':',pulit_gwas$POS,'_',
                                                         pulit_gwas$Other_Allele,'/',
                                                         pulit_gwas$Tested_Allele,sep = '')
                                         )
    return(reformatted_pulit_gwas)
  }
  
  clean_mahajan_gwas <- function(mahajan_gwas) {
    reformatted_mahajan_gwas <- data.frame(CHR = mahajan_gwas$Chr,
                                         POS = mahajan_gwas$Pos,
                                         EA = mahajan_gwas$EA,
                                         NEA = mahajan_gwas$NEA,
                                         EAF = mahajan_gwas$EAF,
                                         BETA = mahajan_gwas$Beta,
                                         VARBETA = (mahajan_gwas$SE)^2,
                                         P = mahajan_gwas$Pvalue,
                                         N = mahajan_gwas$Neff,
                                         CHR_POS = paste(mahajan_gwas$Chr,':',mahajan_gwas$Pos,'_',
                                                         mahajan_gwas$NEA,'/',
                                                         mahajan_gwas$EA,sep = '')
    )
    return(reformatted_mahajan_gwas)
  }
  
  coloc_analysis <- function(bfd_gwas, bfd_gwas_name, t2d_gwas, t2d_gwas_name){
    setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/results')
    
    dir.create("coloc_abf", showWarnings = FALSE)
    dir.create("coloc_abf/bayesian_factors", showWarnings = FALSE)
    dir.create("hyprcoloc", showWarnings = FALSE)

    summary_statistics_coloc_abf <- data.frame(CHR = NA,
                                     START = NA,
                                     STOP = NA,
                                     LEAD_SNP_BF = NA,
                                     LEAD_SNP_P_COMBINED = NA,
                                     nsnps = NA,
                                     PPH0 = NA,
                                     PPH1 = NA,
                                     PPH2 = NA,
                                     PPH3 = NA,
                                     PPH4 = NA)
    summary_statistics_hyprcoloc <- data.frame(CHR = NA,
                                     START = NA,
                                     STOP = NA,
                                     LEAD_SNP_P_COMBINED = NA,
                                     CANDIDATE_SNP = NA,
                                     POSTERIOR_PROB = NA,
                                     REGIONAL_PROB = NA,
                                     POSTERIOR_EXPLAINED_SNP = NA,
                                     DROPPED_TRAIT = NA)
    chromosomes <- c(1:22)
    for (chromosome in chromosomes) {
      chromosome_t2d_gwas <- subset(t2d_gwas, CHR == chromosome)
      chromosome_bfd_gwas <- subset(bfd_gwas, CHR == chromosome)
      chromosome_sig_gwas <- subset(chromosome_bfd_gwas, P < 5e-08)
      chromosome_sig_gwas_2 <- subset(chromosome_t2d_gwas, P < 5e-08)
      chromosome_sig_gwas_2 <- subset(chromosome_bfd_gwas, SNP %in% chromosome_sig_gwas_2$SNP)
      chromosome_sig_gwas <- rbind.data.frame(chromosome_sig_gwas,chromosome_sig_gwas_2)
      chromosome_sig_gwas <- chromosome_sig_gwas[(!duplicated(chromosome_sig_gwas$SNP)),] #unique df including all genome-wide significant SNPS for both GWAS.
      chromosomal_end <- max(chromosome_bfd_gwas$POS)
      window_length = 200000 #bases
      stop_position = 0
      while (stop_position < chromosomal_end & !is.na(stop_position) & nrow(chromosome_sig_gwas) > 0) {
        center_position = min(chromosome_sig_gwas$POS)
        start_position = center_position - window_length/2
        stop_position = center_position + window_length/2
        chromosome_sig_gwas <- subset(chromosome_sig_gwas, POS >= stop_position)
        
        while (nrow(chromosome_sig_gwas) > 0 & (stop_position + window_length/2 > min(chromosome_sig_gwas$POS))) 
        {
          stop_position = stop_position + window_length/2
          chromosome_sig_gwas <- subset(chromosome_sig_gwas, POS > stop_position)
        }
        
        t2d_locus <- subset(chromosome_t2d_gwas, POS > start_position & POS < stop_position)
        whr_locus <- subset(chromosome_bfd_gwas, POS > start_position & POS < stop_position)
        t2d_locus = t2d_locus[(t2d_locus$SNP %in% whr_locus$SNP),]
        whr_locus = whr_locus[(whr_locus$SNP %in% t2d_locus$SNP),]
        t2d_locus <- t2d_locus[(!duplicated(t2d_locus$SNP)),]
        whr_locus <- whr_locus[(!duplicated(whr_locus$SNP)),]
        lead_snp_p = whr_locus
        lead_snp_p$P_combined <- whr_locus$P + t2d_locus$P
        lead_snp_p <- lead_snp_p[(order(lead_snp_p$P_combined)),]
        lead_snp_p <- lead_snp_p$SNP[1]
        start_position = min(t2d_locus$POS)
        stop_position = max(t2d_locus$POS)
        if (nrow(t2d_locus) > 1 & nrow(t2d_locus) == nrow(whr_locus)) {
          
          ## coloc_abf
          datasetA = list(beta = t2d_locus$BETA, snp = t2d_locus$SNP, N = t2d_locus$N, MAF = t2d_locus$EAF, pvalues = t2d_locus$P, varbeta = t2d_locus$VARBETA, type = "cc", s = .09)
          datasetB = list(beta = whr_locus$BETA, MAF = whr_locus$EAF, snp = whr_locus$SNP, N = whr_locus$N, pvalues = whr_locus$P, varbeta = whr_locus$VARBETA, type = "quant")
          coloc_results_filename <- paste("coloc_abf/bayesian_factors/chr",chromosome,"_",start_position,"_",stop_position,'_',bfd_gwas_name,'_',t2d_gwas_name,"_coloc_abf.txt", sep = "")
          coloc_results <- coloc.abf(datasetA, datasetB, p12 = 5e-06)
          coloc_abf_results <- coloc_results$results
          write.delim(coloc_abf_results, coloc_results_filename)
          lead_snp_abf <- coloc_abf_results[(order(coloc_abf_results$SNP.PP.H4)),]
          lead_snp_abf <- lead_snp_abf$snp[nrow(lead_snp_abf)]
          locus_summary <- t(data.frame(coloc_results$summary))
          locus_summary <- c(chromosome,start_position,stop_position,lead_snp_abf,lead_snp_p,locus_summary)
          summary_statistics_coloc_abf <- rbind(summary_statistics_coloc_abf,locus_summary)
          number_of_snps <- nrow(t2d_locus)
          ##hyprcoloc
          if(number_of_snps < 1.5e4) {
            trait_names <- c(bfd_gwas_name,t2d_gwas_name)
            rsid <- whr_locus$SNP
            BETAs <- data.frame(WHRadjBMI = whr_locus$BETA, 
                                T2DadjBMI = t2d_locus$BETA)
            rownames(BETAs) <- rsid
            BETAs <- as.matrix.data.frame(BETAs)
            SEs <- data.frame(WHRadjBMI = (whr_locus$VARBETA)^0.5, T2DadjBMI = (t2d_locus$VARBETA)^.5)
            rownames(SEs) <- rsid
            SEs <- as.matrix.data.frame(SEs)
            hyprcoloc_results <- hyprcoloc(effect.est = BETAs,
                                           effect.se = SEs,
                                           prior.12 = 5e-06,
                                           snp.id = rsid)
            hyprcoloc_results <- hyprcoloc_results$results
            for (row_number in 1:nrow(hyprcoloc_results)) {
              candidate_snp <- hyprcoloc_results$candidate_snp[row_number]
              posterior_prob <- hyprcoloc_results$posterior_prob[row_number]
              regional_prob <- hyprcoloc_results$regional_prob[row_number]
              posterior_explained_by_snp <- hyprcoloc_results$posterior_explained_by_snp[row_number]
              dropped_trait <- hyprcoloc_results$dropped_trait[row_number]
              locus_summary <- c(chromosome,start_position,stop_position,lead_snp_p,candidate_snp,posterior_prob,regional_prob,posterior_explained_by_snp,dropped_trait)
              summary_statistics_hyprcoloc <- rbind.data.frame(summary_statistics_hyprcoloc,locus_summary)
            }
          }
        }
      }
    }
    summary_statistics_coloc_abf <- na.omit(summary_statistics_coloc_abf)
    summary_statistics_hyprcoloc[1,] <- NULL
    write.delim(summary_statistics_coloc_abf,"coloc_abf/summary_statistics_coloc_abf.txt")
    write.delim(summary_statistics_hyprcoloc,"hyprcoloc/summary_statistics_hyprcoloc.txt")
    print("Whole genome is complete.")
  }
  
  locuscompare <- function(in_fn1, in_fn2, marker_col1 = "rsid", pval_col1 = "pval",
                          title1 = "eQTL",marker_col2 = "rsid", pval_col2 = "pval", title2 = "GWAS",
                          snp = NULL, population = "EUR", combine = TRUE, legend = TRUE,
                          legend_position = c('bottomright','topright','topleft'),
                          lz_ylab_linebreak = FALSE, genome = c('hg19','hg38')) {
    d1 = read_metal(in_fn1, marker_col1, pval_col1)
    d2 = read_metal(in_fn2, marker_col2, pval_col2)
    
    merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
    genome = match.arg(genome)
    merged = get_position(merged, genome)
    
    chr = unique(merged$chr)
    if (length(chr) != 1){
      mode_chromosome <- Mode(merged$chr)
      merged <- subset(merged, chr == mode_chromosome)
    }
    
    snp = get_lead_snp(merged, snp)
    ld = retrieve_LD(chr, snp, population)
    p = make_combined_plot(merged, title1, title2, ld, chr, snp, combine,
                           legend, legend_position, lz_ylab_linebreak)
    return(p)
  }
  
  locuscompare_t2d_bfd <- function(colocalization_results, bfd_gwas, bfd_gwas_name, t2d_gwas, t2d_gwas_name) {
    setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/results')
    dir.create("locuscompare", showWarnings = FALSE)
    setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/results/locuscompare')
    t2d_gwas <- t2d_gwas[(!duplicated(t2d_gwas$SNP)),]
    bfd_gwas <- bfd_gwas[(!duplicated(bfd_gwas$SNP)),]
    for (locus in 1:nrow(colocalization_results)) {
      chromosome <- colocalization_results$CHR[locus]
      start_position <- colocalization_results$START[locus]
      stop_position <- colocalization_results$STOP[locus]
      locus_bfd_gwas <- subset(bfd_gwas, CHR == chromosome & POS > start_position & POS < stop_position)
      locus_t2d_gwas <- subset(t2d_gwas, CHR == chromosome & POS > start_position & POS < stop_position)
      locus_bfd_gwas <- data.frame(rsid = locus_bfd_gwas$SNP, pval = locus_bfd_gwas$P)
      locus_t2d_gwas <- data.frame(rsid = locus_t2d_gwas$SNP, pval = locus_t2d_gwas$P)
      locuscompare(in_fn1 = locus_bfd_gwas,
                   in_fn2 = locus_t2d_gwas,
                   title1 = bfd_gwas_name,
                   title2 = t2d_gwas_name,
                   population = "EUR",
                   genome = 'hg19')
      ggsave(filename = paste(chromosome,'_',start_position,'_',stop_position,'_',bfd_gwas_name,'_',t2d_gwas_name,'.jpeg',sep = ''), 
             width = 7.5, height = 5)
    }
  }
  
  get_effect_directions <- function(snp_list, t2d_gwas, bfd_gwas, bmi_gwas) {
    t2d_gwas <- subset(t2d_gwas, SNP %in% snp_list)
    bfd_gwas <- subset(bfd_gwas, SNP %in% snp_list)
    bmi_gwas <- subset(bmi_gwas, SNP %in% snp_list)
    
    bfd_gwas$BETA[bfd_gwas$P > 5e-03] <- 0
    t2d_gwas$BETA[t2d_gwas$P > 5e-03] <- 0
    bmi_gwas$BETA[bmi_gwas$P > 5e-05] <- 0  
    
    ordered_heatmap <- data.frame(T2DadjBMI = sign(t2d_gwas$BETA),
                                  WHRadjBMI = sign(bfd_gwas$BETA),
                                  BMI = sign(bmi_gwas$BETA))
    rownames(ordered_heatmap) <- t2d_gwas$SNP
    
    protective_SNPs <- ordered_heatmap$T2DadjBMI < 0
    ordered_heatmap$T2DadjBMI[protective_SNPs] <- -1*ordered_heatmap$T2DadjBMI[protective_SNPs]
    ordered_heatmap$WHRadjBMI[protective_SNPs] <- -1*ordered_heatmap$WHRadjBMI[protective_SNPs]
    ordered_heatmap$BMI[protective_SNPs] <- -1*ordered_heatmap$BMI[protective_SNPs]
    ordered_heatmap <- ordered_heatmap[(order(-ordered_heatmap$BMI)),]
    ordered_heatmap <- ordered_heatmap[(order(-ordered_heatmap$WHRadjBMI)),]
    return(ordered_heatmap)
  }
  
}

## clean data
{
  T2DadjBMI <- clean_mahajan_gwas(T2DadjBMI)
  WHRadjBMI <- clean_pulit_gwas(WHRadjBMI)
  BMI <- clean_pulit_gwas(pulit_gwas = BMI)
  
  rm(clean_pulit_gwas,clean_mahajan_gwas)
  T2DadjBMI <- T2DadjBMI[(!duplicated(T2DadjBMI$CHR_POS)),]
  WHRadjBMI <- WHRadjBMI[(!duplicated(WHRadjBMI$CHR_POS)),]
  BMI <- BMI[(!duplicated(BMI$CHR_POS)),]
  BMI <- subset(BMI, EAF != 0 & EAF != 1)
  T2DadjBMI <- subset(T2DadjBMI, EAF != 0 & EAF != 1)
  WHRadjBMI <- subset(WHRadjBMI, EAF != 0 & EAF != 1)
  T2DadjBMI <- subset(T2DadjBMI, CHR_POS %in% WHRadjBMI$CHR_POS)
  WHRadjBMI <- subset(WHRadjBMI, CHR_POS %in% T2DadjBMI$CHR_POS)
  BMI <- subset(BMI, CHR_POS %in% T2DadjBMI$CHR_POS & CHR_POS %in% WHRadjBMI$CHR_POS)
  BMI <- BMI[(order(BMI$CHR_POS)),]
  WHRadjBMI <- WHRadjBMI[(order(WHRadjBMI$CHR_POS)),]
  T2DadjBMI <- T2DadjBMI[(order(T2DadjBMI$CHR_POS)),]
  WHRadjBMI <- na.omit(WHRadjBMI)
  T2DadjBMI <- na.omit(T2DadjBMI)
  BMI <- na.omit(BMI)
  T2DadjBMI$SNP <- WHRadjBMI$SNP
  setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/formatted_data')
  write.delim(T2DadjBMI,'T2DadjBMI_GWAS.txt')
  write.delim(WHRadjBMI,'WHRadjBMI_GWAS.txt')
  write.delim(BMI,'BMI_GWAS.txt')
}

## colocalization analysis
{
  T2DadjBMI <- vroom('T2DadjBMI_GWAS.txt')
  WHRadjBMI <- vroom('WHRadjBMI_GWAS.txt')
  BMI <- vroom('BMI_GWAS.txt')
  coloc_analysis(bfd_gwas = WHRadjBMI, t2d_gwas = T2DadjBMI, bfd_gwas_name = "WHRadjBMI", t2d_gwas_name = "T2DadjBMI")
  summary_statistics_coloc_abf <- vroom('coloc_abf/summary_statistics_coloc_abf.txt')
  summary_statistics_coloc_abf <- subset(summary_statistics_coloc_abf, PPH4 + PPH3 > 0.7)
  summary_statistics_hyprcoloc <- vroom('hyprcoloc/summary_statistics_hyprcoloc.txt')
  summary_statistics_hyprcoloc <- subset(summary_statistics_hyprcoloc, REGIONAL_PROB > 0.70)
  significant_locus_coordinates <- data.frame(CHR = c(summary_statistics_coloc_abf$CHR,summary_statistics_hyprcoloc$CHR),
                                              START = c(summary_statistics_coloc_abf$START,summary_statistics_hyprcoloc$START),
                                              STOP = c(summary_statistics_coloc_abf$STOP,summary_statistics_hyprcoloc$STOP))
  significant_locus_coordinates$COORDINATES <- paste(significant_locus_coordinates$CHR,'_', significant_locus_coordinates$START,'_',
                                                     significant_locus_coordinates$STOP,sep = '')

  significant_locus_coordinates <- subset(significant_locus_coordinates, !duplicated(COORDINATES))
  write.delim(df = significant_locus_coordinates, file = "significant_locus_coordinates.txt")
  locuscompare_t2d_bfd(colocalization_results = significant_locus_coordinates, bfd_gwas = WHRadjBMI,
                       t2d_gwas = T2DadjBMI, bfd_gwas_name = "WHRadjBMI", t2d_gwas_name = "T2DadjBMI")
  
  summary_statistics_coloc_abf$COORDINATES <- paste(summary_statistics_coloc_abf$CHR,'_', summary_statistics_coloc_abf$START,'_',
                                                    summary_statistics_coloc_abf$STOP,sep = '')
  summary_statistics_hyprcoloc$COORDINATES <- paste(summary_statistics_hyprcoloc$CHR,'_', summary_statistics_hyprcoloc$START,'_',
                                                    summary_statistics_hyprcoloc$STOP,sep = '')
  
  visual_inspection_results <- vroom('locuscompare/visual_inspection.txt')
  visual_inspection_results <- subset(visual_inspection_results, VISUAL_INSPECTION == 'PASS')
  setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/figures')
  
  venn.diagram(x= list('COLOC_ABF' = summary_statistics_coloc_abf$COORDINATES,
                       HYPRCOLOC = summary_statistics_hyprcoloc$COORDINATES,
                       VISUAL_INSPECTION = visual_inspection_results$COORDINATES),
               filename = "venn_diagram_colocalization_analysis.jpeg")
  counts <- data.frame(COORDINATES = significant_locus_coordinates$COORDINATES,
                       COUNTS_COLOC = significant_locus_coordinates$COORDINATES %in% summary_statistics_coloc_abf$COORDINATES,
                       COUNTS_HYPRCOLOC = significant_locus_coordinates$COORDINATES %in% summary_statistics_hyprcoloc$COORDINATES,
                       COUNTS_VISUAL_INSPECTION = significant_locus_coordinates$COORDINATES %in% visual_inspection_results$COORDINATES)
  counts$COUNTS_TOTAL <- counts$COUNTS_COLOC+counts$COUNTS_HYPRCOLOC+counts$COUNTS_VISUAL_INSPECTION
  counts <- subset(counts, COUNTS_TOTAL > 1)
  significant_locus_coordinates <- subset(significant_locus_coordinates, COORDINATES %in% counts$COORDINATES)
  write.delim(significant_locus_coordinates, "significant_loci_overlap.txt")
}

## find associated effects at lead variants
{
  setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_colocalization_analysis/results/')
  significant_locus_coordinates <- vroom('significant_loci_overlap.txt')
  summary_statistics_coloc_abf <- vroom('coloc_abf/summary_statistics_coloc_abf.txt')
  summary_statistics_hyprcoloc <- vroom('hyprcoloc/summary_statistics_hyprcoloc.txt')
  visual_inspection_results <- vroom('locuscompare/visual_inspection.txt')
  summary_statistics_coloc_abf$COORDINATES <- paste(summary_statistics_coloc_abf$CHR,'_', summary_statistics_coloc_abf$START,'_',
                                                    summary_statistics_coloc_abf$STOP,sep = '')
  summary_statistics_hyprcoloc$COORDINATES <- paste(summary_statistics_hyprcoloc$CHR,'_', summary_statistics_hyprcoloc$START,'_',
                                                    summary_statistics_hyprcoloc$STOP,sep = '')
  significant_locus_coordinates <- significant_locus_coordinates[(order(significant_locus_coordinates$START)),]
  significant_locus_coordinates <- significant_locus_coordinates[(order(significant_locus_coordinates$CHR)),]
  summary_statistics_hyprcoloc <- summary_statistics_hyprcoloc[(order(summary_statistics_hyprcoloc$START)),]
  summary_statistics_hyprcoloc <- summary_statistics_hyprcoloc[(order(summary_statistics_hyprcoloc$CHR)),]
  summary_statistics_coloc_abf <- summary_statistics_coloc_abf[(order(summary_statistics_coloc_abf$START)),]
  summary_statistics_coloc_abf <- summary_statistics_coloc_abf[(order(summary_statistics_coloc_abf$CHR)),]
  
  summary_statistics_coloc_abf <- subset(summary_statistics_coloc_abf, COORDINATES %in% significant_locus_coordinates$COORDINATES)
  summary_statistics_hyprcoloc <- subset(summary_statistics_hyprcoloc, COORDINATES %in% significant_locus_coordinates$COORDINATES)
  visual_inspection_results <- subset(visual_inspection_results, COORDINATES %in% significant_locus_coordinates$COORDINATES)
  significant_locus_coordinates$COLOC_ABF_LEAD_SNP <- summary_statistics_coloc_abf$LEAD_SNP_BF
  significant_locus_coordinates$HYPRCOLOC_LEAD_SNP <- summary_statistics_hyprcoloc$CANDIDATE_SNP
  significant_locus_coordinates$LOCUSCOMPARE_LEAD_SNP <- summary_statistics_coloc_abf$LEAD_SNP_P_COMBINED
  
  ordered_heatmap_COLOC_ABF_LEAD_SNP <- get_effect_directions(snp_list = significant_locus_coordinates$COLOC_ABF_LEAD_SNP,
                                t2d_gwas = T2DadjBMI,
                                bfd_gwas = WHRadjBMI,
                                bmi_gwas = BMI)
  ordered_heatmap_LOCUSCOMPARE_LEAD_SNP <- get_effect_directions(snp_list = significant_locus_coordinates$LOCUSCOMPARE_LEAD_SNP,
                                                              t2d_gwas = T2DadjBMI,
                                                              bfd_gwas = WHRadjBMI,
                                                              bmi_gwas = BMI)
  View(subset(ordered_heatmap_COLOC_ABF_LEAD_SNP, sign(T2DadjBMI) == -sign(WHRadjBMI)))
  View(subset(ordered_heatmap_LOCUSCOMPARE_LEAD_SNP, sign(T2DadjBMI) == -sign(WHRadjBMI)))
  write.delim(ordered_heatmap_LOCUSCOMPARE_LEAD_SNP,'ordered_heatmap.txt', row.names = TRUE)
}

setwd('/Users/ya8eb/Documents/Research/gwas')
variants_EAF_WHRadjBMI = subset(WHRadjBMI, Freq_Tested_Allele == 0 | Freq_Tested_Allele == 1)
variants_EAF_T2DadjBMI = subset(T2DadjBMI, EAF == 0 | EAF == 1)
variants_EAF_T2DadjBMI$SNP = paste(variants_EAF_T2DadjBMI$Chr,':',variants_EAF_T2DadjBMI$Pos, sep = '')
write.delim(variants_EAF_WHRadjBMI, 'variants_EAF_WHRadjBMI.txt')
write.delim(variants_EAF_T2DadjBMI, 'variants_EAF_T2DadjBMI.txt')
