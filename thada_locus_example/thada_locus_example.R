# load libraries ----------------------------------------------------------
library(GenomicRanges)
#library(genomation)
library(ensembldb)
library(regioneR)
library(EnsDb.Hsapiens.v86)
library(Gviz)
library(rtracklayer)
library(trackViewer)
library(Sushi)
library(locuscomparer)

# THADA eQTL LocusZoom Plot  ---------------------------------------------------------------
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_qtl_colocalization/formatted_data/STARNET')
THADA_AS_eQTL <- read.delim('VAF.txt')
THADA_AS_eQTL <- subset(THADA_AS_eQTL, GENE == 'ENSG00000234936')
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/gwas_qtl_colocalization/formatted_data/GWAS')
T2DadjBMI <- read.delim('t2dadjbmi_gwas_in_loci.txt')
T2DadjBMI <- subset(T2DadjBMI, RSID %in% THADA_AS_eQTL$SNP & !duplicated(RSID))
T2DadjBMI <- data.frame(rsid = T2DadjBMI$RSID,
                        pval = T2DadjBMI$P)
WHRadjBMI <- read.delim('whradjbmi_gwas_in_loci.txt')
WHRadjBMI <- subset(WHRadjBMI, RSID %in% THADA_AS_eQTL$SNP & !duplicated(RSID))
WHRadjBMI <- data.frame(rsid = WHRadjBMI$RSID,
                        pval = WHRadjBMI$P)
THADA_AS_eQTL <- subset(THADA_AS_eQTL, SNP %in% WHRadjBMI$rsid)
THADA_AS_eQTL <- data.frame(rsid = THADA_AS_eQTL$SNP,
                            pval = THADA_AS_eQTL$P)
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/thada_locus_example/figures')
locuscomparer::locuscompare(in_fn1 = THADA_AS_eQTL, in_fn2 = WHRadjBMI, snp = 'rs6752964', title1 = 'THADA-AS VAF eQTL', title2 = 'WHRadjBMI GWAS', genome = 'hg19', population = 'EUR')
ggplot2::ggsave(filename = 'THADA-AS_VAF_WHRadjBMI.png',
                units = 'in',
                height = 6,
                width = 7)
locuscomparer::locuscompare(in_fn1 = THADA_AS_eQTL, in_fn2 = T2DadjBMI, snp = 'rs6752964', title1 = 'THADA-AS VAF eQTL', title2 = 'T2DadjBMI GWAS', genome = 'hg19', population = 'EUR')
ggplot2::ggsave(filename = 'THADA-AS_VAF_T2DadjBMI.png',
                units = 'in',
                height = 6,
                width = 7)
locuscomparer::locuscompare(in_fn1 = T2DadjBMI, in_fn2 = WHRadjBMI, snp = 'rs6752964', title1 = 'T2DadjBMI GWAS', title2 = 'WHRadjBMI GWAS', genome = 'hg19', population = 'EUR')
ggplot2::ggsave(filename = 'T2DadjBMI_WHRadjBMI.png',
                units = 'in',
                height = 6,
                width = 7)


# rs6752964 genome track plot  --------------------------------------------------------------
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/tactical/raw_data/tissue_annotation')
lead_snp <- data.frame(SNP = 'rs6752964',
                       CHR = 2,
                       POS = 43754905)
bed_files <- list.files()
start_position = 43.2e06
stop_position =44e06
Adipose <- read.delim(file = bed_files[1], header = F)
Adipose <- subset(Adipose, V1 == 'chr2' & V3 > start_position & V2 < stop_position)
HippocampusMiddle <- read.delim(file = bed_files[2], header = F)
HippocampusMiddle <- subset(HippocampusMiddle, V1 == 'chr2' & V3 > start_position & V2 < stop_position)
Islets <- read.delim(file = bed_files[3], header = F)
Islets <- subset(Islets, V1 == 'chr2' & V3 > start_position & V2 < stop_position)
Liver <- read.delim(file = bed_files[4], header = F)
Liver <- subset(Liver, V1 == 'chr2' & V3 > start_position & V2 < stop_position)
SkeletalMuscle <- read.delim(file = bed_files[5], header = F)
SkeletalMuscle <- subset(SkeletalMuscle, V1 == 'chr2' & V3 > start_position & V2 < stop_position)

gene_ids_of_interest <- c('ENSG00000234936','ENSG00000152518','ENSG00000115970')
genome_track <- GenomeAxisTrack()
ideogram_track <- IdeogramTrack(chromosome = "chr2", genome = "hg19")
gene_track <- BiomartGeneRegionTrack(genome = "hg19",
                                     name = "Genes",
                                     chromosome = 'chr2',
                                     start = start_position,
                                     end = stop_position,
                                     collapseTranscripts ='gene',
                                     transcriptAnnotation = 'gene')
highlight_track <- HighlightTrack(trackList = list(genome_track, gene_track),
                     start = lead_snp$POS, width = 1,
                     chromosome = lead_snp$CHR)
plotTracks(list(ideogram_track,highlight_track), from = start_position, to = stop_position)

Adipose <- read.delim(file = bed_files[1], header = F)
Adipose <- subset(Adipose, V1 == 'chr2' & V3 > 43754905 & V2 < 43754905)
HippocampusMiddle <- read.delim(file = bed_files[2], header = F)
HippocampusMiddle <- subset(HippocampusMiddle, V1 == 'chr2' & V3 > 43754905 & V2 < 43754905)
Islets <- read.delim(file = bed_files[3], header = F)
Islets <- subset(Islets, V1 == 'chr2' & V3 > 43754905 & V2 < 43754905)
Liver <- read.delim(file = bed_files[4], header = F)
Liver <- subset(Liver, V1 == 'chr2' & V3 > 43754905 & V2 < 43754905)
SkeletalMuscle <- read.delim(file = bed_files[5], header = F)
SkeletalMuscle <- subset(SkeletalMuscle, V1 == 'chr2' & V3 > 43754905 & V2 < 43754905)
repeats <- nrow(Adipose) + nrow(HippocampusMiddle) + nrow(Islets) + nrow(Liver) + nrow(SkeletalMuscle)
annotation_track <- AnnotationTrack(start = c(Adipose$V2,HippocampusMiddle$V2,Islets$V2,Liver$V2,SkeletalMuscle$V2),
                                    width = c(Adipose$V3-Adipose$V2,HippocampusMiddle$V3-HippocampusMiddle$V2,Islets$V3-Islets$V2,Liver$V3-Liver$V2,SkeletalMuscle$V3-SkeletalMuscle$V2),
                                    chromosome = 'chr2',
                                    strand = rep('*',repeats),
                                    group = rep(c(paste('Adipose',Adipose$V4),'Hippocampus','Pancreatic Islets','Liver','Skeletal Muscle'),
                                                c(nrow(Adipose),nrow(HippocampusMiddle),nrow(Islets),nrow(Liver),nrow(SkeletalMuscle))),
                                    genome = 'hg19',
                                    name = 'Chromatin States')
highlight_track <- HighlightTrack(trackList = list(genome_track, gene_track,annotation_track),
                                  start = lead_snp$POS, width = 1,
                                  chromosome = lead_snp$CHR)

plotTracks(list(ideogram_track,highlight_track), from = lead_snp$POS-10000, to = lead_snp$POS+10000,featureAnnotation = "group")



# implement functions -----------------------------------------------------
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/thada_locus_example/raw_data')
adipocyte_peak <- read.delim('GSE110734_SGBS_adipocyte_representative_peaks.narrowPeak')

adipocyte_peak <- data.frame(chrom = adipocyte_peak$X.chrom,
                             start = adipocyte_peak$chromStart,
                             end = adipocyte_peak$chromEnd,
                             value = adipocyte_peak$score)
mapped_SNP <- subset(adipocyte_peak, start < 43754905 & end > 43754905 & chrom == 'chr2')
chrom = 'chr2'
chromstart = 43754905-650000
chromend = 43754905+350000
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/thada_locus_example/figures')
jpeg(filename = "ATAC_Seq_THADA_Adipocyte.png",
     units = 'in',
     width = 5,
     height = 2.5,
     res = 300)
plotBedgraph(signal = adipocyte_peak,
        chrom = chrom,
        chromstart = chromstart,
        chromend = chromend,
        ymax = 1.04)
mtext("Adipocyte ATAC-Seq Consensus Peaks")
axis(side = 2, las=2, tcl=.2)
dev.off()

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/thada_locus_example/raw_data')
genes_in_locus <- read.delim('genes_in_locus.txt')
index <- genes_in_locus$Gene.name == 'AC010883.5'
genes_in_locus$Gene.name[index] = 'THADA-AS'
genes_in_locus <- subset(genes_in_locus, Gene.name %in% c('THADA-AS','THADA'))

genes_in_locus <- data.frame(chrom = genes_in_locus$Chromosome.scaffold.name,
                             start = genes_in_locus$Exon.region.start..bp.,
                             stop = genes_in_locus$Exon.region.end..bp.,
                             gene = genes_in_locus$Gene.name,
                             score = ".",
                             strand = genes_in_locus$Strand)

setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/thada_locus_example/raw_data')
preadipocyte_peak <- read.delim('GSE110734_SGBS_preadipocyte_representative_peaks.narrowPeak')

preadipocyte_peak <- data.frame(chrom = preadipocyte_peak$X.chrom,
                             start = preadipocyte_peak$chromStart,
                             end = preadipocyte_peak$chromEnd,
                             value = preadipocyte_peak$score)
mapped_SNP <- subset(preadipocyte_peak, start < 43754905 & end > 43754905 & chrom == 'chr2')

chrom = 'chr2'
chromstart = 43754905-650000
chromend = 43754905+350000
setwd('/Users/ya8eb/Documents/Research/colocalization_analysis_whradjbmi_t2dadjbmi/thada_locus_example/figures')
jpeg(filename = "ATAC_Seq_THADA_Preadipocyte.png",
     units = 'in',
     width = 5,
     height = 2.5,
     res = 300)
plotBedgraph(signal = preadipocyte_peak,
             chrom = chrom,
             chromstart = chromstart,
             chromend = chromend,
             ymax = 3)
mtext("Preadipocyte ATAC-Seq Consensus Peaks")
axis(side = 2, las=2, tcl=.2)
dev.off()

jpeg(filename = "Genes_in_THADA_locus.png",
     units = 'in',
     width = 5,
     height = 3,
     res = 300)
plotGenes(
  genes_in_locus,
  chrom = chrom,
  chromstart = chromstart,
  chromend = chromend,
  plotgenetype = 'box',
  bentline = F,
  types = 'exon')
labelgenome(chrom,chromstart,chromend,n=5,scale="Mb")
dev.off()
