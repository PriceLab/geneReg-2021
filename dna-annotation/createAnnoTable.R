library(GenomicRanges)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

all.anno <- grep("hg38", builtin_annotations(), v=TRUE)
#  [1] "hg38_genes_1to5kb"               "hg38_genes_promoters"           
#  [3] "hg38_genes_cds"                  "hg38_genes_5UTRs"               
#  [5] "hg38_genes_exons"                "hg38_genes_firstexons"          
#  [7] "hg38_genes_introns"              "hg38_genes_intronexonboundaries"
#  [9] "hg38_genes_exonintronboundaries" "hg38_genes_3UTRs"               
# [11] "hg38_genes_intergenic"           "hg38_cpg_islands"               
# [13] "hg38_cpg_shores"                 "hg38_cpg_shelves"               
# [15] "hg38_cpg_inter"                  "hg38_enhancers_fantom"          
# [17] "hg38_lncrna_gencode"             "hg38_basicgenes"
# [19] "hg38_cpgs"                      

aoi <- all.anno[c(1:11, 16, 19)]

annotations <- build_annotations(genome="hg38", annotations=aoi)
length(annotations) # 7569764
tbl.anno <- as.data.frame(annotations)
colnames(tbl.anno)[1] <- "chrom"
tbl.anno$chrom <- as.character(tbl.anno$chrom)
filename <- "tbl.anno-hg38-13-categories.RData"
printf("saving %d annotations, tbl.anno, to %s", nrow(tbl.anno), filename)
save(tbl.anno, file=filename)

