library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
load("human_RB_binding_sites-granges.RData")
length(gr)  # 31190710
load("annotations.RData")
length(annotations) # 3551273
gr.annoResults <- annotate_regions(regions=gr, annotations=annotations, ignore.strand=TRUE, quiet=FALSE)
save(gr.annoResults, file="gr.annoResults.RData")

