library(GenomicRanges)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  # annots <- "hg38_basicgenes"
  # annotations <- build_annotations(genome="hg38", annotations=annots)
  # save(annotations, file="annotations.RData")  # 27M

load("annotations.RData")
tbl <- get(load("human_RB_binding_sites.RData"))
dim(tbl) # 31190710       10
head(tbl, n=30)
tbl.test <- tbl[sort(sample(nrow(tbl), size=10)),]
gr.test <- GRanges(as.data.frame(tbl.test, row.names=NULL))
gr.test.annoResults <- annotate_regions(regions=gr.test, annotations=annotations, ignore.strand=TRUE, quiet=FALSE)
tbl.test <- as.data.frame(gr.test.annoResults, row.names=NULL)
dim(tbl.test)
tbl.freq <- as.data.frame(table(tbl.test$annot.symbol, tbl.test$gene))
tbl.freq <- subset(tbl.freq, Freq > 0)
dim(tbl.freq)
save(gr.test.annoResults, file="gr.test.annoResults.RData")

gr <- GRanges(as.data.frame(tbl, row.names=NULL))
gr.annoResults <- annotate_regions(regions=gr, annotations=annotations, ignore.strand=TRUE, quiet=FALSE)
length(gr.annoResults) # 233693995
save(gr.annoResults, file="gr.annoResults.RData")

tbl.anno <- as.data.frame(gr.annoResults, row.names=NULL)
dim(tbl.anno) # 233693995        21
save(tbl.anno, file="tbl.anno.results.RData")

coi <- c("seqnames", "start", "end", "gene", "annot.symbol")
dups <- which(duplicated(tbl.anno[, coi]))
tbl.anno.nodups <- tbl.anno[-dups,]
save(tbl.anno.nodups, file="tbl.anno.nodups.RData")
