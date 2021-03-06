load("gr.annoResults.RData")
length(gr.annoResults)
tbl.anno <- as.data.frame(gr.annoResults, row.names=NULL)
dim(tbl.anno)
coi <- c("seqnames", "start", "end", "gene", "annot.symbol")
dups <- which(duplicated(tbl.anno[, coi]))
length(dups)
coi <- c("chrom", "start", "end", "width", "strand", "gene", "method", "celltype", "accession", "score", "annot.symbol", "annot.type")
tbl.anno.nodups <- tbl.anno[coi, -dups,]
save(tbl.anno.nodups, file="tbl.anno.nodups.RData")
write.table(tbl.anno.nodups, file="../fill-2021/tblAnno-ready.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
