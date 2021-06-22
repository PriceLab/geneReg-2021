library(data.table)
f <- "eqtls/Brain_Cortex.v8.EUR.signif_pairs.txt"
tbl <- fread(f, header=TRUE)
dim(tbl)
indices <- sort(sample(seq_len(nrow(tbl)), 100))
tbl.sub <- tbl[indices,]
fwrite(tbl.sub, file="brainSmallTest.csv", row.names=TRUE)


f <- "eqtls/Brain_Cortex.v8.EUR.signif_pairs.txt"
tbl <- fread(f, header=TRUE)
dim(tbl)
indices <- sort(sample(seq_len(nrow(tbl)), 1000))
tbl.sub <- tbl[indices,]
fwrite(tbl.sub, file="brain-1k-Test.csv", row.names=TRUE)
