library(RPostgreSQL)
library(EnsDb.Hsapiens.v79)

dir <- "../incoming/metaClustersByTF"
bed.files <- list.files(dir, "*.bed")
length(bed.files) # 1361

uniprot.names <- unlist(lapply(strsplit(bed.files, "_"), function(tokens) tokens[1]))
tbl.files <- data.frame(file=bed.files, name=uniprot.names, stringsAsFactors=FALSE)

tbl.ens <- select(EnsDb.Hsapiens.v79,
                  key=uniprot.names,
                  keytype="UNIPROTID",
                  columns=c("SYMBOL"))
tbl.ens$SYMBOL[match("T", tbl.ref$SYMBOL)] <- "TBXT"

tbl.ref <- merge(tbl.files, tbl.ens, by.x="name", by.y="UNIPROTID", all.x=TRUE, no.dups=TRUE)
dim(tbl.ref)
length(bed.files)
length(uniprot.names)  # 1361
dups <- which(duplicated(tbl.ref$file))
if(length(dups) > 0)
    tbl.ref <- tbl.ref[-dups,]
dim(tbl.ref)   # 1361 3

nas <- which(is.na(tbl.ref$SYMBOL))
length(nas)
tbl.ref$SYMBOL[nas] <- tbl.ref$name[nas]
nas <- which(is.na(tbl.ref$SYMBOL))
length(nas)

db <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="genereg2021", host="khaleesi")

   # probably want an empty table to start

dbGetQuery(db, "select count(*) from grtd")

max <- length(bed.files)

for(i in seq_len(max)){
    filename = bed.files[i]
    full.path <- file.path(dir, filename)
    file.exists(full.path)
    tbl <- read.table(full.path, sep="\t", header=FALSE, as.is=TRUE)[, 1:3]
    colnames(tbl) <- c("chrom", "start", "endpos")
    tbl$width <- 1 + tbl$end - tbl$start
    tf <- subset(tbl.ref, file==filename)$SYMBOL
    printf("--- %5d: %s", i, tf)
    tbl$tf <- tf
    dbWriteTable(db, "grtd", tbl, row.names=FALSE, append=TRUE)
    } # for file, 6 files to start


dbGetQuery(db, "select count(*) from grtd")      # 110,458,580
dbGetQuery(db, "select * from grtd limit 6;")
