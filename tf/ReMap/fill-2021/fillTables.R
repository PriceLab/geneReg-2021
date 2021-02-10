library(RPostgreSQL)

tbl.anno <- get(load("../incoming/tbl.experimentAnno.RData"))
dim(tbl.anno)  # 5798 6
colnames(tbl.anno) <- c("id", "tf", "celltype", "treatment", "experiment", "url")
db <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="genereg2021", host="khaleesi")

   # probably want an empty table to start

dbGetQuery(db, "select count(*) from remapAnno")

   # try a small fill, just 6 rows

tbl.sub <- head(tbl.anno)
dbSendQuery(db, "copy remapAnno from stdin")
postgresqlCopyInDataframe(db, tbl.sub)

dbGetQuery(db, "select count(*) from remapAnno")
dbGetQuery(db, "select * from remapAnno limit 6;")

    # now delete those 6 rows, load them all
dbGetQuery(db, "delete from remapAnno")
dbGetQuery(db, "select count(*) from remapAnno")
dbSendQuery(db, "copy remapAnno from stdin")
postgresqlCopyInDataframe(db, tbl.anno)
dbGetQuery(db, "select count(*) from remapAnno")
dbDiscconect(db)


tbl.regions <- get(load("../incoming/remap2020_all_macs2_hg38_v1_0.RData"))
dim(tbl.regions)

colnames(tbl.regions)
coi <- c("chrom", "start", "end", "name", "score", "peakStart", "peakEnd")
tbl.regions.prepped <- tbl.regions[, coi]
colnames(tbl.regions.prepped) <- coi <- c("chrom", "start", "end", "name", "score", "peakstart", "peakend")


   # probably want an empty table to start


   # try a small fill, just 6 rows

db <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="genereg2021", host="khaleesi")
dbGetQuery(db, "select count(*) from remapRegions")
tbl.sub <- head(tbl.regions.prepped)
dbSendQuery(db, "copy remapRegions from stdin")
postgresqlCopyInDataframe(db, tbl.sub)

dbGetQuery(db, "select count(*) from remapRegions")
dbGetQuery(db, "select * from remapRegions")

    # now delete those 6 rows, load them all
dbGetQuery(db, "delete from remapRegions")
dbGetQuery(db, "select count(*) from remapRegions")
dbSendQuery(db, "copy remapRegions from stdin")
postgresqlCopyInDataframe(db, tbl.regions.prepped)
dbGetQuery(db, "select count(*) from remapRegions")  # 164,732,372
dbDiscconect(db)
