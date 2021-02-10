names <- get(load("experimentNames.RData"))
length(names)    # 5798
head(names)

#  experiment  tf   celltype
# "ENCSR000AHD.CTCF.MCF-7"
# "ENCSR000AHF.TAF1.MCF-7"
# "ENCSR000AKB.CTCF.GM12878"
# "ENCSR000AKO.CTCF.K-562"
# "ENCSR000ALA.CTCF.endothelial_umbilical-vein"
# "ENCSR000ALJ.CTCF.keratinocyte"

names.sub <- sample(names, 5)
strsplit(names.sub, ".", fixed=TRUE)

tokens <- strsplit(names, ".", fixed=TRUE)

tbls <- lapply(tokens, function(x){
           data.frame(experiment=x[1], tf=x[2], cellName=x[3], stringsAsFactors=FALSE)
         })
tbl.names <- do.call(rbind, tbls)
dim(tbl.names)                        # 5798 3
head(tbl.names)

length(unique(tbl.names$experiment))  # 2855
length(unique(tbl.names$tf))          # 1165
length(unique(tbl.names$cellName))    # 1924

tbl.names$id <- names
coi <- c("id", "experiment", "tf", "cellName")
tbl.names <- tbl.names[, coi]

cellNames <- sort(unique(tbl.names$cellName))
# tokens <- strsplit(cellTypes, "_")

tbls.cellTypeAndTreatment <- lapply(cellNames, function(cellName){
    x <- strsplit(cellName, "_")[[1]]
    type = x[1]
    treatment <- NA
    if(length(x) > 1)
       treatment <- paste(x[2:length(x)], collapse=";")
    #browser()
    data.frame(cellTypeID=cellName, cellType=x[1], treatment=treatment)
    })

tbl.cellTypeAndTreatment <- do.call(rbind, tbls.cellTypeAndTreatment)
dim(tbl.cellTypeAndTreatment)   # 1924 3
dim(tbl.names)                  # 5798 4

head(tbl.cellTypeAndTreatment)
head(tbl.names)

tbl.x <- head(tbl.names)
tbl.x

head(merge(tbl.x, tbl.cellTypeAndTreatment, by.x="cellName", by.y="cellTypeID", all.y=FALSE))
tbl.out <- merge(tbl.names, tbl.cellTypeAndTreatment, by.x="cellName", by.y="cellTypeID")
dim(tbl.out)  # 5798 6
coi <- c("id", "tf", "cellType", "treatment", "experiment")
tbl.final <- tbl.out[, coi]
dim(tbl.final)   # 5798  5

tbl.final[sample(1:nrow(tbl.final), size=10),]

# construct urls for each of the experiment ids?

# https://www.encodeproject.org/experiments/ENCSR000EFI/
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61944
# https://ddbj.nig.ac.jp/DRASearch/study?acc=ERP000209
# https://www.google.com/search?q=ERR063469

erp <- grep("^ERP", tbl.final$experiment)
length(erp)  # 158
enc <- grep("^ENC", tbl.final$experiment)
length(enc) #  1960
gse <- grep("^GSE", tbl.final$experiment)
length(gse) # 3679

 # just one experiment falls outside these three groups: ERR063469
 # use the google url for that
tbl.final[setdiff(seq_len(nrow(tbl.final)), sort(c(erp, enc, gse))),]

err <- grep("^ERR", tbl.final$experiment)

urls <- rep("", nrow(tbl.final))

urls[erp] <- sprintf("https://ddbj.nig.ac.jp/DRASearch/study?acc=%s", tbl.final$experiment[erp])
urls[enc] <- sprintf("https://www.encodeproject.org/experiments/%s/", tbl.final$experiment[enc])
urls[gse] <- sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s", tbl.final$experiment[gse])
urls[err] <- sprintf("https://www.google.com/search?q=%s", tbl.final$experiment[err])

which(nchar(urls) == 0)

tbl.final$url <- urls

tbl.experimentAnno <- tbl.final
save(tbl.experimentAnno, file="tbl.experimentAnno.RData")

