library(RUnit)
source("../prep-gtex.R")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test ctor"))

    f <- "../../incoming/brainSmallTest.csv"
    checkTrue(file.exists(f))
    gd <- gtex.digester$new(input.filename=f,
                            projectName="GTEx",
                            tissue="ctx",
                            assay="unknown",
                            output.file.basename="gtex-brain-cortex",
                            verbose=TRUE)
    checkTrue(all(c("gtex.digester", "R6") %in% class(gd)))

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_extractCoreColumns <- function()
{
    message(sprintf("--- test_extractCoreColumns"))

    f <- "../../incoming/brainSmallTest.csv"
    checkTrue(file.exists(f))
    gd <- gtex.digester$new(input.filename=f,
                            projectName="GTEx",
                            tissue="ctx",
                            assay="unknown",
                            output.file.basename="gtex-brain-cortex",
                            verbose=TRUE)
    gd$extractCoreColumns()

    tbl <- gd$getCurrentTable()
    checkEquals(dim(tbl), c(100, 6))
    checkEquals(colnames(tbl), c("chrom", "hg38", "pvalue", "project", "tissue", "assay"))
    checkEquals(as.character(lapply(tbl, class)), 
                c("character", "integer", "numeric", "character", "character", "character"))
    checkEquals(unique(tbl$project), "GTEx")
    checkEquals(unique(tbl$tissue), "ctx")
    checkEquals(unique(tbl$assay), "unknown")

} # test_extractCoreColumns
#----------------------------------------------------------------------------------------------------
test_add.hg19.snplocs.rsid <- function()
{
    message(sprintf("--- test_add.hg19.snplocs.rsid"))

    f <- "../../incoming/brainSmallTest.csv"
    checkTrue(file.exists(f))
    gd <- gtex.digester$new(input.filename=f,
                            projectName="GTEx",
                            tissue="ctx",
                            assay="unknown",
                            output.file.basename="gtex-brain-cortex",
                            verbose=TRUE)
    gd$extractCoreColumns()
    gd$add.hg19.snplocs.rsid()
    
    tbl <- gd$getCurrentTable()
    checkEquals(dim(tbl), c(101, 8))
    checkEquals(colnames(tbl), c("chrom", "hg19", "hg38", "rsid", "pvalue", "project", "tissue", "assay"))
    checkEquals(as.character(lapply(tbl, class)), 
       c("character", "integer", "integer", "character", "numeric", "character", "character", "character"))
    checkEquals(length(which(is.na(tbl$rsid))), 9)

      # next up: use liftover on the failed hg19 lookups
    
} # test_add.hg19.snplocs.rsid
#----------------------------------------------------------------------------------------------------
