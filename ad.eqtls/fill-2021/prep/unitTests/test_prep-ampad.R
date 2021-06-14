library(RUnit)
source("../prep-ampad.R")
#----------------------------------------------------------------------------------------------------
outputFile <- "ampadDigested"
ad <- ampAD.digester$new(input.filename="~/github/geneReg-2021/ad.eqtls/fill-2021/smallSamples/tcx_mayo.csv",
                         projectName="ampad-mayo",
                         tissue="tcx",
                         output.file.basename=outputFile,
                         verbose=TRUE)

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_extractCoreColumns()
    test_add.hg38()
    try_bigFile()
    
} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))
    checkTrue(all(c("ampAD.digester", "R6") %in% class(ad)))

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_extractCoreColumns <- function()
{
    message(sprintf("--- test_extractCorColumns"))

    ad$extractCoreColumns()
    intermediate.filename <- sprintf("%s.csv", outputFile)
    checkTrue(file.exists(intermediate.filename))
    tbl.core <- read.table(intermediate.filename, sep=",", as.is=TRUE, header=TRUE)
    checkEquals(dim(tbl.core), c(99, 6))
    checkEquals(colnames(tbl.core), c("chromosome", "snpLocation", "snpid", "gene", "geneSymbol", "pvalue"))

} # test_extractCoreColumns
#----------------------------------------------------------------------------------------------------
test_add.hg38 <- function()
{
    message(sprintf("--- test_add.hg38"))

    ad$add.hg38()
    final.output.file <- sprintf("%s.RData", outputFile)
    tbl <- get(load(final.output.file))
    checkEquals(dim(tbl), c(99, 9))
    checkEquals(colnames(tbl),
                c("chrom", "hg19", "hg38", "rsid", "pvalue", "ensg", "geneSymbol", "project", "tissue"))
    checkEquals(unique(tbl$geneSymbol), "WASH7P")
    checkEquals(unique(tbl$project), "ampad-mayo")
    checkEquals(unique(tbl$tissue), "tcx")
    checkEquals(length(unique(tbl$rsid)), nrow(tbl))
    checkEquals(length(unique(tbl$hg19)), nrow(tbl))
    checkEquals(length(unique(tbl$hg38)), nrow(tbl))
    checkEquals(lapply(tbl, class),
                list(chrom="character",
                     hg19="integer",
                     hg38="integer",
                     rsid="character",
                     pvalue="numeric",
                     ensg="character",
                     geneSymbol="character",
                     project="character",
                     tissue="character"))
} # test_add.hg38
#----------------------------------------------------------------------------------------------------
try_bigFile <- function()
{
    message(sprintf("--- try_bigFile"))
   
   raw.file <- "/proj/price4/cory/AMP-AD_eQTLs/TCX_Mayo_cis_eQTL_release.csv"
   ad <- ampAD.digester$new(input.filename=raw.file,
                            projectName="ampad-mayo",
                            tissue="tcx",
                            output.file.basename="big",
                            verbose=TRUE)


   ad$extractCoreColumns()
   ad$add.hg38()
   tbl.transformed <- get(load("big.RData"))
   checkEquals(dim(tbl.transformed), c(63606555, 9))
   checkEquals(lapply(tbl.transformed, class),
               list(chrom="character",
                    hg19="integer",
                    hg38="integer",
                    rsid="character",
                    pvalue="numeric",
                    ensg="character",
                    geneSymbol="character",
                    project="character",
                    tissue="character"))

} # try_bigFile
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
