library(RUnit)
source("../prep-ampad.R")
#----------------------------------------------------------------------------------------------------
outputFile <- "ampadDigested"
small.sample.file <- "~/github/geneReg-2021/ad.eqtls/fill-2021/smallSamples/tcx_mayo.csv"
stopifnot(file.exists(small.sample.file))
small.sample.file.2 <- "~/github/geneReg-2021/ad.eqtls/fill-2021/smallSamples/tcx_mayo-multipleChromosomes.csv"
stopifnot(file.exists(small.sample.file.2))
ad <- ampAD.digester$new(input.filename=small.sample.file.2,
                         projectName="ampad-mayo",
                         tissue="tcx",
                         assay="unknown",
                         output.file.basename=outputFile,
                         verbose=TRUE)

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_extractCoreColumns()
    test_add.hg38.snplocs()

    test_rs28525262()

    #try_bigFile()
    
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
    intermediate.filename <- sprintf("%s-tmp.csv", outputFile)
    checkTrue(file.exists(intermediate.filename))
    tbl.core <- read.table(intermediate.filename, sep=",", as.is=TRUE, header=TRUE)
    checkEquals(dim(tbl.core), c(100, 6))
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
test_add.hg38.snplocs <- function()
{
    message(sprintf("--- test_add.hg38.locs"))

    outputFile <- "snplocs.100.variants"
    intermediate.file <- sprintf("%s-tmp.csv", outputFile)

    if(file.exists(intermediate.file)) # ensure a clean start
       unlink(intermediate.file)

    ad <- ampAD.digester$new(input.filename=small.sample.file.2,
                             projectName="ampad-mayo",
                             tissue="tcx",
                             assay="unknown",
                             output.file.basename=outputFile,
                             verbose=TRUE)
    ad$extractCoreColumns()
    t <- system.time(ad$add.hg38.snplocs())
    print(t)
    final.output.file <- sprintf("%s.RData", outputFile)
    tbl <- get(load(final.output.file))
    checkEquals(dim(tbl), c(100, 10))
    checkEquals(colnames(tbl),
                c("chrom", "hg19", "hg38", "rsid", "pvalue", "ensg", "geneSymbol", "project", "tissue", "assay"))
    top.10.expected.genes <- c("", "ANAPC4", "ANKRD65", "APOL1", "B9D1", "BLM",
                               "BORCS5", "BTG2", "C9orf139", "CALHM2")
    checkEquals(top.10.expected.genes, head(sort(unique(tbl$geneSymbol)), n=10))
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
                     tissue="character",
                     assay="character"))

        # one rsid in this bunch has no hg38 location
        # seems to be retired. dnsnp search says both
        # "No results" and
        # "G/T single-nucleotide variation in the ATP10B gene on chromosome 5"

     checkTrue(is.na(subset(tbl, rsid=="rs113887609")$hg38))

} # test_add.hg38.snplocs
#----------------------------------------------------------------------------------------------------
test_rs28525262 <- function()
{
   message(sprintf("--- test_rs28525262"))
    
   ad <- ampAD.digester$new(input.filename="tbl.rs28525262.csv",
                            projectName="testing",
                            tissue="tcx",
                            assay="unknown",
                            output.file.basename="tbl.rs28525262",
                            verbose=TRUE)
   ad$extractCoreColumns()
   ad$add.hg38.snplocs()
    
   tbl <- get(load("tbl.rs28525262.RData"))

} #     test_rs28525262
#----------------------------------------------------------------------------------------------------
try_bigFile <- function()
{
   message(sprintf("--- try_bigFile"))
   
   raw.file <- "/proj/price4/cory/AMP-AD_eQTLs/TCX_Mayo_cis_eQTL_release.csv"
   ad <- ampAD.digester$new(input.filename=raw.file,
                            projectName="ampad-mayo",
                            tissue="tcx",
                            assay="unknown",
                            output.file.basename="big",
                            verbose=TRUE)


   ad$extractCoreColumns()
   ad$add.hg38.snplocs()
   #print(t)

   tbl.transformed <- get(load("big.RData"))
   checkEquals(dim(tbl.transformed), c(63676897, 10))
   checkEquals(lapply(tbl.transformed, class),
               list(chrom="character",
                    hg19="integer",
                    hg38="integer",
                    rsid="character",
                    pvalue="numeric",
                    ensg="character",
                    geneSymbol="character",
                    project="character",
                    tissue="character",
                    assay="character"))

} # try_bigFile
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
