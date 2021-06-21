library(RUnit)
source("../prep-gtex.R")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_extractCoreColumns()
    test_add.hg19.snplocs.rsid()
    
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
    checkEquals(dim(tbl), c(100, 8))
    checkEquals(colnames(tbl), c("chrom", "hg38", "pvalue", "ensg", "geneSymbol", "project", "tissue", "assay"))
    checkEquals(as.character(lapply(tbl, class)), 
                c("character", "integer", "numeric", "character", "character", "character", "character", "character"))
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
    checkEquals(dim(tbl), c(101, 10))
    checkEquals(colnames(tbl), c("chrom", "hg19", "hg38", "rsid", "pvalue",
                                 "ensg", "geneSymbol", "project", "tissue", "assay"))
    checkEquals(as.character(lapply(tbl, class)), 
                c("character", "integer","integer","character","numeric","character","character",
                  "character","character","character"))

    checkEquals(length(which(is.na(tbl$rsid))), 9)
    checkEquals(length(which(tbl$hg19==-1)), 16)

      #------------------------------------------------------------
      # use liftover on the failed hg19 lookups
      #------------------------------------------------------------
    
    gd$liftover.missing.hg19.pos()
    tbl <- gd$getCurrentTable()
    checkEquals(dim(tbl), c(101, 10))
    checkEquals(colnames(tbl), c("chrom", "hg19", "hg38", "rsid", "pvalue",
                                 "ensg", "geneSymbol", "project", "tissue", "assay"))
    checkEquals(as.character(lapply(tbl, class)), 
                c("character", "integer","integer","character","numeric","character","character",
                  "character","character","character"))

    checkEquals(length(which(is.na(tbl$rsid))), 9)
    checkEquals(length(which(tbl$hg19==-1)), 0)


      #------------------------------------------------------------------
      # do a reverse liftover, from hg19 locs to hg39.  everything work?
      #------------------------------------------------------------------

    gr.hg19 <- GRanges(data.frame(chrom=tbl$chrom, start=tbl$hg19, end=tbl$hg19, stringsAsFactors=FALSE))
    if(!file.exists("hg19ToHg38.over.chain")){
       message(sprintf("downloading liftover chain file"))
       system("curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz")
       system("gunzip hg19ToHg38.over.chain.gz")
       }
    chain <- import.chain("hg19ToHg38.over.chain")
    gr.list <-liftOver(gr.hg19, chain)
    gr.hg38 <- unlist(gr.list)
    tbl.hg38 <- as.data.frame(gr.hg38)[, c("seqnames", "start")]

      #------------------------------------------------------------------
      # 3 snps fail the reverse mapping
      #------------------------------------------------------------------
    tbl[which(tbl$hg38 != tbl.hg38$start),]

       #    chrom      hg19      hg38       rsid      pvalue            ensg    geneSymbol project tissue   assay      signature
       # 42  chr9  44195357  42775762 rs28586892 9.08484e-09 ENSG00000170165    CR848007.2    GTEx    ctx unknown  chr9:42775762
       # 60  chr1 144044255 121315228 rs56399439 5.05606e-06 ENSG00000226067  CH17-118O6.2    GTEx    ctx unknown chr1:121315228
       # 93  chr1 149702931 121039576       <NA> 7.82282e-06 ENSG00000274642 CH17-472G23.1    GTEx    ctx unknown chr1:121039576

} # test_add.hg19.snplocs.rsid
#----------------------------------------------------------------------------------------------------
test_somePuzzlingResults <- function()
{

    # "58",ENSG00000170165.5,chr9_42775762_C_A_b38,-125476,0.142077,47,52,9.08484e-09,-0.79125,0.129563,0.000145102,2.23868e-09,2.10386e-06
    #    chrom     hg38      pvalue            ensg geneSymbol project tissue   assay
    # 58     9 42775762 9.08484e-09 ENSG00000170165 CR848007.2    GTEx    ctx unknown
    x <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP150.GRCh38, GRanges(seqnames="9", IRanges(start=42775762)))
    # x$RefSNP_id  [1] "rs28586892"
    subset(tbl, hg38==42775762)
    #    chrom     hg19     hg38       rsid      pvalue            ensg geneSymbol project tissue   assay     signature
    # 42  chr9 44195357 42775762 rs28586892 9.08484e-09 ENSG00000170165 CR848007.2    GTEx    ctx unknown chr9:42775762
    checkEquals(start(snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, "rs28586892")), 44195357)

    gr.hg19 <- GRanges(seqnames="chr9", IRanges(start=44195357))
    if(!file.exists("hg19ToHg38.over.chain")){
       message(sprintf("downloading liftover chain file"))
       system("curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz")
       system("gunzip hg19ToHg38.over.chain.gz")
       }
    chain <- import.chain("hg19ToHg38.over.chain")
    gr.list <-liftOver(gr.hg19, chain)
    gr.hg38 <- unlist(gr.list)
    checkEquals(start(gr.hg38), 42518796)   # but we expect the hg38 loc we started with, 42775762
    
} # test_somePuzzlingResults
#----------------------------------------------------------------------------------------------------
