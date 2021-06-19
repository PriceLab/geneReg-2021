library(R6)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
#---------------------------------------------------------------------------------------

gtex.digester = R6Class("gtex.digester",
       #--------------------------------------------------------------------------------
       private = list(input.filename=NULL,
                      projectName=NULL,
                      tissue=NULL,
                      assay=NULL,
                      tbl=NULL,
                      output.file.basename=NULL,
                      intermediate.filename=NULL,
                      final.output.filename=NULL,
                      verbose=NULL
                      ),
       #--------------------------------------------------------------------------------
       public = list(
           initialize = function(input.filename, projectName, tissue, assay,
                                 output.file.basename,verbose=FALSE){
               stopifnot(file.exists(input.filename))
            
               private$input.filename <- input.filename
               private$output.file.basename <- output.file.basename

               private$intermediate.filename <- sprintf("%s-tmp.csv", private$output.file.basename)
               private$final.output.filename <- sprintf("%s.RData",   private$output.file.basename)

               private$projectName <- projectName
               private$tissue <- tissue
               private$assay  <- assay
               private$verbose <- verbose
               private$tbl <- fread(input.filename)
               },

           extractCoreColumns = function(){
               coi <- c("phenotype_id", "variant_id", "pval_nominal")
               tbl <- as.data.frame(private$tbl)[, coi]
               loc.tokens <- strsplit(tbl$variant_id, "_")
               chrom <- sub("chr", "", unlist(lapply(loc.tokens, "[", 1)))
               hg38 <- as.integer(unlist(lapply(loc.tokens, "[", 2)))
               tbl <- data.frame(chrom=chrom, hg38=hg38, pvalue=tbl$pval_nominal,
                                 project=private$projectName,
                                 tissue=private$tissue,
                                 assay=private$assay,
                                 stringsAsFactors=FALSE)
               private$tbl <- tbl
               },

           getCurrentTable = function(){
               private$tbl
               }, 

           add.hg19.snplocs.rsid = function(){
               if(private$verbose)
                   message(sprintf("freading %s", private$outputFilename))
               gr <- GRanges(seqnames=private$tbl$chrom, IRanges(private$tbl$hg38))
               tbl.snpLocs <- as.data.frame(snpsByOverlaps(SNPlocs.Hsapiens.dbSNP150.GRCh38, gr))
               tbl.snpLocs$signature <- paste(tbl.snpLocs$seqnames, tbl.snpLocs$pos, sep=":")
               tbl <- private$tbl
               tbl$signature <- paste(tbl$chrom, tbl$hg38, sep=":")
               tbl.01 <- merge(tbl, tbl.snpLocs, by="signature", all.x=TRUE)               
               tbl.new <- tbl.01[, c("chrom", "hg38", "RefSNP_id", "pvalue", "project", "tissue", "assay")]
               colnames(tbl.new)[3] <- "rsid"
               rsids <- tbl.new$rsid
               na.rsids <- which(is.na(rsids))
               print(length(na.rsids))
               if(length(na.rsids) > 0){
                   printf("--- removing %d na.rsids", length(na.rsids))
                   rsids <- rsids[(-na.rsids)]
                   }
               x <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, rsids, ifnotfound="drop")
               tbl.hg19 <- as.data.frame(x)
               colnames(tbl.hg19)[grep("pos", colnames(tbl.hg19))] <- "hg19"
               colnames(tbl.hg19)[grep("RefSNP_id", colnames(tbl.hg19))] <- "rsid"
               tbl.merged <- merge(tbl.new, tbl.hg19[, c("hg19", "rsid")], by="rsid", all.x=TRUE)
               coi <- c("chrom", "hg19", "hg38", "rsid", "pvalue", "project", "tissue", "assay")
               tbl <- tbl.merged[, coi]
               browser()
               missing.hg19.pos <- which(is.na(tbl$hg19))
               printf("--- missing.hg19.pos count: %d", length(missing.hg19.pos))
               if(length(missing.hg19.pos) > 0)
                   tbl$hg19[missing.hg19.pos] <- -1
               tbl$hg19 <- as.integer(tbl$hg19)
               private$tbl <- tbl
               } # add.hg19.snplocs

           ) # public

       ) # class gtex.digester
   #--------------------------------------------------------------------------------