library(R6)
library(GenomicRanges)
library(rtracklayer)
library(data.table)

ampAD.digester = R6Class("ampAD.digester",
       #--------------------------------------------------------------------------------
       private = list(input.filename=NULL,
                      projectName=NULL,
                      tissue=NULL,
                      output.file.basename=NULL,
                      intermediate.filename=NULL,
                      final.output.filename=NULL,
                      rowCount=NULL,
                      verbose=NULL
                      ),
       #--------------------------------------------------------------------------------
       public = list(
           initialize = function(input.filename, projectName, tissue, output.file.basename,
                                 verbose=FALSE){
               stopifnot(file.exists(input.filename))
            
               private$input.filename <- input.filename
               private$output.file.basename <- output.file.basename

               private$intermediate.filename <- sprintf("%s-tmp.csv", private$output.file.basename)
               private$final.output.filename <- sprintf("%s.RData",   private$output.file.basename)

               private$projectName <- projectName
               private$tissue <- tissue
               private$verbose <- verbose
               },

           extractCoreColumns = function(){
               intermediate.filename <- sprintf("%s-tmp.csv", private$output.file.basename)
               cmd <- sprintf("cat %s  | cut -d ',' -f 1-5,7 > %s", private$input.filename,
                              intermediate.filename)
               if(private$verbose)
                   message(sprintf("system(%s)", cmd))
               system(cmd)
               invisible(cmd)
               },

           add.hg38 = function(){
               if(private$verbose)
                   message(sprintf("freading %s", private$outputFilename))
               tbl <- fread(private$intermediate.filename, sep=",", showProgress=FALSE)
               colnames(tbl) <- c("chrom", "start", "rsid", "ensg", "geneSymbol", "pvalue")
               tbl$chrom <- paste0("chr", tbl$chrom)
               tbl$end <- tbl$start
               tbl$hg19 <- tbl$start
               gr <- GRanges(tbl)
               if(!file.exists("hg19ToHg38.over.chain")){
                  message(sprintf("downloading liftover chain file"))
                  system("curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz")
                  system("gunzip hg19ToHg38.over.chain.gz")
                  }
               chain <- import.chain("hg19ToHg38.over.chain")
               if(private$verbose) message(sprintf("--- running liftover on %d rows", nrow(tbl)))
               gr.list <-liftOver(gr, chain)
               gr.hg38 <- unlist(gr.list)
               tbl <- as.data.frame(gr.hg38)
               colnames(tbl) <- c("chrom", "hg38", "drop.end", "drop.width", "drop.strand", "rsid",
                                  "ensg", "geneSymbol", "pvalue", "hg19")
               coi <- c("chrom", "hg19", "hg38", "rsid", "pvalue", "ensg", "geneSymbol")
               tbl <- tbl[, coi]
               tbl$chrom <- as.character(tbl$chrom)
               tbl$project <- private$projectName
               tbl$tissue  <- private$tissue
               if(private$verbose)
                   message(sprintf("--- writing %d rows, %d columns to %s",
                                   nrow(tbl), ncol(tbl), private$final.output.filename))
               save(tbl, file=private$final.output.filename)
               }
           ) # public

       ) # class ampAD.digester
   #--------------------------------------------------------------------------------
