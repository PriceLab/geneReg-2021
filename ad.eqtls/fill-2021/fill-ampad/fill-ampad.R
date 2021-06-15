source("src/prep-ampad.R")
directory <- "/proj/price4/cory/AMP-AD_eQTLs/"
files <- c("CER_Mayo_cis_eQTL_release.csv",
           "DLPFC_ROSMAP_cis_eQTL_release.csv",
           "TCX_Mayo_cis_eQTL_release.csv")

for(file in files){
   file.base <- sub(".csv", "", file, fixed=TRUE)
   tissue <- strsplit(file, "_")[[1]][1]
   full.path <- file.path(directory, file)
   ad <- ampAD.digester$new(input.filename=full.path,
                         projectName="ampad-mayo",
                         tissue=tissue,
                         output.file.basename=file.base,
                         verbose=TRUE)
   ad$extractCoreColumns()
   ad$add.hg38()
   }
