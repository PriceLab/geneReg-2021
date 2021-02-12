library(EnsDb.Hsapiens.v79)
bb.files <- dir(pattern=".bb")
length(bb.files)
uniprot.names <- unlist(lapply(strsplit(bb.files, "_"), function(tokens) tokens[1]))

tbl.geneSymbols <- select(EnsDb.Hsapiens.v79,
                          key=uniprot.names,
                          keytype="UNIPROTID",
                          columns=c("SYMBOL"))
tbl.geneSymbols$SYMBOL[match("T", tbl.geneSymbols$SYMBOL)] <- "TBXT"

for(file in bb.files){
    out.file <- sub(".bb", ".bed", file)
    cmd <- sprintf("./bigBedToBed %s %s", file, out.file)
    printf("creating %s", out.file)
    system(cmd)
    } # for file

bed.files <- dir(pattern=".bed")
length(bb.files)   # 1361
length(bed.files)  # 1361

