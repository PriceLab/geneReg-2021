library(data.table)
tbl <- fread(file="human_RBP_binding_sites.txt", sep="\t", header=FALSE, nrow=-1)
dim(tbl) # 31,190,710 x 11
tbl <- tbl[, -5]
colnames(tbl) <- c("chrom", "start", "end", "name", "strand", "gene", "method", "celltype", "accession", "score")
dim(tbl) #  31,190,710       10
fivenum(tbl$score) # 0.5   -  4,769,520

subset(tbl, score > 1000000)
#    chrom    start      end                      name score1 strand   gene     method celltype           accession  score2
# 1: chr11 45826744 45826745  human_RBP_CLIPdb_2408095      0      - YTHDC1 iCLIP,CIMS  HEK293T GSE78030,GSM2064708 4769520
# 2: chr11 45826744 45826745  human_RBP_CLIPdb_2408097      0      - YTHDF1 iCLIP,CIMS  HEK293T GSE78030,GSM2064705 3950281
# 3: chr11 45826744 45826745  human_RBP_CLIPdb_2408098      0      - YTHDF2 iCLIP,CIMS  HEK293T GSE78030,GSM2064706 1017741
# 4: chr16 34160673 34160674 human_RBP_CLIPdb_10232039      0      + YTHDC1 iCLIP,CIMS  HEK293T GSE78030,GSM2064708 1174987
# 5: chr16 34160673 34160674 human_RBP_CLIPdb_10232040      0      + YTHDC2 iCLIP,CIMS  HEK293T GSE78030,GSM2064709 1789309

save(tbl, file="human_RB_binding_sites.RData")
