library(annotatr)
load("annotations.RData")
tbl <- get(load("human_RB_binding_sites.RData"))
dim(tbl)
gata2.loc <- list(chrom="chr3", start=128479303, end=128493361, string="chr3:128,479,303-128,493,361")
tbl.gata2 <- subset(tbl, chrom=="chr3" & start >= gata2.loc$start & end <= gata2.loc$end)
length(unique(tbl.gata2$gene))  #  105
fivenum(tbl.gata2$score)  # [1]    0.5533157    4.0622183    7.0667659   13.6456034 1446.0000000

iCLIP,CIMS                         HeLa GSE47976,GSM1163973  11$
 18983895  chr3 128479553 128479554 human_RBP_CLIPdb_18944619      -   UPF1            iCLIP,CIMS                         HeLa GSE47976,GSM1163974    $
 > table(tbl.gata2.hi$celltype)

                       HEK293                      HEK293T                         HeLa HeLa,treated_with_puromycoin                        HepG2     $
                           23                          115                           27                            4                           27     $
       TREX_FLP-In_293T_cells
                            2
as.data.frame(sort(table(tbl.gata2.hi$celltype)))
hi <- 10

tbl.gata2.hi <- subset(tbl.gata2, score >= hi & celltype=="K562")
dim(tbl.gata2.hi)
length(unique(tbl.gata2.hi$gene)) 
as.data.frame(sort(table(tbl.gata2.hi$celltype)))
#                           Var1 Freq
#  1       TREX_FLP-In_293T_cells    2
#  2 HeLa,treated_with_puromycoin    4
#  3                      SH-SY5Y   20
#  4                       HEK293   23
#  5                         HeLa   27
#  6                        HepG2   27
#  7                      HEK293T  115
#  8                         K562  325
# The K562 leukemia cell has properties of self-renewal and pluripotency similar to those of the hematopoietic stem cell. (1981)



tbl.freq <- as.data.frame(sort(table(tbl.gata2.hi$gene)))
#     Var1 Freq
# 1   PUM2    3
# 2   UPF1    3
# 3 YTHDF1    3
# 4 YTHDF2    3
# 5  RBM15    5
# 6 RBM15B   10

mtx <- get(load(system.file(package="TrenaProjectErythropoiesis", "extdata","expression", "brandLabDifferentiationTimeCourse-27171x28.RData")))
cor(mtx["GATA2",], mtx["RBM15B",]) # [1] 0.879514

correlations <- lapply(tbl.freq$Var1, function(g) cor(mtx["GATA2",], mtx[g,]))
names(correlations) <- tbl.freq$Var1
tbl.cor <- with(correlations, data.frame(rbp=names(correlations), cor=as.numeric(correlations), stringsAsFactors=FALSE))
tbl.cor2 <- merge(tbl.cor, tbl.freq, by.x="rbp", by.y="Var1")
tbl.cor2 <- tbl.cor2[order(abs(tbl.cor2$cor), decreasing=TRUE),]

# pnas 2020, knockdown of AKAP8L suppressed the commitment of HSCs to
# erythroid lineage and cell proliferation and delayed differentiation
# of colony-forming unit-erythroid (CFU-E) to the proerythroblast
# stage (ProE). In contrast, the knockdown of TERF2IP and RNF10
# delayed differentiation of PolyE to OrthoE stage.

#        rbp         cor Freq
# 2   AKAP8L  0.90053214    6
# 28   SF3B4  0.88736206    9
# 16  GTF2F1  0.86238025    5
# 23    PUM2  0.86190922   11
# 15   GPKOW  0.85260924    1
# 12 FAM120A  0.80718394    3
# 22   PRPF8  0.80666300   18
# 29    SLTM  0.78954978   16
# 18 IGF2BP2  0.78633508    1
# 11   EWSR1  0.76965081   34
# 7    DHX30  0.76270780    6
# 3    BUD13  0.73830350    1
# 4   CSTF2T  0.71107711   24
# 20   LSM11  0.67619199    2
# 9   EFTUD2  0.66490277    1
# 21    NONO  0.63853038   25
# 38    YBX3  0.63014224    6
# 25   RBM15  0.59541195   21
# 37    XRN2 -0.44746049    1
# 34   U2AF2  0.43589307    3
# 10  EIF4G2  0.43499812   30
# 24  RBFOX2  0.38310617   10
# 31   TBRG4 -0.24374598    3
# 27  SERBP1  0.22765232    1
# 26   RBM27 -0.18536479    1
# 14  GEMIN5 -0.17579199    6
# 13 FASTKD2  0.14230294    1
# 36    UPF1  0.13622940   29
# 6     DDX6  0.12421667    5
# 1     AARS  0.11155971    2
# 32    TIA1 -0.11013400   12
# 30  TARDBP  0.10977076   12
# 35   UCHL5  0.08501428    4
# 33   U2AF1  0.05787377    1
# 5    DDX3X          NA    4
# 8   DROSHA          NA    2
# 17  HNRNPM          NA    4
# 19   KHSRP          NA    4
# 
