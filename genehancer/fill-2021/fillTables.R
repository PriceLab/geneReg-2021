library(RPostgreSQL)

db <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="gh50", host="khaleesi")

    #----------------------------
    # first, the elements table  
    #----------------------------

tbl <- read.table("../incoming/ISB/GeneHancer_elements.txt", sep="\t", as.is=TRUE, header=TRUE)
dim(tbl)  # 394565 6
colnames(tbl)
lapply(tbl, class)
colnames(tbl) <- c("chr", "element_start", "element_end", "ghid", "is_elite", "type")
dbWriteTable(db, "elements", tbl, row.names=FALSE, append=TRUE)
dbGetQuery(db, "select count(*) from elements")      # 394,565

    #----------------------------
    # now the TFBS table
    #----------------------------

tbl <- read.table("../incoming/ISB/GeneHancer_TFBSs.txt", sep="\t", as.is=TRUE, header=TRUE)
dim(tbl)  #  7,657,543       3
colnames(tbl)
lapply(tbl, class)
colnames(tbl) <- c("ghid", "tf", "tissues")
dbWriteTable(db, "tfbs", tbl, row.names=FALSE, append=TRUE)
dbGetQuery(db, "select count(*) from tfbs")      # 7,657,543


    #----------------------------
    # the tissues table
    #----------------------------

tbl <- read.table("../incoming/ISB/GeneHancer_tissues.txt", sep="\t", as.is=TRUE, header=TRUE)
dim(tbl)  #  5,725,390       4
colnames(tbl)
lapply(tbl, class)
colnames(tbl) <- c("ghid", "source", "tissue", "category")
dbWriteTable(db, "tissues", tbl, row.names=FALSE, append=TRUE)
dbGetQuery(db, "select count(*) from tissues")      # 5,725,390

    #----------------------------
    # gene associations
    #----------------------------

tbl <- read.table("../incoming/ISB/GeneHancer_gene_associations.txt", sep="\t", as.is=TRUE, header=TRUE)
dim(tbl)  #  2199479      10
colnames(tbl)
lapply(tbl, class)
colnames(tbl) <- c("ghid", "symbol", "eqtl_score", "erna_score", "chic_score",
                   "expression_score", "distance_score", "tss_proximity",
                   "combined_score", "is_elite")
dbWriteTable(db, "associations", tbl, row.names=FALSE, append=TRUE)
dbGetQuery(db, "select count(*) from associations")      # 2,199,479




