##############

setwd("C:/Users/user/Google Drive/Work/Semantic research framework/RFBR 2014, Aging gene ontology on invertebrates/GSEA/GSEA_vc01/Data/p02/fly_counts_FBtr")

X1C0M1 <- read.table(file="1C0M1.counts", sep="")
X2C0M1 <- read.table(file="2C0M1.counts", sep="")
X1C1M1 <- read.table(file="1C1M1.counts", sep="")
X2C1M1 <- read.table(file="2C1M1.counts", sep="")
X1C2M1 <- read.table(file="1C2M1.counts", sep="")
X2C2M1 <- read.table(file="2C2M1.counts", sep="")
X1C3M1 <- read.table(file="1C3M1.counts", sep="")
X2C3M1 <- read.table(file="2C3M1.counts", sep="")
X2F1M1 <- read.table(file="2F1M1.counts", sep="")
X3F1M1 <- read.table(file="3F1M1.counts", sep="")
X2F2M1 <- read.table(file="2F2M1.counts", sep="")
X3F2M1 <- read.table(file="3F2M1.counts", sep="")
X2F3M1 <- read.table(file="2F3M1.counts", sep="")
X3F3M1 <- read.table(file="3F3M1.counts", sep="")
X1R0 <- read.table(file="1R0.counts", sep="")
X2R0 <- read.table(file="2R0.counts", sep="")
X1R1 <- read.table(file="1R1.counts", sep="")
X2R1 <- read.table(file="2R1.counts", sep="")
X1R2 <- read.table(file="1R2.counts", sep="")
X2R2 <- read.table(file="2R2.counts", sep="")
X1R3 <- read.table(file="1R3.counts", sep="")
X2R3 <- read.table(file="2R3.counts", sep="")
X1S0M1 <- read.table(file="1S0M1.counts", sep="")
X2S0M1 <- read.table(file="2S0M1.counts", sep="")
X1S1M1 <- read.table(file="1S1M1.counts", sep="")
X2S1M1 <- read.table(file="2S1M1.counts", sep="")

countTable <- data.frame("TRANSCRIPTS"=X1C0M1[,1],
                         "X1C0M1"=X1C0M1[,2],
                         "X2C0M1"=X2C0M1[,2],
                         "X1C1M1"=X1C1M1[,2],
                         "X2C1M1"=X2C1M1[,2],
                         "X1C2M1"=X1C2M1[,2],
                         "X2C2M1"=X2C2M1[,2],
                         "X1C3M1"=X1C3M1[,2],
                         "X2C3M1"=X2C3M1[,2],
                         "X2F1M1"=X2F1M1[,2],
                         "X3F1M1"=X3F1M1[,2],
                         "X2F2M1"=X2F2M1[,2],
                         "X3F2M1"=X3F2M1[,2],
                         "X2F3M1"=X2F3M1[,2],
                         "X3F3M1"=X3F3M1[,2],
                         "X1R0"=X1R0[,2],
                         "X2R0"=X2R0[,2],
                         "X1R1"=X1R1[,2],
                         "X2R1"=X2R1[,2],
                         "X1R2"=X1R2[,2],
                         "X2R2"=X2R2[,2],
                         "X1R3"=X1R3[,2],
                         "X2R3"=X2R3[,2],
                         "X1S0M1"=X1S0M1[,2],
                         "X2S0M1"=X2S0M1[,2],
                         "X1S1M1"=X1S1M1[,2],
                         "X2S1M1"=X2S1M1[,2])

setwd("C:/Users/user/Google Drive/Work/Semantic research framework/RFBR 2014, Aging gene ontology on invertebrates/GSEA/GSEA_vc01")
write.table(countTable, file = "./Data/p02/p02_all.csv",sep=",",col.names=TRUE,row.names=FALSE)

fbgbtr <- read.table(file="./Data/fbgn_fbtr_fbpp_fb_2014_06.csv", sep="", fill=TRUE)


#############
rm(list=ls())
loadPackages <- function() {
    if(!try(require("DESeq"))) stop("package affy is missing")
    if(!try(require("affy"))) stop("package affy is missing")
    if(!try(require("affyPLM"))) stop("package affyPLM is missing")
    if(!try(require("plier"))) stop("package plier is missing")
    if(!try(require("limma"))) stop("package limma is missing")
    if(!try(require("biomaRt"))) stop("package biomaRt is missing")
    if(!try(require("org.Dm.eg.db"))) stop("package org.Dm.eg.db is missing")
    if(!try(require("org.Hs.eg.db"))) stop("package org.Hs.eg.db is missing")
    if(!try(require("igraph"))) stop("igraph")
    if(!try(require("marray"))) stop("marray")
    if(!try(require("AnnotationDbi"))) stop("AnnotationDbi")
    if(!try(require("gplots"))) stop("gplots")
    if(!try(require("GO.db"))) stop("GO.db")
    if(!try(require("KEGG.db"))) stop("KEGG.db")
    if(!try(require("data.table"))) stop("data.table")
    if(!try(require("gProfileR"))) stop("package gProfileR is missing")
    if(!try(require("RColorBrewer"))) stop("package RColorBrewer is missing")
    if(!try(require("Category"))) stop("package Category is missing")
    if(!try(require("GOstats"))) stop("package GOstats is missing")
    if(!try(require("Matrix"))) stop("package Matrix is missing")
    
    source("./R code/loadData.R")
    source("./R code/extractFactors.R")
    source("./R code/plotPCA.R")
    source("./R code/plotBox.R")
    source("./R code/plotHist.R")
    source("./R code/exprDiff.R")
    source("./R code/loadSets.R")
    source("./R code/useKEGGDrivenOntology.R")
    source("./R code/plotVenn.R")
    source("./R code/enrichSets.R")
    source("./R code/plotNetwork.R")
    source("./R code/plotHeatmap.R")
    
    #     source("./R code/runQC.R")
    #     source("./R code/diffExp.R")
    source("./R code/checkLoadArg.R")
    source("./R code/GSCstatBatch.R")
}



setwd("C:/Users/user/gitrep/gsea")
loadPackages()


#############
# load data file, currently from the 
if (file.exists("./Data/datafile_p02.rds")) {
    l_da <- readRDS("./Data/datafile_p02.rds")
} else {
    #that's the function
    l_da <- loadData(datadir = getwd(), countfile = "/Data/p02_all.csv", setup="/Data/setup.txt")
    saveRDS(l_da,"./Data/datafile_p02.rds")
}


# quality checks and initial analysis
extractFactors(l_da)
plotPCA(l_da)
plotBox(l_da)
plotHist(l_da)

contrasts = c("temperature_C25_control - temperature_C4_H24",
              "temperature_C25_control - temperature_C0_H24",
              "temperature_C25_control - temperature_C-4_H24",
              "fungus_none_control - fungus_min_H24",
              "fungus_none_control - fungus_max_H24",
              "radiation_G0_control - radiation_G200_H48",
              "radiation_G0_control - radiation_G500_H48",
              "radiation_G0_control - radiation_G1200_H48",
              "starvation_none_control - starvation_exists_H16")


# express differences
dexp <- exprDiff(l_da, contrasts=contrasts)

# load gene set
# "biological_process" "molecular_function" "cellular_component" "rfbr_orthologs"  "kegg" "aging_diseases" "aging_pathways"
gsetgroup <- "biological_process"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "molecular_function"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "cellular_component"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "rfbr_orthologs"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "kegg"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "aging_diseases"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "aging_pathways"; gsets <- loadSets(dexp=dexp, group=gsetgroup)


# vennPlot
rm(plotVenn); source("./R code/plotVenn.R")

contrasts = c("temperature_C25_control - temperature_C4_H24","temperature_C25_control - temperature_C0_H24","temperature_C25_control - temperature_C-4_H24")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(data = dexp, transcript2gene = gsets, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

contrasts = c("fungus_none_control - fungus_min_H24", "fungus_none_control - fungus_max_H24")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(data = dexp, transcript2gene = gsets, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

contrasts = c("radiation_G0_control - radiation_G200_H48","radiation_G0_control - radiation_G500_H48","radiation_G0_control - radiation_G1200_H48")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(data = dexp, transcript2gene = gsets, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

contrasts = c("starvation_none_control - starvation_exists_H16")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(data = dexp, transcript2gene = gsets, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

# enrich sets
contrasts = c("temperature_C25_control - temperature_C4_H24","temperature_C25_control - temperature_C0_H24","temperature_C25_control - temperature_C-4_H24",
              "fungus_none_control - fungus_min_H24","fungus_none_control - fungus_max_H24",
              "radiation_G0_control - radiation_G200_H48","radiation_G0_control - radiation_G500_H48","radiation_G0_control - radiation_G1200_H48",
              "starvation_none_control - starvation_exists_H16")
dexp <- exprDiff(l_da, contrasts=contrasts)
rm(plotNetwork); source("./R code/plotNetwork.R")
pcut <- 0.1; type="padj"
for (contrast in contrasts) {
    if (file.exists(paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcut,".rds",sep=""))) {
        gsah<- readRDS(paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcut,".rds",sep=""))
    }
    else {
        gsah <- enrichSets(dexp = dexp, 
                           pvalues = dexp$padj[[contrast]], 
                           sets = gsets, 
                           pcutoff = pcut)
        saveRDS(gsah,paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcut,".rds",sep=""))
    }
    print(contrast)
    tryCatch ({
        plotNetwork(dexp = dexp, 
                    contrast = contrast,
                    sets = gsets,
                    enrset = gsah, 
                    pcutoff = pcut,
                    type = "padj", # either "pval" or "padj"
                    gsetgroup = gsetgroup) 
    }, error=function(e){cat("No network diagram for ", contrast, ", ", gsetgroup,"\n")})
}

control<-"temperature_C25_control - temperature_C4_H24"; gsah<- readRDS(paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcut,".rds",sep="")); dexp = dexp; contrast = contrast; sets = gsets; enrset = gsah; pcutoff = pcut; type = "padj"; gsetgroup = gsetgroup


gsetgroup <- "biological_process"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "molecular_function"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "cellular_component"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "rfbr_orthologs"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "kegg"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "aging_diseases"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "aging_pathways"; gsets <- loadSets(dexp=dexp, group=gsetgroup)

rm(plotHeatmap); source("./R code/plotHeatmap.R")
plotHeatmap(dexp = dexp, contrasts = contrasts, sets = gsets, gsetgroup = gsetgroup, type = "padj", pcutoff = pcut)

dexp = dexp; contrasts = contrasts; sets = gsets; gsetgroup = gsetgroup; type = "padj"; pcutoff = pcut


contrasts = c("temperature_C25_control - temperature_C4_H24","temperature_C25_control - temperature_C0_H24","temperature_C25_control - temperature_C-4_H24")







dexp = dexp; contrast = contrast; sets = gsets; enrset = gsah; pcutoff = pcut; type = "padj"; gsetgroup = gsetgroup


pcutoff = pcut <- 0.1; type = "padj"

gsetgroup <- "rfbr_orthologs"; sets <- gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "aging_diseases"; sets <- gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "aging_pathways"; sets <- gsets <- loadSets(dexp=dexp, group=gsetgroup)


contrast <- "temperature_C25_control - temperature_C4_H24"; gsah <- enrichSets(dexp = dexp, pvalues = dexp$padj[[contrast]], sets = gsets, pcutoff = pcut)
contrast <- "temperature_C25_control - temperature_C0_H24"; gsah <- enrichSets(dexp = dexp, pvalues = dexp$padj[[contrast]], sets = gsets, pcutoff = pcut)
contrast <- "temperature_C25_control - temperature_C-4_H24"; gsah <- enrichSets(dexp = dexp, pvalues = dexp$padj[[contrast]], sets = gsets, pcutoff = pcut)
contrast <- "fungus_none_control - fungus_min_H24"; gsah <- enrichSets(dexp = dexp, pvalues = dexp$padj[[contrast]], sets = gsets, pcutoff = pcut)
contrast <- "fungus_none_control - fungus_max_H24"; gsah <- enrichSets(dexp = dexp, pvalues = dexp$padj[[contrast]], sets = gsets, pcutoff = pcut)
contrast <- "radiation_G0_control - radiation_G200_H48"; gsah <- enrichSets(dexp = dexp, pvalues = dexp$padj[[contrast]], sets = gsets, pcutoff = pcut)
contrast <- "radiation_G0_control - radiation_G500_H48"; gsah <- enrichSets(dexp = dexp, pvalues = dexp$padj[[contrast]], sets = gsets, pcutoff = pcut)
contrast <- "radiation_G0_control - radiation_G1200_H48"; gsah <- enrichSets(dexp = dexp, pvalues = dexp$padj[[contrast]], sets = gsets, pcutoff = pcut)
contrast <- "starvation_none_control - starvation_exists_H16"; gsah <- enrichSets(dexp = dexp, pvalues = dexp$padj[[contrast]], sets = gsets, pcutoff = pcut)


pvalues <- dexp$padj[[contrast]]; pvalues[pvalues==0] <- 1e-10
dt <- data.frame("FBtranscriptID"=dexp$name, "padj"=pvalues, "lgfc"=dexp$lgFC[[contrast]])
dt$FBtranscriptID <- as.character(dt$FBtranscriptID)
dt[is.na(dt$padj),2] <- 1; dt[is.na(dt$lgfc),3] <- 0; sets[is.na(sets$ENTREZID)==TRUE,2] <- "unknown"; sets[is.na(sets$GSET)==TRUE,3] <- "unknown"
dtd <- merge(dt, sets,by="FBtranscriptID")
lst_n <- dtd[which(dtd$GSET %in% enrset[enrset$padj <= pcutoff,"gset"] & dtd$padj <= pcutoff & dtd$lgfc < 0),c(5,4)]; lst_p <- dtd[which(dtd$GSET %in% enrset[enrset$padj <= pcutoff,"gset"] & dtd$padj <= pcutoff & dtd$lgfc > 0),c(5,4)]
mlst_p <- merge(lst_p,lst_p,by="ENTREZID"); mlst_n <- merge(lst_n,lst_n,by="ENTREZID")
am_p <- get.adjacency(graph.edgelist(as.matrix(mlst_p[,2:3]), directed=FALSE)); am_n <- get.adjacency(graph.edgelist(as.matrix(mlst_n[,2:3]), directed=FALSE))
if (length(lst_p[,1])>0) {dev.new(); plot(graph.adjacency(am_p, mode="min",weighted=TRUE,add.rownames = TRUE,diag=FALSE),main=paste(contrast, "up-regulated", sep=", "))}
if (length(lst_n[,1])>0) {dev.new(); plot(graph.adjacency(am_n, mode="min",weighted=TRUE,add.rownames = TRUE,diag=FALSE),main=paste(contrast, "down-regulated", sep=", "))}




