##############
# unfortunatelly, the original Piano is not fully suitable for my case
# so I re-used Piano scripts making them less generic yet fitting my data better
source("http://bioconductor.org/biocLite.R")
biocLite(c("DESeq", "affy", "affyPLM", "plier", "limma", "biomaRt", "org.Dm.eg.db", "igraph", "marray", "gplots"))

source("http://bioconductor.org/biocLite.R")
biocLite("KEGGREST")
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

#############
# contrasts: "base - change" or "control - modification" or "wild - mutant"
contrasts=c("Control_fm - Irradiation_fm", 
            "Control_fm - toluene_fm",
            "Control_fm - Dioxin_fm",
            "Control_fm - formaldehyde_fm",
            "Control_m - Irradiation_m", 
            "Control_m - toluene_m",
            "Control_m - Dioxin_m",
            "Control_m - formaldehyde_m")


contrasts=c("Control_fm - Irradiation_fm", 
            "Control_fm - toluene_fm",
            "Control_fm - Dioxin_fm",
            "Control_fm - formaldehyde_fm")
#############
# load data file, currently from the 
if (file.exists("./Data/datafile.rds")) {
    l_da <- readRDS("./Data/datafile.rds")
} else {
    #that's the function
    l_da <- loadData(datadir = getwd(), countfile = "/Data/all.csv", setup="/Data/setup.txt")
    saveRDS(l_da,"./Data/datafile.rds")
}



setwd("C:/Users/user/Google Drive/Work/Semantic research framework/RFBR 2014, Aging gene ontology on invertebrates/GSEA/GSEA_vc01")
loadPackages()

rm("loadData"); source("./R code/loadData.R")
if (file.exists("./Data/datafile_new.rds")) {
    l_da <- readRDS("./Data/datafile_new.rds")
} else {
    l_da <- loadData(datadir = getwd(), countfile = "/Data/all_stress.csv", setup="/Data/setup.txt")
    saveRDS(l_da,"./Data/datafile_new.rds")
}


# quality checks and initial analysis
extractFactors(l_da); plotPCA(l_da); plotBox(l_da); plotHist(l_da)


# express differences
rm("exprDiff"); source("./R code/exprDiff.R")
dexp <- exprDiff(l_da, contrasts=contrasts)

# load gene set
# "biological_process" "molecular_function" "cellular_component" "rfbr_orthologs"  "kegg" "aging_diseases" "aging_pathways"
gsetgroup <- "molecular_function"

rm("loadSets"); source("./R code/loadSets.R")
gsets <- loadSets(dexp=dexp, group=gsetgroup)


# vennPlot
rm("plotVenn"); source("./R code/plotVenn.R")
contrasts=c("Control_fm - Irradiation_fm", "Control_fm - toluene_fm", "Control_fm - Dioxin_fm", "Control_fm - formaldehyde_fm")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(data = dexp, transcript2gene = gsets, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

contrasts=c("Control_m - Irradiation_m", "Control_m - toluene_m", "Control_m - Dioxin_m", "Control_m - formaldehyde_m")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(data = dexp, transcript2gene = gsets, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)


# enrich sets
rm("enrichSets"); source("./R code/enrichSets.R")
rm("plotNetwork"); source("./R code/plotNetwork.R")
pcut <- 0.1
for (contrast in contrasts) {
    print(contrast)
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
    tryCatch ({
        dev.new()
        plotNetwork(dexp = dexp, 
                    contrast = contrast,
                    sets = gsets,
                    enrset = gsah, 
                    pcutoff = pcut,
                    type = "padj", # either "pval" or "padj"
                    gsetgroup = gsetgroup) 
    }, error=function(e){cat("ERROR in ", contrast, ": ", conditionMessage(e),"\n")})
}

rm("plotHeatmap"); source("./R code/plotHeatmap.R")
plotHeatmap(dexp = dexp,
            contrasts = contrasts,
            sets = gsets,
            gsetgroup = gsetgroup,
            pcutoff = pcut)

rm("plotNetwork"); source("./R code/plotNetwork.R")
pltN <- plotNetwork(dexp = dexp, contrast = contrast, sets = gsets, enrset = gsah,  pcutoff = pcut, type = "pval", gsetgroup = gsetgroup) 

















