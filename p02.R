##############

path = "/home/antonkulaga/denigma/gsea/"
#path = "C:/Users/user/Google Drive/Work/Semantic research framework/RFBR 2014, Aging gene ontology on invertebrates/GSEA/GSEA_vc01/"
dataFolder = paste(path,"Data/",sep ="")
flyData = paste(dataFolder,"p02/fly_counts_FBtr",sep ="")

setwd(flyData)

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

setwd(path)
write.table(countTable, file = "./Data/p02/p02_all.csv",sep=",",col.names=TRUE,row.names=FALSE)

fbgbtr <- read.table(file="./Data/fbgn_fbtr_fbpp_fb_2014_06.csv", sep="", fill=TRUE)


#############
rm(list=ls())
path = "/home/antonkulaga/denigma/gsea/" #variable for base folder
#path = "C:/Users/user/Google Drive/Work/Semantic research framework/RFBR 2014, Aging gene ontology on invertebrates/GSEA/GSEA_vc01/"
dataFolder = paste(path,"Data/",sep ="") #variable for data folder

#vector of packages names
packages <- c("DESeq","affy","affyPLM","plier","limma",
              "biomaRt","org.Dm.eg.db","org.Hs.eg.db","igraph",
              "marray","AnnotationDbi","gplots","GO.db","KEGG.db",
              "data.table","gProfileR","RColorBrewer","Category","GOstats","Matrix")
src = paste(path,"Rcode/",sep="")

#function that loads packages 
#and tries to download those that are not installed
loadPackages <- function(pkgs,src) {
  source("http://bioconductor.org/biocLite.R") #enable biocLite
  library()  
  for(lib in pkgs){     
    if( !require(character.only = TRUE,package = lib, quietly = TRUE) ) {
      print(paste("trying to install",lib,"..."))
      biocLite(lib, dependencies=TRUE)
      require(character.only = TRUE,package = lib)
    }
  }              
}

loadSources = function(){
  source("./Rcode/loadData.R")
  source("./Rcode/extractFactors.R")
  source("./Rcode/plotPCA.R")
  source("./Rcode/plotBox.R")
  source("./Rcode/plotHist.R")
  source("./Rcode/exprDiff.R")
  source("./Rcode/loadSets.R")
  source("./Rcode/useKEGGDrivenOntology.R")
  source("./Rcode/plotVenn.R")
  source("./Rcode/enrichSets.R")
  source("./Rcode/plotNetwork.R")
  source("./Rcode/plotHeatmap.R")
  
  #     source("./R code/runQC.R")
  #     source("./R code/diffExp.R")
  source("./Rcode/checkLoadArg.R")
  source("./Rcode/GSCstatBatch.R")  
}


#setwd("C:/Users/user/gitrep/gsea")
setwd(path)
loadPackages(packages,src)
loadSources()
conc = function(one,two) {paste(one,two,sep="")} #shorcut function for concatenation

#############

loadDataFile = function(name, originalName, setupName = "setup.txt") {
  f = conc(dataFolder,name)
  if(file.exists(f)) return (readRDS(f)) else  {
    d <- loadData(datadir = dataFolder, countfile = originalName, setup=setupName)
    saveRDS(d,f)
    return (d)
  }
}

l_da = loadDataFile("datafile_p02.rds","p02_all.csv")

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
rm(plotVenn); source("./Rcode/plotVenn.R")

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



gsetgroup <- "biological_process"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "molecular_function"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "cellular_component"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "rfbr_orthologs"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "kegg"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "aging_diseases"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsetgroup <- "aging_pathways"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
pcut <- 0.1; type="padj"


rm(plotNetwork); source("./Rcode/plotNetwork.R")
rm(plotHeatmap); source("./Rcode/plotHeatmap.R")


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




plotHeatmap(dexp = dexp, contrasts = contrasts, sets = gsets, gsetgroup = gsetgroup, type = "padj", pcutoff = pcut)




gsetgroup <- "molecular_function"; gsets <- loadSets(dexp=dexp, group=gsetgroup)
gsah<- readRDS(paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcut,".rds",sep=""))
dexp = dexp; contrast = contrast; sets = gsets; enrset = gsah; pcutoff = pcut; type = "padj"; gsetgroup = gsetgroup



