# cf1 <- read.table(file="./Data/Cf.merged.1.counts", sep="")
# cf2 <- read.table(file="./Data/Cf.merged.2.counts", sep="")
# cm1 <- read.table(file="./Data/Cm.merged.1.counts", sep="")
# cm2 <- read.table(file="./Data/Cm.merged.2.counts", sep="")
# df1 <- read.table(file="./Data/Df.merged.1.counts", sep="")
# df2 <- read.table(file="./Data/Df.merged.2.counts", sep="")
# dm1 <- read.table(file="./Data/Dm.merged.1.counts", sep="")
# dm2 <- read.table(file="./Data/Dm.merged.2.counts", sep="")
# ff1 <- read.table(file="./Data/Ff.merged.1.counts", sep="")
# ff2 <- read.table(file="./Data/Ff.merged.2.counts", sep="")
# fm1 <- read.table(file="./Data/Fm.merged.1.counts", sep="")
# fm2 <- read.table(file="./Data/Fm.merged.2.counts", sep="")
# rf1 <- read.table(file="./Data/Rf.merged.1.counts", sep="")
# rf2 <- read.table(file="./Data/Rf.merged.2.counts", sep="")
# rm1 <- read.table(file="./Data/Rm.merged.1.counts", sep="")
# rm2 <- read.table(file="./Data/Rm.merged.2.counts", sep="")
# tf1 <- read.table(file="./Data/Tf.ferged.1.counts", sep="")
# tf2 <- read.table(file="./Data/Tf.ferged.2.counts", sep="")
# tm1 <- read.table(file="./Data/Tm.merged.1.counts", sep="")
# tm2 <- read.table(file="./Data/Tm.merged.2.counts", sep="")
# 
# countTable <- data.frame("TRANSCRIPTS"=ff1[,1],
#                          "cf1"=cf1[,2],
#                          "cf2"=cf2[,2],
#                          "cm1"=cm1[,2],
#                          "cm2"=cm2[,2],
#                          "df1"=df1[,2],
#                          "df2"=df2[,2],
#                          "dm1"=dm1[,2],
#                          "dm2"=dm2[,2],
#                          "ff1"=ff1[,2],
#                          "ff2"=ff2[,2],
#                          "fm1"=fm1[,2],
#                          "fm2"=fm2[,2],
#                          "rf1"=rf1[,2],
#                          "rf2"=rf2[,2],
#                          "rm1"=rm1[,2],
#                          "rm2"=rm2[,2],
#                          "tf1"=tf1[,2],
#                          "tf2"=tf2[,2],
#                          "tm1"=tm1[,2],
#                          "tm2"=tm2[,2])
# write.table(countTable, file = "./Data/p03_01_data.csv",sep=",",col.names=TRUE,row.names=FALSE)

# p0301 <- read.table(file="./Data/p03_01_data.csv",sep=",",header=TRUE,fill=TRUE)
# p0302 <- read.table(file="./Data/p03_02_data.csv",sep=",",header=TRUE,fill=TRUE)
# p0301[,1] <- as.character(p0301[,1]); p0302[,1] <- as.character(p0302[,1])
# p03 <- merge(p0301,p0302,by="TRANSCRIPTS")
# dim(p0301); dim(p0302); dim(p03)
# 
# write.table(p03, file = "./Data/p03_data.csv",sep=",",col.names=TRUE,row.names=FALSE)

#############
rm(list=ls())
# path = "/home/antonkulaga/denigma/gsea/" #variable for base folder
path <- "C:/Users/user/gitrep/gsea"


dataFolder = paste(path,"/Data/",sep ="") #variable for data folder

#vector of packages names
packages <- c("DESeq","affy","affyPLM","plier","limma",
              "biomaRt","org.Dm.eg.db","org.Hs.eg.db","igraph",
              "marray","AnnotationDbi","gplots","GO.db","KEGG.db",
              "data.table","gProfileR","RColorBrewer","Category","GOstats","Matrix")
src = paste(path,"/Rcode/",sep="")

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
    source("./Rcode/plotManualNetwork.R")
    source("./Rcode/plotHeatmap.R")
}


setwd(path)
loadPackages(packages,src)
loadSources()
conc = function(one,two) {paste(one,two,sep="")} #shorcut function for concatenation

#############

l_da <- loadData(datadir = dataFolder, countfile = "p03_data.csv", setup="setup.txt")
saveRDS(l_da,"./Data/datafile_p03.rds")

loadDataFile = function(name, originalName, setupName = "setup.txt") {
    f = conc(dataFolder,name)
    if(file.exists(f)) return (readRDS(f)) else  {
        d <- loadData(datadir = dataFolder, countfile = originalName, setup=setupName)
        saveRDS(d,f)
        return (d)
    }
}

l_da = loadDataFile("datafile_p03.rds","p03_data.csv")

# quality checks and initial analysis
extractFactors(l_da)
plotPCA(l_da)
plotBox(l_da)
plotHist(l_da)


# express differences
contrasts = c("control_f - dioxin_f", "control_f - formaldehyde_f", "control_f - radiation_f", "control_f - toluene_f",
              "control_m - dioxin_m", "control_m - formaldehyde_m", "control_m - radiation_m", "control_m - toluene_m",
              "temperature_C25_control - temperature_C4_H24", "temperature_C25_control - temperature_C0_H24", "temperature_C25_control - temperature_C-4_H24",
              "fungus_none_control - fungus_min_H24", "fungus_none_control - fungus_max_H24",
              "radiation_G0_control - radiation_G200_H48", "radiation_G0_control - radiation_G500_H48", "radiation_G0_control - radiation_G1200_H48",
              "starvation_none_control - starvation_exists_H16")
dexp <- exprDiff(l_da, contrasts=contrasts)

# load gene set
# "biological_process" "molecular_function" "cellular_component" "rfbr_orthologs"  "kegg" "aging_diseases" "aging_pathways"

# vennPlot
contrasts = c"control_f - dioxin_f", "control_f - formaldehyde_f", "control_f - radiation_f", "control_f - toluene_f")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(dexp = dexp, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

contrasts = c("control_m - dioxin_m", "control_m - formaldehyde_m", "control_m - radiation_m", "control_m - toluene_m")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(dexp = dexp, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

contrasts = c("temperature_C25_control - temperature_C4_H24","temperature_C25_control - temperature_C0_H24","temperature_C25_control - temperature_C-4_H24")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(dexp = dexp, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

contrasts = c("fungus_none_control - fungus_min_H24", "fungus_none_control - fungus_max_H24")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(dexp = dexp, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

contrasts = c("control_m - radiation_m", "radiation_G0_control - radiation_G200_H48","radiation_G0_control - radiation_G500_H48","radiation_G0_control - radiation_G1200_H48")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(dexp = dexp, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

contrasts = c("starvation_none_control - starvation_exists_H16")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(dexp = dexp, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)


# enrich sets
# express differences
contrasts = c("control_f - dioxin_f", "control_f - formaldehyde_f", "control_f - radiation_f", "control_f - toluene_f",
              "control_m - dioxin_m", "control_m - formaldehyde_m", "control_m - radiation_m", "control_m - toluene_m",
              "temperature_C25_control - temperature_C4_H24", "temperature_C25_control - temperature_C0_H24", "temperature_C25_control - temperature_C-4_H24",
              "fungus_none_control - fungus_min_H24", "fungus_none_control - fungus_max_H24",
              "radiation_G0_control - radiation_G200_H48", "radiation_G0_control - radiation_G500_H48", "radiation_G0_control - radiation_G1200_H48",
              "starvation_none_control - starvation_exists_H16")
dexp <- exprDiff(l_da, contrasts=contrasts)



gsetgroup <- "biological_process"; gsets <- loadSets(group=gsetgroup)
gsetgroup <- "molecular_function"; gsets <- loadSets(group=gsetgroup)
gsetgroup <- "cellular_component"; gsets <- loadSets(group=gsetgroup)
gsetgroup <- "rfbr_orthologs"; gsets <- loadSets(group=gsetgroup); gsets <- gsets[gsets[,3]!="missing",]
gsetgroup <- "kegg"; gsets <- loadSets(group=gsetgroup)
gsetgroup <- "aging_diseases"; gsets <- loadSets(group=gsetgroup)
gsetgroup <- "aging_pathways"; gsets <- loadSets(group=gsetgroup)


pcutoff <- 0.1; type<-"padj"; shownodes="all"
pcutoff <- 0.05; type<-"pval"; shownodes="all"

source("./Rcode/plotManualNetwork.R")
plotManualNetwork(dexp = dexp, contrast = contrast, gsets = gsets, enrsets = enrsets, pcutoff = pcutoff, type = type, gsetgroup = gsetgroup, shownodes=shownodes)

# contrast <- "control_f - toluene_f"
# contrast <- "radiation_G0_control - radiation_G1200_H48"

gsetgroups <- c("biological_process", "molecular_function", "cellular_component", "rfbr_orthologs", "kegg", "aging_diseases", "aging_pathways")
multiEnrichSets <- function(gsetgroups) {
    for (gsetgroup in gsetgroups) {
        gsets <- loadSets(group=gsetgroup)
        
        for (contrast in contrasts) {
            if (file.exists(paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcutoff,".rds",sep=""))) {
                enrsets<- readRDS(paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcutoff,".rds",sep=""))
            } else {
                enrsets <- enrichSets(dexp = dexp, 
                                   pvalues = dexp$padj[[contrast]], 
                                   gsets = gsets, 
                                   pcutoff = pcutoff)
                saveRDS(enrsets,paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcutoff,".rds",sep=""))
            }
            print(contrast)
            tryCatch ({
                plotNetwork(dexp = dexp, 
                            contrast = contrast,
                            gsets = gsets,
                            enrsets = enrsets, 
                            pcutoff = pcutoff,
                            type = type, # either "pval" or "padj"
                            gsetgroup = gsetgroup,
                            shownodes=shownodes) 
            }, error=function(e){cat("No network diagram for ", contrast, ", ", gsetgroup,"\n")})
        }
        
    }
}

multiEnrichSets(gsetgroups)

for (contrast in contrasts) {
    if (file.exists(paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcutoff,".rds",sep=""))) {
        enrsets<- readRDS(paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcutoff,".rds",sep=""))
    } else {
        enrsets <- enrichSets(dexp = dexp, 
                           pvalues = dexp$padj[[contrast]], 
                           gsets = gsets, 
                           pcutoff = pcutoff)
        saveRDS(enrsets,paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcutoff,".rds",sep=""))
    }
    print(contrast)
    tryCatch ({
        plotNetwork(dexp = dexp, 
                    contrast = contrast,
                    gsets = gsets,
                    enrsets = enrsets, 
                    pcutoff = pcutoff,
                    type = type, # either "pval" or "padj"
                    gsetgroup = gsetgroup,
                    shownodes="all") 
    }, error=function(e){cat("No network diagram for ", contrast, ", ", gsetgroup,"\n")})
}


contrasts = c("control_f - dioxin_f", "control_f - formaldehyde_f", "control_f - radiation_f", "control_f - toluene_f",
              "control_m - dioxin_m", "control_m - formaldehyde_m", "control_m - radiation_m", "control_m - toluene_m")
dexp <- exprDiff(l_da, contrasts=contrasts)
pcutoff <- 0.1; type<-"padj"; shownodes="all"
gsetgroup <- "biological_process"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "molecular_function"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "cellular_component"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "kegg"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)

pcutoff <- 0.05; type="pval"; shownodes="all"
gsetgroup <- "biological_process"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "molecular_function"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "cellular_component"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "kegg"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)

contrasts = c("temperature_C25_control - temperature_C4_H24", "temperature_C25_control - temperature_C0_H24", "temperature_C25_control - temperature_C-4_H24",
              "fungus_none_control - fungus_min_H24", "fungus_none_control - fungus_max_H24",
              "radiation_G0_control - radiation_G200_H48", "radiation_G0_control - radiation_G500_H48", "radiation_G0_control - radiation_G1200_H48",
              "starvation_none_control - starvation_exists_H16")
dexp <- exprDiff(l_da, contrasts=contrasts)

pcutoff <- 0.1; type<-"padj"
gsetgroup <- "biological_process"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "molecular_function"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "cellular_component"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "kegg"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)

pcutoff <- 0.05; type="pval"
gsetgroup <- "biological_process"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "molecular_function"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "cellular_component"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "kegg"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)




contrasts = c("control_m - dioxin_m", "control_m - formaldehyde_m", "control_m - radiation_m", "control_m - toluene_m",
              "temperature_C25_control - temperature_C4_H24", "temperature_C25_control - temperature_C0_H24", "temperature_C25_control - temperature_C-4_H24",
              "fungus_none_control - fungus_min_H24", "fungus_none_control - fungus_max_H24",
              "radiation_G0_control - radiation_G200_H48", "radiation_G0_control - radiation_G500_H48", "radiation_G0_control - radiation_G1200_H48",
              "starvation_none_control - starvation_exists_H16")
dexp <- exprDiff(l_da, contrasts=contrasts)

pcutoff <- 0.1; type<-"padj"
gsetgroup <- "biological_process"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "molecular_function"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "cellular_component"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "kegg"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)

pcutoff <- 0.05; type="pval"
gsetgroup <- "biological_process"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "molecular_function"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "cellular_component"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "kegg"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)








gsetgroup <- "rfbr_orthologs"; gsets <- loadSets(group=gsetgroup); gsets <- gsets[gsets[,3]!="missing",]; plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "aging_diseases"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)
gsetgroup <- "aging_pathways"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, gsets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcutoff)




contrasts = c("control_f - dioxin_f", "control_f - formaldehyde_f", "control_f - radiation_f", "control_f - toluene_f",
              "control_m - dioxin_m", "control_m - formaldehyde_m", "control_m - radiation_m", "control_m - toluene_m",
              "temperature_C25_control - temperature_C4_H24", "temperature_C25_control - temperature_C0_H24", "temperature_C25_control - temperature_C-4_H24",
              "fungus_none_control - fungus_min_H24", "fungus_none_control - fungus_max_H24",
              "radiation_G0_control - radiation_G200_H48", "radiation_G0_control - radiation_G500_H48", "radiation_G0_control - radiation_G1200_H48",
              "starvation_none_control - starvation_exists_H16")

# express differences
dexp <- exprDiff(l_da, contrasts=contrasts)

df_table <- data.frame("FBtr"=dexp$name,
                       "pval_dioxin_f" = dexp$pval[[1]],
                       "pval_formaldehyde_f" = dexp$pval[[2]],
                       "pval_radiation_f" = dexp$pval[[3]],
                       "pval_toluene_f" = dexp$pval[[4]],
                       "pval_dioxin_m" = dexp$pval[[5]],
                       "pval_formaldehyde_m" = dexp$pval[[6]],
                       "pval_radiation_m" = dexp$pval[[7]],
                       "pval_toluene_m" = dexp$pval[[8]],
                       "pval_temperature_C4_H24" = dexp$pval[[9]],
                       "pval_temperature_C0_H24" = dexp$pval[[10]],
                       "pval_temperature_C-4_H24" = dexp$pval[[11]],
                       "pval_fungus_min_H24" = dexp$pval[[12]],
                       "pval_fungus_max_H24" = dexp$pval[[13]],
                       "pval_radiation_G200_H48" = dexp$pval[[14]],
                       "pval_radiation_G500_H48" = dexp$pval[[15]],
                       "pval_radiation_G1200_H48" = dexp$pval[[16]],
                       "pval_starvation_exists_H16" = dexp$pval[[17]],
                       "padj_dioxin_f" = dexp$padj[[1]],
                       "padj_formaldehyde_f" = dexp$padj[[2]],
                       "padj_radiation_f" = dexp$padj[[3]],
                       "padj_toluene_f" = dexp$padj[[4]],
                       "padj_dioxin_m" = dexp$padj[[5]],
                       "padj_formaldehyde_m" = dexp$padj[[6]],
                       "padj_radiation_m" = dexp$padj[[7]],
                       "padj_toluene_m" = dexp$padj[[8]],
                       "padj_temperature_C4_H24" = dexp$padj[[9]],
                       "padj_temperature_C0_H24" = dexp$padj[[10]],
                       "padj_temperature_C-4_H24" = dexp$padj[[11]],
                       "padj_fungus_min_H24" = dexp$padj[[12]],
                       "padj_fungus_max_H24" = dexp$padj[[13]],
                       "padj_radiation_G200_H48" = dexp$padj[[14]],
                       "padj_radiation_G500_H48" = dexp$padj[[15]],
                       "padj_radiation_G1200_H48" = dexp$padj[[16]],
                       "padj_starvation_exists_H16" = dexp$padj[[17]],
                       "lgfc_dioxin_f" = dexp$lgFC[[1]],
                       "lgfc_formaldehyde_f" = dexp$lgFC[[2]],
                       "lgfc_radiation_f" = dexp$lgFC[[3]],
                       "lgfc_toluene_f" = dexp$lgFC[[4]],
                       "lgfc_dioxin_m" = dexp$lgFC[[5]],
                       "lgfc_formaldehyde_m" = dexp$lgFC[[6]],
                       "lgfc_radiation_m" = dexp$lgFC[[7]],
                       "lgfc_toluene_m" = dexp$lgFC[[8]],
                       "lgfc_temperature_C4_H24" = dexp$lgFC[[9]],
                       "lgfc_temperature_C0_H24" = dexp$lgFC[[10]],
                       "lgfc_temperature_C-4_H24" = dexp$lgFC[[11]],
                       "lgfc_fungus_min_H24" = dexp$lgFC[[12]],
                       "lgfc_fungus_max_H24" = dexp$lgFC[[13]],
                       "lgfc_radiation_G200_H48" = dexp$lgFC[[14]],
                       "lgfc_radiation_G500_H48" = dexp$lgFC[[15]],
                       "lgfc_radiation_G1200_H48" = dexp$lgFC[[16]],
                       "lgfc_starvation_exists_H16" = dexp$lgFC[[17]]
)
df_table$FBtr <- as.character(df_table$FBtr)
write.table(df_table,file=paste("./Figures/datatable.csv",sep=""),sep=";",na="",row.names=FALSE,col.names=TRUE)

v_rname <- pv <- av <- vector()
for(i in 1:dim(dexp$pval)[2]) {
    df_table[is.na(df_table[,(1+i)]),(1+i)] <- 1
    df_table[is.na(df_table[,(11+i)]),(11+i)] <- 1
    df_table[is.na(df_table[,(19+i)]),(19+i)] <- 0
    
    #     pv <- df_table[which(df_table[,(1+i)] <=0.5 & df_table[,(19+i)] != 0),1]
    #     av <- df_table[which(df_table[,(11+i)]<=0.1 & df_table[,(19+i)] != 0),1]
    #     v_rname <- c(v_rname,pv,av)
}

trnscrNames <- mappedkeys(org.Dm.egENSEMBLTRANS2EG)
cols <- c("ENTREZID", "SYMBOL")
anno1 <- select(org.Dm.eg.db, keys = trnscrNames, columns = cols, keytype = "ENSEMBLTRANS")

ensembl = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
anno2 <- getBM(attributes=c("flybase_transcript_id", "entrezgene", "flybasename_gene"), mart=ensembl)

length(anno1[is.na(anno1[,1]),1]); length(anno1[is.na(anno1[,2]),2]); length(anno1[is.na(anno1[,3]),3])
length(anno2[is.na(anno2[,1]),1]); length(anno2[is.na(anno2[,2]),2]); length(anno2[is.na(anno2[,3]),3])
length(anno[is.na(anno[,1]),1]); length(anno[is.na(anno[,2]),2]); length(anno[is.na(anno[,3]),3])

length(anno1[,1]); length(anno2[,1]); length(anno[,1])

anno <- merge(anno1,anno2,by.x="FBtr",by.y="FBtr")

length(anno[is.na(anno[,1]),1]); length(anno[is.na(anno[,2]),2]); length(anno[is.na(anno[,3]),3]); length(anno[is.na(anno[,4]),4]); length(anno[is.na(anno[,5]),5])




df_table2 <- merge(anno2,df_table,by="FBtr")
write.table(df_table2,file=paste("./Figures/datatable.csv",sep=""),sep=";",na="",row.names=FALSE,col.names=TRUE)




colnames(anno1) <- colnames(anno2) <- c("FBtr", "ENTREZID", "GNAME")
anno <- rbind(anno1, anno2)
anno <- anno[!duplicated(anno),]











