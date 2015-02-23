

# p0301 <- read.table(file="./Data/p03_01_data.csv",sep="",header=TRUE,fill=TRUE)
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
    source("./Rcode/plotHeatmap.R")
}


setwd(path)
loadPackages(packages,src)
loadSources()
conc = function(one,two) {paste(one,two,sep="")} #shorcut function for concatenation

#############

# l_da <- loadData(datadir = dataFolder, countfile = "p03_data.csv", setup="setup.txt")
# saveRDS(l_da,"./Data/datafile_p03.rds")

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

contrasts = c("Control_m - Irradiation_m", 
              "Control_m - toluene_m",
              "Control_m - Dioxin_m",
              "Control_m - formaldehyde_m",
              "temperature_C25_control - temperature_C4_H24",
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

# vennPlot
rm(plotVenn); source("./Rcode/plotVenn.R")

contrasts = c("temperature_C25_control - temperature_C4_H24","temperature_C25_control - temperature_C0_H24","temperature_C25_control - temperature_C-4_H24")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(dexp = dexp, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

contrasts = c("fungus_none_control - fungus_min_H24", "fungus_none_control - fungus_max_H24")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(dexp = dexp, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

contrasts = c("Control_m - Irradiation_m", "radiation_G0_control - radiation_G200_H48","radiation_G0_control - radiation_G500_H48","radiation_G0_control - radiation_G1200_H48")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(dexp = dexp, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)

contrasts = c("starvation_none_control - starvation_exists_H16")
dexp <- exprDiff(l_da, contrasts=contrasts)
plotVenn(dexp = dexp, contrasts = contrasts, pvalue=0.05, adj.pvalue=0.1)


# enrich sets
contrasts = c("temperature_C25_control - temperature_C4_H24","temperature_C25_control - temperature_C0_H24","temperature_C25_control - temperature_C-4_H24",
              "fungus_none_control - fungus_min_H24","fungus_none_control - fungus_max_H24",
              "radiation_G0_control - radiation_G200_H48","radiation_G0_control - radiation_G500_H48","radiation_G0_control - radiation_G1200_H48",
              "starvation_none_control - starvation_exists_H16")
dexp <- exprDiff(l_da, contrasts=contrasts)



gsetgroup <- "biological_process"; gsets <- loadSets(group=gsetgroup)
gsetgroup <- "molecular_function"; gsets <- loadSets(group=gsetgroup)
gsetgroup <- "cellular_component"; gsets <- loadSets(group=gsetgroup)
gsetgroup <- "rfbr_orthologs"; gsets <- loadSets(group=gsetgroup)
gsetgroup <- "kegg"; gsets <- loadSets(group=gsetgroup)
gsetgroup <- "aging_diseases"; gsets <- loadSets(group=gsetgroup)
gsetgroup <- "aging_pathways"; gsets <- loadSets(group=gsetgroup)

pcut <- 0.1; type<-"padj"
pcut <- 0.05; type="pval"


for (contrast in contrasts) {
    if (file.exists(paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcut,".rds",sep=""))) {
        gsah<- readRDS(paste("./Data/GSA_hyper_",contrast,"_",gsetgroup,"_",pcut,".rds",sep=""))
    } else {
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
                    type = type, # either "pval" or "padj"
                    gsetgroup = gsetgroup) 
    }, error=function(e){cat("No network diagram for ", contrast, ", ", gsetgroup,"\n")})
}



plotHeatmap(dexp = dexp, contrasts = contrasts, sets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcut)


pcut <- 0.1; type<-"padj"
gsetgroup <- "biological_process"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, sets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcut)
gsetgroup <- "molecular_function"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, sets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcut)
gsetgroup <- "cellular_component"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, sets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcut)
gsetgroup <- "kegg"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, sets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcut)

pcut <- 0.05; type="pval"
gsetgroup <- "biological_process"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, sets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcut)
gsetgroup <- "molecular_function"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, sets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcut)
gsetgroup <- "cellular_component"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, sets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcut)
gsetgroup <- "kegg"; gsets <- loadSets(group=gsetgroup); plotHeatmap(dexp = dexp, contrasts = contrasts, sets = gsets, gsetgroup = gsetgroup, type = type, pcutoff = pcut)





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

df_table <- data.frame("FBtr"=dexp$name,
                       "pval_temperature_C4_H24" = dexp$pval[[1]],
                       "pval_temperature_C0_H24" = dexp$pval[[2]],
                       "pval_temperature_C-4_H24" = dexp$pval[[3]],
                       "pval_fungus_min_H24" = dexp$pval[[4]],
                       "pval_fungus_max_H24" = dexp$pval[[5]],
                       "pval_radiation_G200_H48" = dexp$pval[[6]],
                       "pval_radiation_G500_H48" = dexp$pval[[7]],
                       "pval_radiation_G1200_H48" = dexp$pval[[8]],
                       "pval_starvation_exists_H16" = dexp$pval[[9]],
                       "padj_temperature_C4_H24" = dexp$padj[[1]],
                       "padj_temperature_C0_H24" = dexp$padj[[2]],
                       "padj_temperature_C-4_H24" = dexp$padj[[3]],
                       "padj_fungus_min_H24" = dexp$padj[[4]],
                       "padj_fungus_max_H24" = dexp$padj[[5]],
                       "padj_radiation_G200_H48" = dexp$padj[[6]],
                       "padj_radiation_G500_H48" = dexp$padj[[7]],
                       "padj_radiation_G1200_H48" = dexp$padj[[8]],
                       "padj_starvation_exists_H16" = dexp$padj[[9]],
                       "lgfc_temperature_C4_H24" = dexp$lgFC[[1]],
                       "lgfc_temperature_C0_H24" = dexp$lgFC[[2]],
                       "lgfc_temperature_C-4_H24" = dexp$lgFC[[3]],
                       "lgfc_fungus_min_H24" = dexp$lgFC[[4]],
                       "lgfc_fungus_max_H24" = dexp$lgFC[[5]],
                       "lgfc_radiation_G200_H48" = dexp$lgFC[[6]],
                       "lgfc_radiation_G500_H48" = dexp$lgFC[[7]],
                       "lgfc_radiation_G1200_H48" = dexp$lgFC[[8]],
                       "lgfc_starvation_exists_H16" = dexp$lgFC[[9]]
)
df_table$FBtr <- as.character(df_table$FBtr)

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











