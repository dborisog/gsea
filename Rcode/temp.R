
#########################################################
#########################################################
#########################################################
#########################################################

dmT <- readRDS("./Data/dmT_lg.rds")
dmT
rm(trnscrNames,trnscrNames2,ensembl)
trnscrNames <- mappedkeys(org.Dm.egENSEMBLTRANS2EG)
ensembl = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
trnscrNames2 <- getBM(attributes=c("flybase_transcript_id"), mart=ensembl)
trnscrNames <- unique(c(as.character(trnscrNames), as.character(trnscrNames2)))

# for either of GO: "biological_process" "molecular_function" "cellular_component" 

x <- getBM(attributes=c("flybase_transcript_id", "entrezgene","flybasename_gene"),values=trnscrNames, mart=ensembl)
dmT2 <- dmT[,c(1,4)]
df_orth <- merge(x, dmT2, by.x="entrezgene",by.y="ENTREZID",all.x=TRUE)
df_orth[is.na(df_orth[,4]),4] <- "missing"
df_orth <- df_orth[!is.na(df_orth[,1]),]
df_orth <- df_orth[,c(2,1,4,3)]
colnames(df_orth) <- c("FBtranscriptID", "ENTREZID","GSET","GNAME")

saveRDS(df_orth,paste("./Data/dmT2_lg.rds",sep=""))



#########################################################
#########################################################
#########################################################
#########################################################

library("AnnotationDbi")
library("GenomicFeatures")
library("DESeq")
library("org.Dm.eg.db")
library("GO.db")
library(data.table) 

source("http://bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")


#########################################################
#########################################################
#########################################################
#########################################################



setwd("C:/Users/user/Google Drive/Work/Semantic research framework/RFBR 2014, Aging gene ontology on invertibrates/GSEA/GSEA_vc01")

countTable <- read.csv("./Data/all.txt", header = TRUE, row.names=1, sep = "",  dec = ",", fill = TRUE)
condition <- factor(c("Control_m","Control_m","Control_fm","Control_fm","Irradiation_m","Irradiation_m","Irradiation_fm","Irradiation_fm","Dioxin_m","Dioxin_m","Dioxin_fm","Dioxin_fm","formaldehyde_m","formaldehyde_m","formaldehyde_fm","formaldehyde_fm","toluene_m","toluene_m","toluene_fm","toluene_fm"))
cds <- newCountDataSet(countTable, condition)
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds )
res <- nbinomTest(cds, "Control_fm","toluene_fm")
topTable = res[order(res$pval),]
# 
# discover & select
keytypes(org.Dm.eg.db)
columns(org.Dm.eg.db)
head(keys(org.Dm.eg.db, keytype="ENSEMBLTRANS"))
fbids <- topTable$id[1:5]
cols <- c("ENTREZID", "GO")
anno <- select(org.Dm.eg.db, keys = fbids, columns = cols, keytype = "ENSEMBLTRANS")
head(anno)
merging top table and annotation
fbids <- topTable$id
# fbids <- res$id
cols <- c("ENTREZID")
anno <- select(x = org.Dm.eg.db, keys = fbids, columns = cols, keytype = "ENSEMBLTRANS")
resAnno <- merge(res, anno, by.x="id", by.y = "ENSEMBLTRANS", all.x = TRUE)

# prepare annotation
# though this line would not work because of non-unique transcript id
rownames(resAnno) <- resAnno$id

    if (fitMethod == "nbinom") {
        dataForLimma <- arrayData$dataCount
        res <- nbinomTest(cds, "Control_fm","toluene_fm")
        
        countTable <- read.csv("./Data/all.txt", header = TRUE, row.names=1, sep = "",  dec = ",", fill = TRUE)
        condition <- factor(c("Control_m","Control_m","Control_fm","Control_fm","Irradiation_m","Irradiation_m","Irradiation_fm","Irradiation_fm","Dioxin_m","Dioxin_m","Dioxin_fm","Dioxin_fm","formaldehyde_m","formaldehyde_m","formaldehyde_fm","formaldehyde_fm","toluene_m","toluene_m","toluene_fm","toluene_fm"))
        cds <- newCountDataSet(countTable[1:100,], condition)
        cds <- estimateSizeFactors( cds )
        cds <- estimateDispersions( cds )
        res <- nbinomTest(cds, "Control_fm","toluene_fm")
        topTable = res[order(res$pval),]
        
    }
	  
  
fct <- rownames(subset(extractFactors(l_da)$factors, factors == "Control_fm" | factors == "Irradiation_fm"))


#########################################################
# Venn diagram
#########################################################

########
# draw a venn diagram for these categories
# > unique(dmG_lg$Influence)
# [1] "dmG_anti" "dmG_pro"  "dmG_none"
# > unique(dmOce_lg$Influence)
# [1] "dmOce_anti"        "dmOce_pro"         "dmOce_unclear"     "dmOce_unannotated"
# > unique(dmOmm_lg$Influence)
# [1] "dmOmm_pro"  "dmOmm_anti"
dmG_anti <- dm_lg$Influence == "dmG_anti"                       # 1
dmG_pro <- dm_lg$Influence == "dmG_pro"                         # 2
dmG_none <- dm_lg$Influence == "dmG_none"                       # 3
dmOce_anti <- dm_lg$Influence == "dmOce_anti"                   # 4
dmOce_pro <- dm_lg$Influence == "dmOce_pro"                     # 5
dmOce_unclear <- dm_lg$Influence == "dmOce_unclear"             # 6
dmOce_unannotated <- dm_lg$Influence == "dmOce_unannotated"     # 7
dmOmm_pro <- dm_lg$Influence == "dmOmm_pro"                     # 8
dmOmm_anti <- dm_lg$Influence == "dmOmm_anti"                   # 9

# prepare massive Venn
clms <- cbind(dmG_anti, dmG_pro, dmG_none, 
              dmOce_anti, dmOce_pro, dmOce_unclear, dmOce_unannotated, 
              dmOmm_pro, dmOmm_anti)
rownames(clms) <- dm_lg$ENTREZID

cPro <- vennCounts(clms[,c(2,5,8)])
vennDiagram(cPro)

cAnti <- vennCounts(clms[,c(1,4,9)])
vennDiagram(cAnti)

cNem <- vennCounts(clms[,c(1,2,4,5,6)])
vennDiagram(cNem)

cMouse <- vennCounts(clms[,c(1,2,3,8,9)])
vennDiagram(cMouse)

cNM <- vennCounts(clms[,c(4,5,8,9)])
vennDiagram(cNM)

cNM <- vennCounts(clms[,c(4,5,8,9)])
vennDiagram(cNM)

#########################################################
#########################################################
#########################################################
# ANNOTATION
#########################################################

setwd("C:/Users/user/Google Drive/Work/Semantic research framework/RFBR 2014, Aging gene ontology on invertibrates/GSEA/GSEA_vc01")
geneadm <- read.delim(file="./Data/genage_models_export_v2014-12-02.tsv", header=TRUE)

geneadm[,c(3,4,5,6,7,8,9,10,11,12,15,16)] <- data.frame(lapply(geneadm[,c(3,4,5,6,7,8,9,10,11,12,15,16)], as.character), stringsAsFactors=FALSE)

#TRGT SRC     COLUMN
# 01  --      DB
# 02  01      DB Object ID
# 03  03      DB Object Symbol
# 04  --      Qualifier
# 05  02      ENTREZID
# 06  --      DB:Reference
# 07  --      Evidence Code
# 08  --      With (or) From
# 09  04      DB Object Name
# 10  09      DB Object Synonym
# 11  --      DB Object Type
# 12  --      Taxon
# 13  --      Date
# 14  --      Annotation Extension
# 15  --      Gene Product Form ID
# 16  05      Model organism
# 17  --      Tissue
# 18  12      Influence


gatorfbr <- data.frame(DB = "GeneAge",
                       DB.Object.ID = geneadm[,1],
                       DB.Object.Symbol = geneadm[,3],
                       ENTREZID = geneadm[,2],
                       DB.Reference = "",
                       Evidence.Code = "",
                       With.or.From = "",
                       DB.Object.Name = geneadm[,4],
                       DB.Object.Synonym = geneadm[,9],
                       DB.Object.Type = "",
                       Date = "",
                       Annotation.Extension = "",
                       Gene.Product.Form.ID = "",
                       Taxon = geneadm[,5],
                       Tissue = "",
                       Influence = geneadm[,12]
                       )

write.table(gatorfbr, file="./data_from_geneage.csv",
            quote=FALSE, sep="\t",
            row.names=FALSE, col.names=TRUE)


gatorfbr[which(gatorfbr$ENTREZID %in% c(32001, 53581, 31845, 41709, 248194, 31251, 31220, 43904, 40336)),]

gatorfbr[which(gatorfbr$DB.Object.ID %in% c(22, 294, 848, 286, 294, 1756, 1843)),]


#########################################################
#########################################################
#########################################################
#########################################################


getTranscriptGene <- function(data) {
    #### transcript to gene
    #### from org.Dm.eg.db
    fbids <- dexp$resTable[[1]]$ProbesetID

    cols <- c("ENTREZID")
    anno <- select(org.Dm.eg.db, keys = fbids, columns = cols, keytype = "ENSEMBLTRANS")
    
    #### transcript to gene
    #### from biomaRt
    ensembl = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
    flyentr <- getBM(attributes=c("flybase_transcript_id", "entrezgene"),
                     filters = "flybase_transcript_id",
                     values = fbids, mart = ensembl)
    
    #### transcript to gene
    #### combine org.Dm.eg.db with biomaRt
    colnames(flyentr) <- c("ENSEMBLTRANS", "ENTREZID")
    full <- rbind(anno,flyentr)
    full <- unique(full)
    full <- full[complete.cases(full),]
    
    return(full)
}

#############
#############
#############
rm("loadGSC")
source("./R code/loadGSC.R")
# "biological_process" "molecular_function" "cellular_component" "rfbr_orthologs"  "kegg" "aging_diseases" "aging_pathways"
goSubgroup <- "kegg"
myGsc <- loadGSC(dexp, group=goSubgroup)
#############
#############
#############

loadGSC <- function(gdata, group) {
    fbids <- as.character(gdata$resTable[[1]][,1])
    
    if (file.exists(paste("./Data/dmelanogaster_getBM","_",group,".rds",sep=""))) {
        res <- readRDS(paste("./Data/dmelanogaster_getBM","_",group,".rds",sep=""))
    } else {
        
        if (group=="rfbr_orthologs") {
            
            fbids <- as.data.frame(fbids)
            colnames(fbids) <- "ENSEMBLTRANS"
            
            dmT <- readRDS("./Data/dmT_lg.rds")
            dmTm <- merge(fbids,dmT,by="ENSEMBLTRANS",all.x=TRUE)
            
            dmTm$Influence[is.na(dmTm$Influence)] <- "missing"
            genes2genesets <- unique(dmTm[which(complete.cases(dmTm[,c(1,3)])),c(1,3)])
            
            
        } else if(group=="kegg") {
            cols <- c("ENTREZID", "PATH")
            anno <- select(org.Dm.eg.db, keys = fbids, columns = cols, keytype = "ENSEMBLTRANS")
            
            eid <- anno[complete.cases(anno[,3]),c(2,3)]
            # KEGGPATHID2NAME  An annotation data object that maps KEGG pathway identifiers to KEGG pathway names
            yy <- as.data.frame(KEGGPATHID2NAME)
            dm_keggname <- yy[which(yy$path_id %in% eid[,2]),]
            dmT_kegg <- merge(anno,dm_keggname, by.x="PATH", by.y="path_id", all.x=TRUE)
            dm_kegg <- dmT_kegg[,c(2,4)]
            genes2genesets <- unique(dm_kegg[which(complete.cases(dm_kegg)),])
            
        } else if(group=="aging_diseases") {
            ad.name <- c("Type II diabetes mellitus", "Colorectal cancer", "Alzheimer's disease", "Amyotrophic lateral sclerosis (ALS)", 
                         "Endometrial cancer", "Huntington's disease", "Hypertrophic cardiomyopathy (HCM)", "Long-term depression", 
                         "Non-small cell lung cancer", "Pancreatic cancer", "Parkinson's disease", "Pathways in cancer", "Prostate cancer", 
                         "Non-alcoholic fatty liver disease (NAFLD)")
            
            genes2genesets <- useKEGGDrivenOntology(ad.name)
        } else if(group=="aging_pathways") {
            ap.name <- c("Apoptosis", "Cell cycle", "Circadian rhythm - mammal", "Glutathione metabolism", "Insulin signaling pathway",
                         "MAPK signaling pathway", "mTOR signaling pathway", "p53 signaling pathway", "PPAR signaling pathway",
                         "Regulation of autophagy", "Toll-like receptor signaling pathway",
                         "AMPK signaling", "FoxO signaling", "HIF-1 signaling", "PI3K-Akt signaling", "Ras signaling pathway")
            
            genes2genesets <- useKEGGDrivenOntology(ap.name)
        } else {
            ensembl = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
            x <- getBM(attributes=c("flybase_transcript_id", "name_1006", "namespace_1003"),
                       # attributes=c("flybase_transcript_id", "name_1006")
                       values=fbids, 
                       mart=ensembl)
            genes2genesets <- unique(x[which(complete.cases(x) & x[,3]==group),1:2])
        }
        
        tmp <- try(gsc <- as.data.frame(genes2genesets, stringsAsFactors=FALSE), silent=TRUE)
        if(class(tmp) == "try-error") {
            stop("argument gdata could not be converted into a data.frame")
        }
        
        
        
        # Get rid of factors:
        for(i in 1:ncol(gsc)) {
            gsc[,i] <- as.character(gsc[,i])
        }
        
        # Check gsc for two columns:
        if(ncol(gsc)!=2) {
            stop("argument gdata has to contain exactly two columns")  
        }
        
        # Remove redundant rows:
        tmp <- nrow(gsc)
        gsc <- unique(gsc)
        #info$redundantGS <- tmp - nrow(gsc)
        
        # Convert to list object:
        geneSets <- unique(gsc[,2])
        gscList <- list()
        for(iGeneSet in 1:length(geneSets)) {
            gscList[[iGeneSet]] <- gsc[gsc[,2] == geneSets[iGeneSet],1]
        }
        names(gscList) <- geneSets
        gsc <- gscList
        

        
        #********************************
        # Return values:
        #********************************
        
        res <- list(gsc,addInfo)
        names(res) <- c("gsc","addInfo")
        class(res) <- "GSC"
        
        saveRDS(res,paste("./Data/dmelanogaster_getBM","_",group,".rds",sep=""))
        # listAttributes(mart)[c(27:32,37:50,73:74,91:100,123:130,414:419,1004:1009,1017:1024,1079:1085),]
        # listAttributes(mart)[c(27:32,37:50,73:74,91:100,123:130,414:419,1004:1009,1017:1024,1079:1085),]
    }
    
    
    return(res)
    
}

# GSEA: GO-BP, GO-MF, GO-CC; KEGG; RFBR1; RFBR2; aging diseases; aging pathways
# + GO-BP, GO-MF, GO-CC; 
# - KEGG; 
# + RFBR1; 
# - RFBR2; 
# - aging diseases; 
# - aging pathways


#########################################################
#########################################################
# KEGG
#########################################################
source("http://bioconductor.org/biocLite.R")
biocLite("KEGG.db")
library("KEGG.db")


fbids <- dmT$ENTREZID

# fbids <- as.character(gdata$resTable[[1]][,1])

# if (group=="rfbr_orthologs") {
#     
#     fbids <- as.data.frame(fbids)
#     colnames(fbids) <- "ENSEMBLTRANS"
#     
#     dmT <- readRDS("./Data/dmT_lg.rds")
#     dmTm <- merge(fbids,dmT,by="ENSEMBLTRANS",all.x=TRUE)
#     
#     dmTm$Influence[is.na(dmTm$Influence)] <- "missing"
#     genes2genesets <- unique(dmTm[which(complete.cases(dmTm[,c(1,3)])),c(1,3)])
#     
# } else if(group=="kegg") {
    cols <- c("ENTREZID", "PATH")
    anno <- select(org.Dm.eg.db, keys = fbids, columns = cols, keytype = "ENSEMBLTRANS")
    
    eid <- anno[complete.cases(anno[,3]),c(2,3)]
    # KEGGPATHID2NAME  An annotation data object that maps KEGG pathway identifiers to KEGG pathway names
    yy <- as.data.frame(KEGGPATHID2NAME)
    dm_keggname <- yy[which(yy$path_id %in% eid[,2]),]
    dmT_kegg <- merge(anno,dm_keggname, by.x="PATH", by.y="path_id", all.x=TRUE)
    dm_kegg <- dmT_kegg[,c(2,4)]
    genes2genesets <- unique(dm_kegg[which(complete.cases(dm_kegg)),])
# } else {
#     ensembl = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
#     x <- getBM(attributes=c("flybase_transcript_id", "name_1006", "namespace_1003"),
#                # attributes=c("flybase_transcript_id", "name_1006")
#                values=fbids, 
#                mart=ensembl)
#     genes2genesets <- unique(x[which(complete.cases(x) & x[,3]==group),1:2])
# }


#########################################################
#########################################################
# RFBR2;
#########################################################

#########################################################
#########################################################
# aging diseases;
#########################################################
# + 1. have a vector of KEGG pathways associated to human aging disesases
# + 2. transform to fly orthologs
# + 3. ling flybase transcript ids to fly these orthologs

# 1. vektor of KEGG pathways
ad.name <- c("Type II diabetes mellitus", "Colorectal cancer", "Alzheimer's disease", "Amyotrophic lateral sclerosis (ALS)", 
             "Endometrial cancer", "Huntington's disease", "Hypertrophic cardiomyopathy (HCM)", "Long-term depression", 
             "Non-small cell lung cancer", "Pancreatic cancer", "Parkinson's disease", "Pathways in cancer", "Prostate cancer", 
             "Non-alcoholic fatty liver disease (NAFLD)")

genes2genesets <- useKEGGDrivenOntology(ad.name)
    

#########################################################
#########################################################
# aging pathways
#########################################################

# 1. vektor of KEGG pathways
ap.name <- c("Apoptosis", "Cell cycle", "Circadian rhythm - mammal", "Glutathione metabolism", "Insulin signaling pathway",
             "MAPK signaling pathway", "mTOR signaling pathway", "p53 signaling pathway", "PPAR signaling pathway",
             "Regulation of autophagy", "Toll-like receptor signaling pathway",
             "AMPK signaling", "FoxO signaling", "HIF-1 signaling", "PI3K-Akt signaling", "Ras signaling pathway")

genes2genesets <- useKEGGDrivenOntology(ap.name)


#########################################################
#########################################################
#########################################################
#########################################################

#########################################################
#########################################################
# visualization: PCA
#########################################################

###################################
# all normalized data #############
rm("runQC")
source("./R code/runQC.R")
runQC(l_da, hist=FALSE, boxplot=FALSE, pca=TRUE, verbose=FALSE)



###################################
# significant transcripts #########

contrasts = c("Control_fm - Irradiation_fm", 
              "Control_fm - toluene_fm",
              "Control_fm - Dioxin_fm",
              "Control_fm - formaldehyde_fm",
              "Control_m - Irradiation_m", 
              "Control_m - toluene_m",
              "Control_m - Dioxin_m",
              "Control_m - formaldehyde_m")

rm(rownm)
for (i in 1:length(contrasts)) {
    rownm <- c(rownm,dexp$resTable[[i]][which(dexp$resTable[[i]]$P.Value <= 0.05),1])
}
rownm <- sort(unique(rownm))

dataForPCA <- l_da$dataNorm[rownames(l_da$dataNorm) %in% rownm,]
# Centralize data
dataForPCA <- dataForPCA - rowMeans(dataForPCA)

# PCA
dataPrcomp <- prcomp(dataForPCA)
pc1 <- dataPrcomp$rotation[,1]
pc2 <- dataPrcomp$rotation[,2]
clr <- 1:length(pc1)
     
# PCA plot
dev.new()
layout(matrix(c(1,2), nrow=2, ncol=1), heights=c(7,1), widths=c(1,1))
par(mar=c(4,4,2,1))
plot(cbind(pc1,pc2), 
     pch=15, 
     col=1:length(names(pc1)),
     cex=1.5,
     main="PCA plot", 
     xlab="PC1", 
     ylab="PC2")
text(cbind(pc1,pc2), labels=names(pc1), pos=4, cex=0.5)

par(mar=c(0,3,1,1))
plot.new()
title(main="Legend")
legend("bottomleft", legend=summary(dataPrcomp)$importance[2,1:2],ncol=2)

dev.new()
direct.label(xyplot(pc2~pc1,
       cex=1.5,
       pch=15,
       groups=names(pc1)))

#########################################################
#########################################################
# visualization: venn
#########################################################
library("limma")
library("AnnotationDbi")
library("biomaRt")
library("org.Dm.eg.db")

library("dplyr")
library("data.table")

###################################
# simple venn #####################
df.tr2data <- getTranscriptGene(dexp)

contrasts=c("Control_fm - Irradiation_fm", 
            "Control_fm - toluene_fm",
            "Control_fm - Dioxin_fm",
            "Control_fm - formaldehyde_fm")


contrasts=c("Control_m - Irradiation_m", 
            "Control_m - toluene_m",
            "Control_m - Dioxin_m",
            "Control_m - formaldehyde_m")

rm("vennPlot")
source("./R code/vennPlot.R")
vennPlot(dexp, df.tr2data, contrasts, pvalue=0.05, adj.pvalue=0.1)



#########################################################
#########################################################
# visualization: heatmap
#########################################################
###################################
# by No of genes in the group #####
library("AnnotationDbi")
library("GenomicFeatures")
library("DESeq")
library("org.Dm.eg.db")
library("GO.db")
library(data.table) 
library(gplots) 
###################################
# by pValue #######################
contrast <- "Control_fm - toluene_fm"


contrasts = c("Control_fm - Irradiation_fm", 
              "Control_fm - toluene_fm",
              "Control_fm - Dioxin_fm",
              "Control_fm - formaldehyde_fm",
              "Control_m - Irradiation_m", 
              "Control_m - toluene_m",
              "Control_m - Dioxin_m",
              "Control_m - formaldehyde_m")

contrasts = c("Control_fm - Irradiation_fm", 
              "Control_fm - toluene_fm")


# "biological_process" "molecular_function" "cellular_component" "rfbr_orthologs"  "kegg" "aging_diseases" "aging_pathways"
goSubgroup <- "biological_process"
myGsc <- loadGSC(dexp, group=goSubgroup)
GSAsign <- 0.1

#########################
# not regulated
rm("heatmapPlot")
source("./R code/heatmapPlot.R")
heatmapPlot(data = dexp, contrasts = contrasts, gsc = myGsc, GSAsign = 0.1)


#########################################################
#########################################################
#########################################################
#########################################################

# GSEA: transform transcripts to genes, genes -- to genesets

fb.entrez <- unlist(as.list(org.Dm.egALIAS2EG))

