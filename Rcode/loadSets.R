loadSets <- function(dexp, group) {
    trnscrNames <- mappedkeys(org.Dm.egENSEMBLTRANS2EG)
    
    if (file.exists(paste("./Data/dmelanogaster_getBM","_",group,".rds",sep=""))) {
        gsc <- readRDS(paste("./Data/dmelanogaster_getBM","_",group,".rds",sep=""))
    } else {
        
        if (group=="rfbr_orthologs") {
            dmT <- readRDS("./Data/dmT_lg.rds") #see prepareGSets.R for more info
            dmT[is.na(dmT[,1]),1] <- "missing"; dmT[is.na(dmT[,2]),2] <- "missing"; dmT[is.na(dmT[,3]),3] <- "missing"; dmT[is.na(dmT[,4]),4] <- "missing"
            FBtr_entrezid_gset <- unique(dmT[,c(2,1,4)])
            
        } else if(group=="kegg") {
            cols <- c("ENTREZID", "PATH")
            anno <- select(org.Dm.eg.db, keys = trnscrNames, columns = cols, keytype = "ENSEMBLTRANS")
            eid <- anno[complete.cases(anno),]
            # KEGGPATHID2NAME  An annotation data object that maps KEGG pathway identifiers to KEGG pathway names
            yy <- as.data.frame(KEGGPATHID2NAME)
            dm_keggname <- yy[which(yy$path_id %in% eid[,3]),]
            dmT_kegg <- merge(anno,dm_keggname, by.x="PATH", by.y="path_id", all.x=TRUE)
            dm_kegg <- dmT_kegg[,c(2,3,4)]
            dm_kegg[is.na(dm_kegg[,2]),2] <- "missing"; dm_kegg[is.na(dm_kegg[,3]),3] <- "missing"
            FBtr_entrezid_gset <- unique(dm_kegg[which(complete.cases(dm_kegg)),])
            
        } else if(group=="aging_diseases") {
            v.aging <- c("Type II diabetes mellitus", "Colorectal cancer", "Alzheimer's disease", "Amyotrophic lateral sclerosis (ALS)", 
                         "Endometrial cancer", "Huntington's disease", "Hypertrophic cardiomyopathy (HCM)", "Long-term depression", 
                         "Non-small cell lung cancer", "Pancreatic cancer", "Parkinson's disease", "Pathways in cancer", "Prostate cancer", 
                         "Non-alcoholic fatty liver disease (NAFLD)")
            FBtr_entrezid_gset <- useKEGGDrivenOntology(v.aging)
            
        } else if(group=="aging_pathways") {
            v.aging <- c("Apoptosis", "Cell cycle", "Circadian rhythm - mammal", "Glutathione metabolism", "Insulin signaling pathway",
                         "MAPK signaling pathway", "mTOR signaling pathway", "p53 signaling pathway", "PPAR signaling pathway",
                         "Regulation of autophagy", "Toll-like receptor signaling pathway", 
                         "AMPK signaling", "FoxO signaling", "HIF-1 signaling", "PI3K-Akt signaling", "Ras signaling pathway")
            FBtr_entrezid_gset <- useKEGGDrivenOntology(v.aging)
            
        } else {
            # for either of GO: "biological_process" "molecular_function" "cellular_component" 
            ensembl = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
            x <- getBM(attributes=c("flybase_transcript_id", "entrezgene","name_1006", "namespace_1003"),
                       values=trnscrNames, mart=ensembl)
            
            FBtr_entrezid_gset <- unique(x[x[,4]==group,1:3])
        }
        
        gsc <- data.frame()
        tmp <- try(gsc <- as.data.frame(FBtr_entrezid_gset, stringsAsFactors=FALSE), silent=TRUE)
        if(class(tmp) == "try-error") {stop("argument FBtr_entrezid_gset could not be converted into a data.frame")}
        
        # Get rid of factors:
        for(i in 1:ncol(gsc)) {gsc[,i] <- as.character(gsc[,i])}
        
        # Check gsc for two columns:
        if(ncol(gsc)!=3) {stop("argument gsc has to contain exactly three columns")  }
        
        # Remove redundant rows:
        gsc <- unique(gsc)
        colnames(gsc) <- c("FBtranscriptID", "ENTREZID","GSET")
        
        saveRDS(gsc,paste("./Data/dmelanogaster_getBM","_",group,".rds",sep=""))
    }
    return(gsc)
}


