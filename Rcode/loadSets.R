loadSets <- function(group) {
    
    if (file.exists(paste("./Data/dmelanogaster_getBM","_",group,".rds",sep=""))==TRUE) {
        gsc <- readRDS(paste("./Data/dmelanogaster_getBM","_",group,".rds",sep=""))
    } else {
        
        if (group=="rfbr_orthologs") {
            dmT <- readRDS("./Data/dmT_lg.rds") #see prepareGSets.R for more info
            dmT[is.na(dmT[,1]),1] <- "missing"; dmT[is.na(dmT[,2]),2] <- "missing"; dmT[is.na(dmT[,3]),3] <- "missing"; dmT[is.na(dmT[,4]),4] <- "missing"
            FBtr_entrezid_gset <- unique(dmT[,c(2,1,4)])
            
        } else if(group=="kegg") {
            
            rm(trnscrNames,trnscrNames2,ensembl)
            trnscrNames <- mappedkeys(org.Dm.egENSEMBLTRANS2EG)
#             ensembl = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
#             trnscrNames2 <- getBM(attributes=c("flybase_transcript_id"), mart=ensembl)
#             trnscrNames <- unique(c(as.character(trnscrNames), as.character(trnscrNames2)))
            
            cols <- c("ENTREZID", "PATH", "SYMBOL")
            anno <- select(org.Dm.eg.db, keys = trnscrNames, columns = cols, keytype = "ENSEMBLTRANS")
            anno$PATH[is.na(anno$PATH)] <- "missing"

            # KEGGPATHID2NAME  An annotation data object that maps KEGG pathway identifiers to KEGG pathway names
            yy <- as.data.frame(KEGGPATHID2NAME)
            dm_keggname <- yy[which(yy$path_id %in% eid[,3]),]
            dmT_kegg <- merge(anno,dm_keggname, by.x="PATH", by.y="path_id", all.x=TRUE)
            dmT_kegg$path_name[is.na(dmT_kegg$path_name)] <- "missing"

            dm_kegg <- dmT_kegg[,c(2,3,5,4)]
            
            tr_entr_gset_symb <- unique(dm_kegg[which(complete.cases(dm_kegg)),])
            
        } else if(group=="aging_diseases") {
            v.aging <- c("Type II diabetes mellitus", "Colorectal cancer", "Alzheimer's disease", "Amyotrophic lateral sclerosis (ALS)", 
                         "Endometrial cancer", "Huntington's disease", "Hypertrophic cardiomyopathy (HCM)", "Long-term depression", 
                         "Non-small cell lung cancer", "Pancreatic cancer", "Parkinson's disease", "Pathways in cancer", "Prostate cancer", 
                         "Non-alcoholic fatty liver disease (NAFLD)")
            tr_entr_gset_symb <- useKEGGDrivenOntology(v.aging)
            
        } else if(group=="aging_pathways") {
            v.aging <- c("Apoptosis", "Cell cycle", "Circadian rhythm - mammal", "Glutathione metabolism", "Insulin signaling pathway",
                         "MAPK signaling pathway", "mTOR signaling pathway", "p53 signaling pathway", "PPAR signaling pathway",
                         "Regulation of autophagy", "Toll-like receptor signaling pathway", 
                         "AMPK signaling", "FoxO signaling", "HIF-1 signaling", "PI3K-Akt signaling", "Ras signaling pathway")
            tr_entr_gset_symb <- useKEGGDrivenOntology(v.aging)
            
        } else {
            
            rm(trnscrNames,trnscrNames2,ensembl)
            trnscrNames <- mappedkeys(org.Dm.egENSEMBLTRANS2EG)
            ensembl = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
            trnscrNames2 <- getBM(attributes=c("flybase_transcript_id"), mart=ensembl)
            trnscrNames <- unique(c(as.character(trnscrNames), as.character(trnscrNames2)))
            
            # for either of GO: "biological_process" "molecular_function" "cellular_component" 
            
            x <- getBM(attributes=c("flybase_transcript_id", "entrezgene","name_1006", "flybasename_gene", "namespace_1003"),
                       values=trnscrNames, mart=ensembl)
            
            tr_entr_gset_symb <- unique(x[x[,5]==group,1:4])
            tr_entr_gset_symb[is.na(tr_entr_gset_symb[,3]),3] <- "missing"
        }
        
        gsc <- data.frame()
        tmp <- try(gsc <- as.data.frame(tr_entr_gset_symb, stringsAsFactors=FALSE), silent=TRUE)
        if(class(tmp) == "try-error") {stop("argument tr_entr_gset_symb could not be converted into a data.frame")}
        
        # Get rid of factors:
        for(i in 1:ncol(gsc)) {gsc[,i] <- as.character(gsc[,i])}
        
        # Check gsc for two columns:
        if(ncol(gsc)!=4) {stop("argument gsc has to contain exactly three columns")  }
        
        # Remove redundant rows:
        gsc <- unique(gsc)
        colnames(gsc) <- c("FBtranscriptID", "ENTREZID","GSET","GNAME")
        
        saveRDS(gsc,paste("./Data/dmelanogaster_getBM","_",group,".rds",sep=""))
    }
    return(gsc)
}


