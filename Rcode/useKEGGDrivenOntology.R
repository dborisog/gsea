useKEGGDrivenOntology <- function(ad_name) {
    # KEGGPATHNAME2ID     An annotation data object that maps KEGG pathway names to identifiers for the corresponding pathway names used by KEGG
    ## KEGGnames--KEGGID
    aa <- as.data.frame(KEGGPATHNAME2ID)
    ad_kid <- aa[which(aa$path_name %in% ad_name),]
    ad_kid$path_id <- paste("hsa",ad_kid$path_id,sep="")
    
    
    # KEGGPATHID2EXTID    An annotation data object that maps KEGG pathway identifiers to Entrez Gene or Open Reading Frame identifiers
    ## KEGGID--ENTREZID
    bb <- as.data.frame(KEGGPATHID2EXTID)
    ad_eid <- bb[which(bb$pathway_id %in% ad_kid$path_id),]
    colnames(ad_eid) <- c("path_id","ENTREZID")
    
    ## KEGGID--ENTREZID--KEGGname
    ad_id <- merge(ad_eid,ad_kid,by="path_id")
    
    # org.Hs.egSYMBOL     Map between Entrez Gene Identifiers and Gene Symbols
    cc <- as.data.frame(org.Hs.egSYMBOL)
    colnames(cc) <- c("ENTREZID","SYMBOL")
    ad_s <- cc[which(cc$ENTREZID %in% unique(ad_eid$ENTREZID)),]

    ad_all <- merge(ad_id,ad_s,by="ENTREZID")
    
    
    # 2. orthologs
    # gorth               Find orthologs.
    ## alias.number--initial.alias--initial.ensg--ensg.number--target.ensg--target.name--target.description
    ad_orth <- gorth(unique(ad_all$SYMBOL),
                     source_organism = "hsapiens", target_organism = "dmelanogaster",
                     region_query = F,numeric_ns = "", mthreshold = Inf, filter_na = T, df = T)
    ad_sfly <- ad_orth[which(ad_orth$initial.alias %in% unique(ad_all$SYMBOL)),c(2,5,6)] # ????????
    ## initial.alias--target.ensg
    colnames(ad_sfly) <- c("SYMBOL", "FLYBASE","FSYMBOL")
    
    
    # 3. link KEGG to the othologs
    ## SYMBOL--ENTREZID--KEGGID--KEGGname--FLYBASE
    ad_ofly <- merge(ad_all,ad_sfly,by="SYMBOL")
    
    
    ## ENSEMBLTRANS--ENTREZID--KEGGID--FLYBASE
    cols <- c("ENTREZID", "FLYBASE") # 
    tr_orgdm <- mappedkeys(org.Dm.egENSEMBLTRANS2EG)
    anno <- select(org.Dm.eg.db, keys = tr_orgdm, columns = cols, keytype = "ENSEMBLTRANS")
        
    ensembl = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
    anno2 <- getBM(attributes=c("flybase_transcript_id", "entrezgene", "flybase_gene_id"), mart=ensembl)
    
    colnames(anno2) <- c("ENSEMBLTRANS", "ENTREZID", "FLYBASE")
    sets <- rbind(anno, anno2)
    sets <- sets[!duplicated(sets),]

    sets$FLYBASE <- toupper(sets$FLYBASE)
    
    ## FLYBASE--SYMBOL--ENTREZID.x--KEGGID--KEGGname--ENSEMBLTRANS--ENTREZID.y   PATH 
    ad_end <- merge(ad_ofly,sets,by="FLYBASE")
    
    ad_end <- ad_end[,c(7,8,5,6)]
    ad_end <- ad_end[complete.cases(ad_end),]
    # return table: Fly transcript, Fly gene, Human gene group
    
    return(unique(ad_end))
}
# Flybase "egSYMBOL" == biomaRt "flybasename_gene"

# alternative, though keggConv processes 10 entries     
# kglnk <- keggLink("pathway", "hsa")
# kegg_gn_pth <- data.frame(kgene = names(kglnk), kpath = kglnk)
# kegg_gn_pth[,1] <- as.character(kegg_gn_pth[,1])
# entrzid <- keggConv("ncbi-geneid", kegg_gn_pth[,1])

# alternative testing. according to these results org.Hs.egSYMBOL provides better coverage of ENTREZ-SYMBOL combinations
# ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# tr_dt2 <- getBM(attributes=c("entrezgene", "hgnc_symbol"),mart=ensembl)
# length(sort(unique(tr_dt2[,2]))); length(sort(unique(cc[,2])))
# length(unique(tr_dt2[complete.cases(tr_dt2),2])); length(unique(cc[complete.cases(cc),2]))

# ensembl = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
# x <- getBM(attributes=c("flybase_transcript_id", "entrezgene","name_1006", "namespace_1003"), values=trnscrNames, mart=ensembl)
# tr_dt <- getBM(attributes=c("flybase_transcript_id", "entrezgene","flybasename_gene", "name_1006", "namespace_1003"),mart=ensembl)
# symb <- as.data.frame(org.Dm.egSYMBOL)
# length(sort(unique(tr_dt[,3]))); length(sort(unique(symb[,2])))


# tr_orgdm <- mappedkeys(org.Dm.egENSEMBLTRANS2EG)
# tr_mrtdm <- tr_dt[,1]
# length(sort(unique(tr_orgdm))); length(sort(unique(tr_mrtdm)))

