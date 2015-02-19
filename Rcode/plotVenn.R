plotVenn <- function(dexp, contrasts, pvalue=0.05, adj.pvalue=0.1) {
    v_contr <- vector() 
    rm(lgfc, v_g_pvalues, v_g_avalues, m_tr_fc_simpl_pv, m_tr_fc_simpl_av, sets)
    
    
    trnscrNames <- mappedkeys(org.Dm.egENSEMBLTRANS2EG); cols <- c("ENTREZID")
    anno <- select(org.Dm.eg.db, keys = trnscrNames, columns = cols, keytype = "ENSEMBLTRANS")
    
    ensembl = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
    x <- getBM(attributes=c("flybase_transcript_id", "entrezgene"), values=dexp$name, mart=ensembl)
    
    colnames(anno) <- colnames(x) <- c("FBtranscriptID", "ENTREZID")
    sets <- rbind(anno, x)
    sets <- sets[!duplicated(sets),]
    colnames(sets) <- c("FBtranscriptID", "ENTREZID")
    
#     vt <- df_tr_fc_simpl_av[which(df_tr_fc_simpl_av[,1]!=0 | df_tr_fc_simpl_av[,2]!=0 | df_tr_fc_simpl_av[,3]!=0),4]
#     vt[!(vt %in% sets[,1])]
#     
#     vt <- df_tr_fc_simpl_pv[which(df_tr_fc_simpl_pv[,1]!=0 | df_tr_fc_simpl_pv[,2]!=0 | df_tr_fc_simpl_pv[,3]!=0),4]
#     vt[!(vt %in% sets[,1])]
    

    sets[is.na(sets$FBtranscriptID)==TRUE,1] <- "unknown"; sets[is.na(sets$ENTREZID)==TRUE,2] <- "unknown"
    for (contrast in contrasts) {
        a <- unlist(strsplit(contrast, split=" - "))
        v_contr <- append(v_contr, a[2])
        
        
        v_g_pvalues <- dexp$pval[[contrast]]; v_g_pvalues[v_g_pvalues==0] <- 1e-10; v_g_pvalues[is.na(v_g_pvalues)] <- 1
        v_g_avalues <- dexp$padj[[contrast]]; v_g_avalues[v_g_avalues==0] <- 1e-10; v_g_avalues[is.na(v_g_avalues)] <- 1
        
        
        lgfc <- dexp$lgFC[[contrast]]; lgfc[is.na(lgfc)] <- 0
        
        if (exists("m_tr_fc_simpl_pv")==FALSE) {
            m_tr_fc_simpl_pv <- m_tr_fc_simpl_av <- matrix(nrow=length(dexp$name),ncol=length(contrasts))
            colnames(m_tr_fc_simpl_pv) <- colnames(m_tr_fc_simpl_av) <- contrasts
            rownames(m_tr_fc_simpl_pv) <- rownames(m_tr_fc_simpl_av) <- dexp$name
            m_tr_fc_simpl_pv[which(v_g_pvalues <= pvalue & lgfc < 0),contrast] <- -1
            m_tr_fc_simpl_pv[which(v_g_pvalues <= pvalue & lgfc > 0),contrast] <- 1
            m_tr_fc_simpl_pv[which(v_g_pvalues <= pvalue & lgfc == 0),contrast] <- 0
            m_tr_fc_simpl_pv[which(v_g_pvalues > pvalue),contrast] <- 0
            
            m_tr_fc_simpl_av[which(v_g_avalues <= adj.pvalue & lgfc < 0),contrast] <- -1
            m_tr_fc_simpl_av[which(v_g_avalues <= adj.pvalue & lgfc > 0),contrast] <- 1
            m_tr_fc_simpl_av[which(v_g_avalues <= adj.pvalue & lgfc == 0),contrast] <- 0
            m_tr_fc_simpl_av[which(v_g_avalues > adj.pvalue),contrast] <- 0
        } else {
            m_tr_fc_simpl_pv[which(v_g_pvalues <= pvalue & lgfc < 0),contrast] <- -1
            m_tr_fc_simpl_pv[which(v_g_pvalues <= pvalue & lgfc > 0),contrast] <- 1
            m_tr_fc_simpl_pv[which(v_g_pvalues <= pvalue & lgfc == 0),contrast] <- 0
            m_tr_fc_simpl_pv[which(v_g_pvalues > pvalue),contrast] <- 0
            
            m_tr_fc_simpl_av[which(v_g_avalues <= adj.pvalue & lgfc < 0),contrast] <- -1
            m_tr_fc_simpl_av[which(v_g_avalues <= adj.pvalue & lgfc > 0),contrast] <- 1
            m_tr_fc_simpl_av[which(v_g_avalues <= adj.pvalue & lgfc == 0),contrast] <- 0
            m_tr_fc_simpl_av[which(v_g_avalues > adj.pvalue),contrast] <- 0
        }
    }
    # calculate initial values for up- & down-regulated significant gene sets
    df_tr_fc_simpl_pv <- as.data.frame(m_tr_fc_simpl_pv) 
    df_tr_fc_simpl_pv$FBtranscriptID <- rownames(df_tr_fc_simpl_pv)
    #
    df_gn_fc_simpl_pv <- merge(sets,df_tr_fc_simpl_pv,by="FBtranscriptID")
    #
    #
    df_tr_fc_simpl_av <- as.data.frame(m_tr_fc_simpl_av); 
    df_tr_fc_simpl_av$FBtranscriptID <- rownames(df_tr_fc_simpl_av)
    #
    df_gn_fc_simpl_av <- merge(sets,df_tr_fc_simpl_av,by="FBtranscriptID")
    
    # prepare data
    df_gn_fc_simpl_pv$FBtranscriptID <- NULL; df_gn_fc_simpl_av$FBtranscriptID <- NULL
    df_gn_fc_simpl_pv_a <- df_gn_fc_simpl_pv_p <- df_gn_fc_simpl_pv_n <- df_gn_fc_simpl_pv
    df_gn_fc_simpl_av_a <- df_gn_fc_simpl_av_p <- df_gn_fc_simpl_av_n <- df_gn_fc_simpl_av
    #
    df_tr_fc_simpl_pv_a <- df_tr_fc_simpl_pv_p <- df_tr_fc_simpl_pv_n <- df_tr_fc_simpl_pv
    df_tr_fc_simpl_av_a <- df_tr_fc_simpl_av_p <- df_tr_fc_simpl_av_n <- df_tr_fc_simpl_av
    #
    for (i in 2:(length(v_contr)+1)) {
        df_gn_fc_simpl_pv_a[df_gn_fc_simpl_pv_a[,i] != 0,i] <- 1
        #
        df_gn_fc_simpl_pv_p[df_gn_fc_simpl_pv_p[,i] != 1,i] <- 0
        #
        df_gn_fc_simpl_pv_n[df_gn_fc_simpl_pv_n[,i] != -1,i] <- 0
        df_gn_fc_simpl_pv_n[df_gn_fc_simpl_pv_n[,i] == -1,i] <- 1
        #
        #
        df_gn_fc_simpl_av_a[df_gn_fc_simpl_av_a[,i] != 0,i] <- 1
        #
        df_gn_fc_simpl_av_p[df_gn_fc_simpl_av_p[,i] != 1,i] <- 0
        #
        df_gn_fc_simpl_av_n[df_gn_fc_simpl_av_n[,i] != -1,i] <- 0
        df_gn_fc_simpl_av_n[df_gn_fc_simpl_av_n[,i] == -1,i] <- 1
    }
    #
    for (i in 1:length(v_contr)) {
        df_tr_fc_simpl_pv_a[df_tr_fc_simpl_pv_a[,i] != 0,i] <- 1
        #
        df_tr_fc_simpl_pv_p[df_tr_fc_simpl_pv_p[,i] != 1,i] <- 0
        #
        df_tr_fc_simpl_pv_n[df_tr_fc_simpl_pv_n[,i] != -1,i] <- 0
        df_tr_fc_simpl_pv_n[df_tr_fc_simpl_pv_n[,i] == -1,i] <- 1
        #
        #
        df_tr_fc_simpl_av_a[df_tr_fc_simpl_av_a[,i] != 0,i] <- 1
        #
        df_tr_fc_simpl_av_p[df_tr_fc_simpl_av_p[,i] != 1,i] <- 0
        #
        df_tr_fc_simpl_av_n[df_tr_fc_simpl_av_n[,i] != -1,i] <- 0
        df_tr_fc_simpl_av_n[df_tr_fc_simpl_av_n[,i] == -1,i] <- 1
    }
    
    
    countDEtrnscr <- function(df_trnscr_simpl){
        df_trnscr_simpl <- df_trnscr_simpl[!duplicated(df_trnscr_simpl),]
        df_trnscr_simpl$FBtranscriptID <- NULL
#         df_trnscr_simplt <- df_trnscr_simpl[which(rowSums(df_trnscr_simpl)>0),]
        m_trnscr_simpl <- as.matrix(df_trnscr_simpl)
        return(m_trnscr_simpl)
    }
    
    
    countDEgenes <- function(df_gene_simpl){
        
        df_gene_simpl <- df_gene_simpl[!duplicated(df_gene_simpl),]
        df_gene_simpl_count  <- setkey(setDT(df_gene_simpl), ENTREZID)[, lapply(.SD, sum), ENTREZID];   df_gene_simpl_count <- as.data.frame(df_gene_simpl_count)
        
        rownames(df_gene_simpl_count) <- df_gene_simpl_count$ENTREZID
        df_gene_simpl_count$ENTREZID <- NULL
        
#         df_gene_simpl_count <- df_gene_simpl_count[which(rowSums(df_gene_simpl_count)>0),]
        
        m_gene_simpl_count <- as.matrix(df_gene_simpl_count)
        return(m_gene_simpl_count)
        
    }
    
    
    m_tr_counta_pv <- countDEtrnscr(df_tr_fc_simpl_pv_a)
    m_tr_countp_pv <- countDEtrnscr(df_tr_fc_simpl_pv_p)
    m_tr_countn_pv <- countDEtrnscr(df_tr_fc_simpl_pv_n)
    #
    m_tr_counta_av <- countDEtrnscr(df_tr_fc_simpl_av_a)
    m_tr_countp_av <- countDEtrnscr(df_tr_fc_simpl_av_p)
    m_tr_countn_av <- countDEtrnscr(df_tr_fc_simpl_av_n)
    #
    #
    m_gn_counta_pv <- countDEgenes(df_gn_fc_simpl_pv_a)
    m_gn_countp_pv <- countDEgenes(df_gn_fc_simpl_pv_p)
    m_gn_countn_pv <- countDEgenes(df_gn_fc_simpl_pv_n)
    #
    m_gn_counta_av <- countDEgenes(df_gn_fc_simpl_av_a)
    m_gn_countp_av <- countDEgenes(df_gn_fc_simpl_av_p)
    m_gn_countn_av <- countDEgenes(df_gn_fc_simpl_av_n)
    
    
    nm <- LETTERS[seq(from = 1, to = length(v_contr))]
    tx <- paste(nm, ': ',v_contr, '; ', sep="", collapse="")
    ps <- gregexpr("; ", tx)[[1]][2]
    substr(tx,start=ps, stop=(ps+1)) <- "\n"
    png(paste("./Figures/venn_all_",v_contr[1],".png",sep=""),width=2100,height=2100,res=300,pointsize=8)
    par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
    vennDiagram(m_tr_counta_pv,names=nm,cex=1,main=paste("\n\nTranscripts, p-value <= ",pvalue,sep=""),cex.main=1)
    vennDiagram(m_tr_counta_av,names=nm,cex=1,main=paste("\n\nTranscripts, FDR adjusted p-value <= ",adj.pvalue,sep=""),cex.main=1)
    vennDiagram(m_gn_counta_pv,names=nm,cex=1,main=paste("\n\nGenes, p-value <= ",pvalue,sep=""),cex.main=1)
    vennDiagram(m_gn_counta_av,names=nm,cex=1,main=paste("\n\nGenes, FDR adjusted p-value <= ",adj.pvalue,sep=""),cex.main=1)
    par(mfrow=c(1,1), mar=c(0,3,1,1))
    mtext("Differentially expressed transcripts and genes.\n\n",side=2,cex=1.3)
    mtext(tx,side=1, cex=0.95)
    dev.off()
    
    
    nm <- LETTERS[seq(from = 1, to = length(v_contr))]
    tx <- paste(nm, ': ',v_contr, '; ', sep="", collapse="")
    ps <- gregexpr("; ", tx)[[1]][2]
    substr(tx,start=ps, stop=(ps+1)) <- "\n"
    png(paste("./Figures/venn_up_",v_contr[1],".png",sep=""),width=2100,height=2100,res=300,pointsize=8)
    par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
    vennDiagram(m_tr_countp_pv,names=nm,cex=1,main=paste("\n\nTranscripts, p-value <= ",pvalue,sep=""),cex.main=1)
    vennDiagram(m_tr_countp_av,names=nm,cex=1,main=paste("\n\nTranscripts, FDR adjusted p-value <= ",adj.pvalue,sep=""),cex.main=1)
    vennDiagram(m_gn_countp_pv,names=nm,cex=1,main=paste("\n\nGenes, p-value <= ",pvalue,sep=""),cex.main=1)
    vennDiagram(m_gn_countp_av,names=nm,cex=1,main=paste("\n\nGenes, FDR adjusted p-value <= ",adj.pvalue,sep=""),cex.main=1)
    par(mfrow=c(1,1), mar=c(0,3,1,1))
    mtext("Up-regulated transcripts and genes.\n\n",side=2,cex=1.3)
    mtext(tx,side=1, cex=0.95)
    dev.off()


    nm <- LETTERS[seq(from = 1, to = length(v_contr))]
    tx <- paste(nm, ': ',v_contr, '; ', sep="", collapse="")
    ps <- gregexpr("; ", tx)[[1]][2]
    substr(tx,start=ps, stop=(ps+1)) <- "\n"
    png(paste("./Figures/venn_down_",v_contr[1],".png",sep=""),width=2100,height=2100,res=300,pointsize=8)
    par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
    vennDiagram(m_tr_countn_pv,names=nm,cex=1,main=paste("\n\nTranscripts, p-value <= ",pvalue,sep=""),cex.main=1)
    vennDiagram(m_tr_countn_av,names=nm,cex=1,main=paste("\n\nTranscripts, FDR adjusted p-value <= ",adj.pvalue,sep=""),cex.main=1)
    vennDiagram(m_gn_countn_pv,names=nm,cex=1,main=paste("\n\nGenes, p-value <= ",pvalue,sep=""),cex.main=1)
    vennDiagram(m_gn_countn_av,names=nm,cex=1,main=paste("\n\nGenes, FDR adjusted p-value <= ",adj.pvalue,sep=""),cex.main=1)
    par(mfrow=c(1,1), mar=c(0,3,1,1))
    mtext("Down-regulated transcripts and genes.\n\n",side=2,cex=1.3)
    mtext(tx,side=1, cex=0.95)
    dev.off()
}
