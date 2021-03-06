plotHeatmap <- function(dexp,contrasts,gsets,gsetgroup,type, pcutoff) {    
    
    v_gs_sg_names <- v_contr <- vector() 
    rm(lgfc, v_g_pvalues, m_tr_fc_simpl, v_gs_pvalues, m_gs_pvalues)
    gsets[is.na(gsets[,1])==TRUE,1] <- "missing"; gsets[is.na(gsets[,2])==TRUE,2] <- "missing"
    gsets[is.na(gsets[,3])==TRUE,3] <- "missing"; gsets[is.na(gsets[,4])==TRUE,4] <- "missing"
    
    #########################
    # calculate significant gene gsets
    for (contrast in contrasts) {
        a <- unlist(strsplit(contrast, split=" - "))
        v_contr <- append(v_contr, a[2])
        
        if (type=="pval"){
            v_g_pvalues <- dexp$pval[[contrast]]
        } else {
            v_g_pvalues <- dexp$padj[[contrast]]
        }; v_g_pvalues[v_g_pvalues==0] <- 1e-10; v_g_pvalues[is.na(v_g_pvalues)] <- 1
        
        
        # calculate initial values for non-regulated significant gene gsets
        if (file.exists(paste("./Data/GSA_hyper_",gsetgroup,contrast,"_",pcutoff,".rds",sep=""))) {
            gsah <- readRDS(paste("./Data/GSA_hyper_",gsetgroup,contrast,"_",pcutoff,".rds",sep=""))
        } else {
            gsah <- enrichSets(dexp = dexp, 
                               pvalues = v_g_pvalues, 
                               gsets = gsets, 
                               pcutoff = pcutoff)
            saveRDS(gsah,paste("./Data/GSA_hyper_",gsetgroup,contrast,"_",pcutoff,".rds",sep="")) 
        } 
        
        # either p-value or adjusted p-value
        if (type=="pval"){
            v_gs_pvalues <- gsah$pval
        } else {
            v_gs_pvalues <- gsah$padj
        }
        
        v_gs_sg_names <- c(v_gs_sg_names,gsah$gset[v_gs_pvalues <= pcutoff])

        if (exists("m_gs_pvalues")==FALSE) {
            m_gs_pvalues <- matrix(nrow=length(v_gs_pvalues),ncol=length(contrasts))
            colnames(m_gs_pvalues) <- contrasts
            rownames(m_gs_pvalues) <- gsah$gset
            m_gs_pvalues[,contrast] <- v_gs_pvalues
        } else {
            m_gs_pvalues[,contrast] <- v_gs_pvalues
        }
        
        
        # calculate initial values for up- & down-regulated significant gene gsets
        lgfc <- dexp$lgFC[[contrast]]; lgfc[is.na(lgfc)] <- 0
        
        if (exists("m_tr_fc_simpl")==FALSE) {
            m_tr_fc_simpl <- matrix(nrow=length(dexp$name),ncol=length(contrasts))
            colnames(m_tr_fc_simpl) <- contrasts
            rownames(m_tr_fc_simpl) <- dexp$name
            m_tr_fc_simpl[which(v_g_pvalues <= pcutoff & lgfc < 0),contrast] <- -1
            m_tr_fc_simpl[which(v_g_pvalues <= pcutoff & lgfc > 0),contrast] <- 1
            m_tr_fc_simpl[which(v_g_pvalues <= pcutoff & lgfc == 0),contrast] <- 0
            m_tr_fc_simpl[which(v_g_pvalues > pcutoff),contrast] <- 0
        } else {
            m_tr_fc_simpl[which(v_g_pvalues <= pcutoff & lgfc < 0),contrast] <- -1
            m_tr_fc_simpl[which(v_g_pvalues <= pcutoff & lgfc > 0),contrast] <- 1
            m_tr_fc_simpl[which(v_g_pvalues <= pcutoff & lgfc == 0),contrast] <- 0
            m_tr_fc_simpl[which(v_g_pvalues > pcutoff),contrast] <- 0
        }
    }
    # v_gs_sg_names -- vector__genesets__significant__names
    v_gs_sg_names <- sort(unique(v_gs_sg_names))
    
    # calculate initial values for non-regulated significant gene sets
    # df_gs_sga_simpl -- data.frame__genesets__signigicant_all__simplified
    df_gs_sga_simpl <- as.data.frame(m_gs_pvalues[rownames(m_gs_pvalues) %in% v_gs_sg_names,])
    # simplification
    df_gs_sga_simpl[df_gs_sga_simpl > pcutoff] <- 2
    df_gs_sga_simpl[df_gs_sga_simpl <= pcutoff] <- 1
    df_gs_sga_simpl[df_gs_sga_simpl == 2] <- 0

    # calculate initial values for up- & down-regulated significant gene sets
    df_tr_fc_simpl <- as.data.frame(m_tr_fc_simpl); 
    df_tr_fc_simpl$FBtranscriptID <- rownames(df_tr_fc_simpl)
    
    df_sets_gs_simpl <- merge(gsets,df_tr_fc_simpl,by="FBtranscriptID")
    df_sets_gs_sg_simpl <- df_sets_gs_simpl[df_sets_gs_simpl$GSET %in% v_gs_sg_names,]
    
    
    # df_gs_sgp_simpl -- dataframe__geneset__signigicant_positive__simplified
    df_gs_sgp_simpl <- setkey(setDT(df_sets_gs_sg_simpl[,c(3,5:length(df_sets_gs_sg_simpl))]), GSET)[, lapply(.SD, max), GSET]
    df_gs_sgp_simpl <- as.data.frame(df_gs_sgp_simpl)
    #
    df_gs_sgp_simpl[df_gs_sgp_simpl==-1] <- 0
    df_gs_sgp_simpl[df_gs_sgp_simpl==1] <- 1
    df_gs_sgp_simpl[df_gs_sgp_simpl==0] <- 0
    rownames(df_gs_sgp_simpl) <- df_gs_sgp_simpl$GSET
    #
    df_gs_sgp_simpl <- df_gs_sga_simpl + df_gs_sgp_simpl[,2:length(df_gs_sgp_simpl)]
    df_gs_sgp_simpl[df_gs_sgp_simpl == 1] <- 0
    df_gs_sgp_simpl[df_gs_sgp_simpl == 2] <- 1
    #
    # df_gs_sgn_simpl -- dataframe__geneset__signigicant_negative__simplified
    df_gs_sgn_simpl <- setkey(setDT(df_sets_gs_sg_simpl[,c(3,5:length(df_sets_gs_sg_simpl))]), GSET)[, lapply(.SD, min), GSET]
    df_gs_sgn_simpl <- as.data.frame(df_gs_sgn_simpl)
    #
    df_gs_sgn_simpl[df_gs_sgn_simpl==1] <- 0
    df_gs_sgn_simpl[df_gs_sgn_simpl==-1] <- 1
    df_gs_sgn_simpl[df_gs_sgn_simpl==0] <- 0
    rownames(df_gs_sgn_simpl) <- df_gs_sgn_simpl$GSET
    #
    df_gs_sgn_simpl <- df_gs_sga_simpl + df_gs_sgn_simpl[,2:length(df_gs_sgn_simpl)]
    df_gs_sgn_simpl[df_gs_sgn_simpl == 1] <- 0
    df_gs_sgn_simpl[df_gs_sgn_simpl == 2] <- 1
    
    
        
    #########################
    # calculate DE-related genes
    # prepare data
    df_sets_gs_simpl$FBtranscriptID <- NULL
    df_sets_gs_simpl_a <- df_sets_gs_simpl_p <- df_sets_gs_simpl_n <- df_sets_gs_simpl
    #
    for (i in 4:(length(v_contr)+3)) {
        df_sets_gs_simpl_a[df_sets_gs_simpl_a[,i] != 0,i] <- 1
        #
        df_sets_gs_simpl_p[df_sets_gs_simpl_p[,i] != 1,i] <- 0
        #
        df_sets_gs_simpl_n[df_sets_gs_simpl_n[,i] != -1,i] <- 0
        df_sets_gs_simpl_n[df_sets_gs_simpl_n[,i] == -1,i] <- 1
    }

    
    countDEgenes <- function(df_gene_simpl){
        
        df_gene_simpl <- df_gene_simpl[!duplicated(df_gene_simpl),]
        df_gene_simpl$ENTREZID <- df_gene_simpl$GNAME <- NULL
        df_gene_simpl_count  <- setkey(setDT(df_gene_simpl), GSET)[, lapply(.SD, sum), GSET]; df_gene_simpl_count <- as.data.frame(df_gene_simpl_count)
        
        rownames(df_gene_simpl_count) <- df_gene_simpl_count$GSET
        df_gene_simpl_count$GSET <- NULL
        
        df_gene_simpl_count <- df_gene_simpl_count[which(rowSums(df_gene_simpl_count)>0),]
        m_gene_simpl_count <- as.matrix(df_gene_simpl_count)
        
        return(m_gene_simpl_count)
        
    }
    
    
    listGenes <- function (df_sets_gs_simpl, m_gs_count, gsets, v_contr){
        
        v_gs_names <- names(as.matrix(m_gs_count)[,1])
        df_entrList <- data.frame(matrix(NA, nrow = dim(m_gs_count)[1], ncol = (dim(m_gs_count)[2]*2+1)))
        
        for (i in 2:(length(v_contr)+1)){
            v_col <- c(1,2,3,(i+2))
            df_tmpr <- df_sets_gs_simpl[,v_col]
            for (j in 1:length(v_gs_names)) {
                tmp_entr <- paste(unique(df_tmpr[which(df_tmpr[,4]>0 & df_tmpr[,2]==v_gs_names[j]),1]), collapse=", ")
                tmp_name <- paste(unique(df_tmpr[which(df_tmpr[,4]>0 & df_tmpr[,2]==v_gs_names[j]),3]), collapse=", ")
                if ( m_gs_count[j,(i-1)] != 0) {df_entrList[j,i] <- tmp_entr; df_entrList[j,(i+length(v_contr))] <- tmp_name}
            }
        }
        
        colnames(df_entrList) <- c("GSET",paste(v_contr,"entr",sep="_"),paste(v_contr,"name",sep="_")) 
        rownames(df_entrList) <- v_gs_names
        df_entrList$GSET <- v_gs_names
        return(df_entrList)
    }
    
    
    m_gset_gn_counta <- countDEgenes(df_sets_gs_simpl_a)
    m_gset_gn_countp <- countDEgenes(df_sets_gs_simpl_p)
    m_gset_gn_countn <- countDEgenes(df_sets_gs_simpl_n)
    
    # prepare matrices for plotting
    colnames(df_gs_sga_simpl) <- colnames(df_gs_sgp_simpl) <- colnames(df_gs_sgn_simpl) <- v_contr
    m_gs_sga_simpl <- as.matrix(df_gs_sga_simpl); m_gs_sgp_simpl <- as.matrix(df_gs_sgp_simpl); m_gs_sgn_simpl <- as.matrix(df_gs_sgn_simpl); 
    colnames(m_gset_gn_counta) <- colnames(m_gset_gn_countp) <- colnames(m_gset_gn_countn) <- v_contr
    
    
    df_entrList_sa <- listGenes(df_sets_gs_simpl_a, m_gs_sga_simpl, gsets, v_contr)
    df_entrList_sp <- listGenes(df_sets_gs_simpl_p, m_gs_sgp_simpl, gsets, v_contr)
    df_entrList_sn <- listGenes(df_sets_gs_simpl_n, m_gs_sgn_simpl, gsets, v_contr)
    
    df_entrList_dea <- listGenes(df_sets_gs_simpl_a, m_gset_gn_counta, gsets, v_contr)
    df_entrList_dep <- listGenes(df_sets_gs_simpl_p, m_gset_gn_countp, gsets, v_contr)
    df_entrList_den <- listGenes(df_sets_gs_simpl_n, m_gset_gn_countn, gsets, v_contr)

    
    # draw heatmaps
    # sg means those genesets having small p.adj of genesets
    if (length(m_gs_sga_simpl)>0) {
        png(paste("./Figures/heatmap_",gsetgroup,"_sg_all_",type,"_",pcutoff,".png",sep=""),width=1200,height=(50*dim(m_gs_sga_simpl)[1]+1800),res=300,pointsize=8)
        heatmap.2(m_gs_sga_simpl[!(rownames(m_gs_sga_simpl) %in% "missing"),], col=colorRampPalette(brewer.pal(9,"GnBu"))(100),lmat=rbind(c(4,3),c(1,2)),lhei=c(1,10),lwid=c(6,4), margins=c(17,9), Key=TRUE, keysize=0.2, dendrogram="none", density.info="none", trace="none",Rowv = FALSE, Colv=FALSE); mtext(paste("All significant gene sets, ",type, "=",pcutoff,sep=""), side = 1, line=3, cex=1.5) 
        dev.off()
        write.table(df_entrList_sa,file=paste("./Figures/heatmap_",gsetgroup,"_sg_all_",type,"_",pcutoff,".csv",sep=""),sep=";",na="",row.names=FALSE,col.names=TRUE)
    }
    if (length(m_gs_sgp_simpl)>0) {
        png(paste("./Figures/heatmap_",gsetgroup,"_sg_pos_",type,"_",pcutoff,".png",sep=""),width=1200,height=(50*dim(m_gs_sgp_simpl)[1]+1800),res=300,pointsize=8)
        heatmap.2(m_gs_sgp_simpl[!(rownames(m_gs_sgp_simpl) %in% "missing"),], col=colorRampPalette(brewer.pal(9,"GnBu"))(100),lmat=rbind(c(4,3),c(1,2)),lhei=c(1,10),lwid=c(6,4), margins=c(17,9), Key=TRUE, keysize=0.2, dendrogram="none", density.info="none", trace="none",Rowv = FALSE, Colv=FALSE); mtext(paste("Up-regulated significant gene sets, ",type, "=",pcutoff,sep=""), side = 1, line=3, cex=1.5) 
        dev.off()
        write.table(df_entrList_sp,file=paste("./Figures/heatmap_",gsetgroup,"_sg_pos_",type,"_",pcutoff,".csv",sep=""),sep=";",na="",row.names=FALSE,col.names=TRUE)
    }
    if (length(m_gs_sgn_simpl)>0) {
        png(paste("./Figures/heatmap_",gsetgroup,"_sg_neg_",type,"_",pcutoff,".png",sep=""),width=1200,height=(50*dim(m_gs_sgn_simpl)[1]+1800),res=300,pointsize=8)
        heatmap.2(m_gs_sgn_simpl[!(rownames(m_gs_sgn_simpl) %in% "missing"),], col=colorRampPalette(brewer.pal(9,"GnBu"))(100),lmat=rbind(c(4,3),c(1,2)),lhei=c(1,10),lwid=c(6,4), margins=c(18,9), Key=TRUE, keysize=1, dendrogram="none", density.info="none", trace="none",Rowv = FALSE, Colv=FALSE); mtext(paste("Down-regulated significant gene sets, ",type, "=",pcutoff,sep=""), side = 1, line=3, cex=1.5) 
        dev.off()
        write.table(df_entrList_sn,file=paste("./Figures/heatmap_",gsetgroup,"_sg_neg_",type,"_",pcutoff,".csv",sep=""),sep=";",na="",row.names=FALSE,col.names=TRUE)
    }
    
    # de means those genesets having genes with small p value
    if (length(m_gset_gn_counta)>0) {
        png(paste("./Figures/heatmap_",gsetgroup,"_de_all_",type,"_",pcutoff,".png",sep=""),width=1200,height=(50*dim(m_gset_gn_counta)[1]+1800),res=300,pointsize=8)
        heatmap.2(m_gset_gn_counta[!(rownames(m_gset_gn_counta) %in% "missing"),],col=colorRampPalette(brewer.pal(9,"GnBu"))(100),lmat=rbind(c(4,3),c(1,2)),lhei=c(1,10),lwid=c(6,4), margins=c(17,9), Key=TRUE, keysize=1, dendrogram="none", density.info="none", trace="none",Rowv = FALSE, Colv=FALSE); mtext(paste("All DE gene sets, ",type, "=",pcutoff,sep=""), side = 1, line=3, cex=1.5)
        dev.off()
        write.table(df_entrList_dea,file=paste("./Figures/heatmap_",gsetgroup,"_de_all_",type,"_",pcutoff,".csv",sep=""),sep=";",na="",row.names=FALSE,col.names=TRUE)
    }
    if (length(m_gset_gn_countp)>0) {
        png(paste("./Figures/heatmap_",gsetgroup,"_de_pos_",type,"_",pcutoff,".png",sep=""),width=1200,height=(50*dim(m_gset_gn_countp)[1]+1800),res=300,pointsize=8)
        heatmap.2(m_gset_gn_countp[!(rownames(m_gset_gn_countp) %in% "missing"),],col=colorRampPalette(brewer.pal(9,"GnBu"))(100), ,lmat=rbind(c(4,3),c(1,2)),lhei=c(1,10),lwid=c(6,4), margins=c(17,9), Key=TRUE, keysize=0.2, dendrogram="none", density.info="none", trace="none",Rowv = FALSE, Colv=FALSE); mtext(paste("Up-regulated DE gene sets, ",type, "=",pcutoff,sep=""), side = 1, line=3, cex=1.5)
        dev.off()
        write.table(df_entrList_dep,file=paste("./Figures/heatmap_",gsetgroup,"_de_pos_",type,"_",pcutoff,".csv",sep=""),sep=";",na="",row.names=FALSE,col.names=TRUE)
    }
    if (length(m_gset_gn_countn)>0) {
        png(paste("./Figures/heatmap_",gsetgroup,"_de_neg_",type,"_",pcutoff,".png",sep=""),width=1200,height=(50*dim(m_gset_gn_countn)[1]+1800),res=300,pointsize=8)
        heatmap.2(m_gset_gn_countn[!(rownames(m_gset_gn_countn) %in% "missing"),],col=colorRampPalette(brewer.pal(9,"GnBu"))(100),lmat=rbind(c(4,3),c(1,2)),lhei=c(1,10),lwid=c(6,4), margins=c(17,9), Key=TRUE, keysize=0.2, dendrogram="none", density.info="none", trace="none",Rowv = FALSE, Colv=FALSE); mtext(paste("Down-regulated DE gene sets, ",type, "=",pcutoff,sep=""), side = 1, line=3, cex=1.5)
        dev.off()
        write.table(df_entrList_den,file=paste("./Figures/heatmap_",gsetgroup,"_de_neg_",type,"_",pcutoff,".csv",sep=""),sep=";",na="",row.names=FALSE,col.names=TRUE)
    }
    
#     # sgde means those genesets having small p.adj of genesets & showing No of genes with small p-value
#     if (length(m_gset_gn_counta)>0) {
#         png(paste("./Figures/heatmap_",gsetgroup,"_de_all_",type,"_",pcutoff,".png",sep=""),width=1200,height=(50*dim(m_gset_gn_counta)[1]+1800),res=300,pointsize=8)
#         heatmap.2(m_gset_gn_counta[which(!(rownames(m_gset_gn_counta) %in% "missing") & 
#                                           (rownames(m_gset_gn_counta) %in% rownames(m_gs_sga_simpl)))
#                                    ,],col=colorRampPalette(brewer.pal(9,"GnBu"))(100),lmat=rbind(c(4,3),c(1,2)),lhei=c(1,10),lwid=c(6,4), margins=c(17,9), Key=TRUE, keysize=1, dendrogram="none", density.info="none", trace="none",Rowv = FALSE, Colv=FALSE); mtext(paste("All DE gene sets, ",type, "=",pcutoff,sep=""), side = 1, line=3, cex=1.5)
#         dev.off()
# #         write.table(df_entrList_dea,file=paste("./Figures/heatmap_",gsetgroup,"_de_all_",type,"_",pcutoff,".csv",sep=""),sep=";",na="",row.names=FALSE,col.names=TRUE)
#     }
#     if (length(m_gset_gn_countp)>0) {
#         png(paste("./Figures/heatmap_",gsetgroup,"_de_pos_",type,"_",pcutoff,".png",sep=""),width=1200,height=(50*dim(m_gset_gn_countp)[1]+1800),res=300,pointsize=8)
#         heatmap.2(m_gset_gn_countp[which(!(rownames(m_gset_gn_countp) %in% "missing") & 
#                                              (rownames(m_gset_gn_countp) %in% rownames(m_gs_sgp_simpl)) )
#                                          ,],col=colorRampPalette(brewer.pal(9,"GnBu"))(100), ,lmat=rbind(c(4,3),c(1,2)),lhei=c(1,10),lwid=c(6,4), margins=c(17,9), Key=TRUE, keysize=0.2, dendrogram="none", density.info="none", trace="none",Rowv = FALSE, Colv=FALSE); mtext(paste("Up-regulated DE gene sets, ",type, "=",pcutoff,sep=""), side = 1, line=3, cex=1.5)
#         dev.off()
# #         write.table(df_entrList_dep,file=paste("./Figures/heatmap_",gsetgroup,"_de_pos_",type,"_",pcutoff,".csv",sep=""),sep=";",na="",row.names=FALSE,col.names=TRUE)
#     }
#     if (length(m_gset_gn_countn)>0) {
#         png(paste("./Figures/heatmap_",gsetgroup,"_de_neg_",type,"_",pcutoff,".png",sep=""),width=1200,height=(50*dim(m_gset_gn_countn)[1]+1800),res=300,pointsize=8)
#         heatmap.2(m_gset_gn_countn[which(!(rownames(m_gset_gn_countn) %in% "missing") & 
#                                              (rownames(m_gset_gn_countn) %in% rownames(m_gs_sgn_simpl)) 
#                                    ,],col=colorRampPalette(brewer.pal(9,"GnBu"))(100),lmat=rbind(c(4,3),c(1,2)),lhei=c(1,10),lwid=c(6,4), margins=c(17,9), Key=TRUE, keysize=0.2, dendrogram="none", density.info="none", trace="none",Rowv = FALSE, Colv=FALSE); mtext(paste("Down-regulated DE gene sets, ",type, "=",pcutoff,sep=""), side = 1, line=3, cex=1.5)
#         dev.off()
# #         write.table(df_entrList_den,file=paste("./Figures/heatmap_",gsetgroup,"_de_neg_",type,"_",pcutoff,".csv",sep=""),sep=";",na="",row.names=FALSE,col.names=TRUE)
#     }
#     

#     
}




