plotNetwork <- function(dexp, contrast, sets, enrset, pcutoff, type, gsetgroup) {
    
    if (type=="pval"){
        pvalues <- dexp$pval[[contrast]]
    } else {
        pvalues <- dexp$padj[[contrast]]
    }; pvalues[pvalues==0] <- 1e-10
    
    dt <- data.frame("FBtranscriptID"=dexp$name, 
                     "padj"=pvalues, 
                     "lgfc"=dexp$lgFC[[contrast]]) 
    dt$FBtranscriptID <- as.character(dt$FBtranscriptID)
    dt[is.na(dt$padj),2] <- 1; dt[is.na(dt$lgfc),3] <- 0 
    sets[is.na(sets$ENTREZID)==TRUE,2] <- "unknown"; sets[is.na(sets$GSET)==TRUE,3] <- "unknown"
    
    dtd <- merge(dt, sets,by="FBtranscriptID")
    
    if (type=="pval"){
        lst_n <- dtd[which(dtd$GSET %in% enrset[enrset$pval <= pcutoff,"gset"] & dtd$padj <= pcutoff & dtd$lgfc < 0),c(5,4)]
        lst_p <- dtd[which(dtd$GSET %in% enrset[enrset$pval <= pcutoff,"gset"] & dtd$padj <= pcutoff & dtd$lgfc > 0),c(5,4)]
    } else {
        lst_n <- dtd[which(dtd$GSET %in% enrset[enrset$padj <= pcutoff,"gset"] & dtd$padj <= pcutoff & dtd$lgfc < 0),c(5,4)]
        lst_p <- dtd[which(dtd$GSET %in% enrset[enrset$padj <= pcutoff,"gset"] & dtd$padj <= pcutoff & dtd$lgfc > 0),c(5,4)]
    }
    
    createGraph <- function(df_data, ch_title) {
        gr <- simplify(graph.data.frame(df_data, directed=FALSE),remove.multiple=TRUE)
        v_weight <- vector()
        for (i in unique(df_data[,1])) {
            v_weight <- c(v_weight,length(unique(df_data[df_data[,1]==i,2])))
        }
        v_weight <- v_weight + 10
        v_weight <- c(v_weight,rep(1,length=length(unique(df_data[,2]))))
        
        plot(gr, 
             vertex.shape=c(rep("circle",length=length(unique(df_data[,1]))),rep("none",length=length(unique(df_data[,2])))),
             edge.arrow.size=1,
             vertex.size=v_weight,
             main=ch_title
        )
    }
    

    if (length(lst_p[,1])>0) {
        png(paste("./Figures/network_",gsetgroup,"_",contrast,"_up_",type,"_",pcutoff,".png",sep=""),width=1500,height=1500,res=300,pointsize=8)
        createGraph(df_data=lst_p, ch_title=paste(contrast, "up-regulated", sep=", "))
        dev.off()
    }
    
    
    if (length(lst_n[,1])>0) {
        png(paste("./Figures/network_",gsetgroup,"_",contrast,"_down_",type,"_",pcutoff,".png",sep=""),width=1500,height=1500,res=300,pointsize=8)
        createGraph(df_data=lst_n, ch_title=paste(contrast, "down-regulated", sep=", "))
        dev.off()
    }
    
}
