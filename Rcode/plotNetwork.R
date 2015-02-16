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
    
    mlst_p <- merge(lst_p,lst_p,by="ENTREZID")
    mlst_n <- merge(lst_n,lst_n,by="ENTREZID")
    
    
    am_p <- get.adjacency(graph.edgelist(as.matrix(mlst_p[,2:3]), directed=FALSE), type="upper", edges=FALSE)
    am_n <- get.adjacency(graph.edgelist(as.matrix(mlst_n[,2:3]), directed=FALSE), type="upper", edges=FALSE)
    
    gr_p <- simplify(graph.adjacency(am_p, mode="upper",weighted=NULL,diag=FALSE),remove.multiple=TRUE)
    gr_n <- simplify(graph.adjacency(am_n, mode="upper",weighted=NULL,diag=FALSE),remove.multiple=TRUE)
    

    if (length(lst_p[,1])>0) {
        png(paste("./Figures/network_",gsetgroup,"_",contrast,"_up_",type,"_",pcutoff,".png",sep=""),width=1500,height=1500,res=300,pointsize=8)
        plot.igraph(gr_p, edge.arrow.size=1, vertex.size=as.matrix(am_p)+10, layout=layout.circle,main=paste(contrast, "up-regulated", sep=", "))
        dev.off()
    }
    
    
    if (length(lst_n[,1])>0) {
        png(paste("./Figures/network_",gsetgroup,"_",contrast,"_down_",type,"_",pcutoff,".png",sep=""),width=1500,height=1500,res=300,pointsize=8)
        plot.igraph(gr_n, edge.arrow.size=1, vertex.size=as.matrix(am_n)+10, layout=layout.circle,main=paste(contrast, "down-regulated", sep=", "))
        dev.off()
    }
    
}


