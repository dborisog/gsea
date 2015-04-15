plotManualNetwork <- function(dexp, contrast, gsets, enrsets, pcutoff, type, gsetgroup, shownodes) {
    # shownodes are either "all" sets and genes or "sets"
    # shownodes="all"
    
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
    
    
    dtd <- merge(dt, gsets,by="FBtranscriptID")
    
    
    if (type=="pval"){
        lst_n <- dtd[which(dtd$GSET %in% enrsets[enrsets$padj <= pcutoff,"gset"] & dtd$padj <= pcutoff & dtd$lgfc < 0),c(5,6)]
        lst_p <- dtd[which(dtd$GSET %in% enrsets[enrsets$padj <= pcutoff,"gset"] & dtd$padj <= pcutoff & dtd$lgfc > 0),c(5,6)]
    } else {
        lst_n <- dtd[which(dtd$GSET %in% enrsets[enrsets$padj <= pcutoff,"gset"] & dtd$padj <= pcutoff & dtd$lgfc < 0),c(5,6)]
        lst_p <- dtd[which(dtd$GSET %in% enrsets[enrsets$padj <= pcutoff,"gset"] & dtd$padj <= pcutoff & dtd$lgfc > 0),c(5,6)]
    }
    
    createGraph <- function(df_data, ch_title, shownodes) {
        if (shownodes == "all") {
            gr <- simplify(graph.data.frame(df_data, directed=FALSE),remove.multiple=TRUE)
            v_weight <- vector()
            for (i in unique(df_data[,1])) {
                v_weight <- c(v_weight,length(unique(df_data[df_data[,1]==i,2])))
            }
            v_weight <- v_weight + 10
            v_weight <- c(v_weight,rep(1,length=length(unique(df_data[,2]))))
            V(gr)$shape <- c(rep("circle",length=length(unique(df_data[,1]))),rep("none",length=length(unique(df_data[,2]))))
            V(gr)$size <- v_weight
            V(gr)$label.color <- c(rep("black",length=length(unique(df_data[,1]))),rep("red",length=length(unique(df_data[,2]))))
            
            
            id <- tkplot(gr, canvas.width=1000, canvas.height=1000)
            canvas <- tkplot.canvas(id)
            width <- as.numeric(tkcget(canvas, "-width"))
            height <- as.numeric(tkcget(canvas, "-height"))
            tkcreate(canvas, "text", width/2, 25, text=ch_title,justify="center", font=tkfont.create(family="helvetica",size=20,weight="bold"))
        } else {
            mlst <- merge(df_data,df_data,by="GNAME"); mlst <- mlst[,2:3]
            am <- get.adjacency(graph.edgelist(as.matrix(mlst)))
            gr <- simplify(graph.adjacency(am, mode="upper"),remove.multiple=TRUE)
            
            id <- tkplot(gr, canvas.width=1000, canvas.height=1000)
            canvas <- tkplot.canvas(id)
            width <- as.numeric(tkcget(canvas, "-width"))
            height <- as.numeric(tkcget(canvas, "-height"))
            tkcreate(canvas, "text", width/2, 25, text=ch_title,justify="center", font=tkfont.create(family="helvetica",size=20,weight="bold"))
        }
    }
    
    
    if (length(lst_p[,1])>0) {
        createGraph(df_data=lst_p, ch_title=paste(contrast, "up-regulated", sep=", "),shownodes=shownodes)
    }
    
    
    if (length(lst_n[,1])>0) {
        createGraph(df_data=lst_n, ch_title=paste(contrast, "down-regulated", sep=", "),shownodes=shownodes)
    }
    
}