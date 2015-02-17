plotVenn <- function(data, transcript2gene, contrasts, pvalue=0.05, adj.pvalue=0.1) {
    
    transcript2gene <- transcript2gene[,1:2]; colnames(transcript2gene) <- c("ENSEMBLTRANS", "ENTREZID")
    lng <- length(contrasts)+1
    v_contr <- vector() 
    for (contrast in contrasts) {
        a <- unlist(strsplit(contrast, split=" - "))
        v_contr <- append(v_contr, a[2])
    }
    
    
    m.tnp <- matrix(nrow=length(data$pval[[1]]),ncol=length(contrasts))
    m.tna <- matrix(nrow=length(data$padj[[1]]),ncol=length(contrasts))
    colnames(m.tnp) <- colnames(m.tna) <- contrasts
    rownames(m.tnp) <- rownames(m.tna) <- rownames(data$pval)
    for (contrast in contrasts) {
        m.tnp[,contrast] <- data$pval[[contrast]] <= pvalue
        m.tna[,contrast] <- data$padj[[contrast]] <= adj.pvalue
        # according to the info on FBgn0000014 I assume that NA can be replaced with 1
        # http://bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf
        m.tnp[is.na(m.tnp[,contrast]),contrast] <- FALSE
        m.tna[is.na(m.tna[,contrast]),contrast] <- FALSE
    }
    
    df.tnp <- as.data.frame(m.tnp)
    df.tnp$ENSEMBLTRANS <- rownames(m.tnp)
    df.tnp <- merge(df.tnp,transcript2gene,by="ENSEMBLTRANS",all.X=TRUE)[,c(2:(length(contrasts)+2))]
    
    df.tna <- as.data.frame(m.tna)
    df.tna$ENSEMBLTRANS <- rownames(m.tna)
    df.tna <- merge(df.tna,transcript2gene,by="ENSEMBLTRANS",all.X=TRUE)[,c(2:(length(contrasts)+2))]
    
    #### tranform transcripts to genes
    # details are here 
    # http://stackoverflow.com/questions/27483609/any-intelligent-solution-to-transform-matrix-of-transcripts-to-a-matrix-of-genes
    m.gnp <- setkey(setDT(df.tnp), ENTREZID)[, lapply(.SD, sum), ENTREZID]
    m.gnp <- as.matrix(as.data.frame(m.gnp)[,2:lng])
    
    
    m.gna <- setkey(setDT(df.tna), ENTREZID)[, lapply(.SD, sum), ENTREZID]
    m.gna <- as.matrix(as.data.frame(m.gna)[,2:lng])
    
    ####
    #### create venn diagram
    nm <- LETTERS[seq(from = 1, to = length(v_contr))]
    tx <- paste(nm, ': ',v_contr, '; ', sep="", collapse="")
    ps <- gregexpr("; ", tx)[[1]][2]
    substr(tx,start=ps, stop=(ps+1)) <- "\n"
    png(paste("./Figures/venn_all_",v_contr[1],".png",sep=""),width=2100,height=2100,res=300,pointsize=8)
    par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
    vennDiagram(m.tnp,names=nm,cex=1,main=paste("\n\nTranscripts, p-value <= ",pvalue,sep=""),cex.main=1)
    vennDiagram(m.gnp,names=nm,cex=1,main=paste("\n\nGenes, p-value <= ",pvalue,sep=""),cex.main=1)
    vennDiagram(m.tna,names=nm,cex=1,main=paste("\n\nTranscripts, FDR adjusted p-value <= ",adj.pvalue,sep=""),cex.main=1)
    vennDiagram(m.gna,names=nm,cex=1,main=paste("\n\nGenes, FDR adjusted p-value <= ",adj.pvalue,sep=""),cex.main=1)
    par(mfrow=c(1,1), mar=c(0,3,1,1))
    mtext("Significant transcripts and genes.\n\n",side=2,cex=1.3)
    mtext(tx,side=1, cex=0.95)
    dev.off()
    
    
    #######################
    #######################
    # up-regulated
    m.tupp <- m.tupa <-  matrix(nrow=length(data$pval[[1]]),ncol=length(contrasts))
    colnames(m.tupp) <- colnames(m.tupa) <- contrasts
    rownames(m.tupp) <- rownames(m.tupa) <- rownames(data$pval)
    for (contrast in contrasts) {
        m.tupp[,contrast] <- (data$pval[[contrast]] <= pvalue & data$lgFC[[contrast]] > 0)
        m.tupa[,contrast] <- (data$padj[[contrast]] <= adj.pvalue & data$lgFC[[contrast]] > 0)
        # according to the info on FBgn0000014 I assume that NA can be replaced with 1
        # http://bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf
        m.tupp[is.na(m.tupp[,contrast]),contrast] <- FALSE
        m.tupa[is.na(m.tupa[,contrast]),contrast] <- FALSE
    }
    
    df.tupp <- as.data.frame(m.tupp)
    df.tupp$ENSEMBLTRANS <- rownames(m.tupp)
    df.tupp <- merge(df.tupp,transcript2gene,by="ENSEMBLTRANS",all.X=TRUE)[,c(2:(length(contrasts)+2))]
    
    df.tupa <- as.data.frame(m.tupa)
    df.tupa$ENSEMBLTRANS <- rownames(m.tupa)
    df.tupa <- merge(df.tupa,transcript2gene,by="ENSEMBLTRANS",all.X=TRUE)[,c(2:(length(contrasts)+2))]
    
    #### tranform transcripts to genes
    # details are here 
    # http://stackoverflow.com/questions/27483609/any-intelligent-solution-to-transform-matrix-of-transcripts-to-a-matrix-of-genes
    m.gupp <- setkey(setDT(df.tupp), ENTREZID)[, lapply(.SD, sum), ENTREZID]
    m.gupp <- as.matrix(as.data.frame(m.gupp)[,2:lng])
    
    
    m.gupa <- setkey(setDT(df.tupa), ENTREZID)[, lapply(.SD, sum), ENTREZID]
    m.gupa <- as.matrix(as.data.frame(m.gupa)[,2:lng])
    
    ####
    #### create venn diagram
    nm <- LETTERS[seq(from = 1, to = length(v_contr))]
    tx <- paste(nm, ': ',v_contr, '; ', sep="", collapse="")
    ps <- gregexpr("; ", tx)[[1]][2]
    substr(tx,start=ps, stop=(ps+1)) <- "\n"
    png(paste("./Figures/venn_up_",v_contr[1],".png",sep=""),width=2100,height=2100,res=300,pointsize=8)
    par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
    vennDiagram(m.tupp,names=nm,cex=1,main=paste("\n\nTranscripts, p-value <= ",pvalue,sep=""),cex.main=1)
    vennDiagram(m.gupp,names=nm,cex=1,main=paste("\n\nGenes, p-value <= ",pvalue,sep=""),cex.main=1)
    vennDiagram(m.tupa,names=nm,cex=1,main=paste("\n\nTranscripts, FDR adjusted p-value <= ",adj.pvalue,sep=""),cex.main=1)
    vennDiagram(m.gupa,names=nm,cex=1,main=paste("\n\nGenes, FDR adjusted p-value <= ",adj.pvalue,sep=""),cex.main=1)
    par(mfrow=c(1,1), mar=c(0,3,1,1))
    mtext("Up-regulated transcripts and genes.\n\n",side=2,cex=1.3)
    mtext(tx,side=1, cex=0.95)
    dev.off()
    
    #######################
    #######################
    # down-regulated
    m.tdwp <- m.tdwa <-  matrix(nrow=length(data$pval[[1]]),ncol=length(contrasts))
    colnames(m.tdwp) <- colnames(m.tdwa) <- contrasts
    rownames(m.tdwp) <- rownames(m.tdwa) <- rownames(data$pval)
    for (contrast in contrasts) {
        m.tdwp[,contrast] <- (data$pval[[contrast]] <= pvalue & data$lgFC[[contrast]] <= 0)
        m.tdwa[,contrast] <- (data$padj[[contrast]] <= adj.pvalue & data$lgFC[[contrast]] <= 0)
        # according to the info on FBgn0000014 I assume that NA can be replaced with 1
        # http://bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf
        m.tdwp[is.na(m.tdwp[,contrast]),contrast] <- FALSE
        m.tdwa[is.na(m.tdwa[,contrast]),contrast] <- FALSE
    }
    
    df.tdwp <- as.data.frame(m.tdwp)
    df.tdwp$ENSEMBLTRANS <- rownames(m.tdwp)
    df.tdwp <- merge(df.tdwp,transcript2gene,by="ENSEMBLTRANS",all.X=TRUE)[,c(2:(length(contrasts)+2))]
    
    df.tdwa <- as.data.frame(m.tdwa)
    df.tdwa$ENSEMBLTRANS <- rownames(m.tdwa)
    df.tdwa <- merge(df.tdwa,transcript2gene,by="ENSEMBLTRANS",all.X=TRUE)[,c(2:(length(contrasts)+2))]
    
    #### tranform transcripts to genes
    # details are here 
    # http://stackoverflow.com/questions/27483609/any-intelligent-solution-to-transform-matrix-of-transcripts-to-a-matrix-of-genes
    m.gdwp <- setkey(setDT(df.tdwp), ENTREZID)[, lapply(.SD, sum), ENTREZID]
    m.gdwp <- as.matrix(as.data.frame(m.gdwp)[,2:lng])
    
    
    m.gdwa <- setkey(setDT(df.tdwa), ENTREZID)[, lapply(.SD, sum), ENTREZID]
    m.gdwa <- as.matrix(as.data.frame(m.gdwa)[,2:lng])
    
    ####
    #### create venn diagram
    nm <- LETTERS[seq(from = 1, to = length(v_contr))]
    tx <- paste(nm, ': ',v_contr, '; ', sep="", collapse="")
    ps <- gregexpr("; ", tx)[[1]][2]
    substr(tx,start=ps, stop=(ps+1)) <- "\n"
    png(paste("./Figures/venn_down_",v_contr[1],".png",sep=""),width=2100,height=2100,res=300,pointsize=8)
    par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
    vennDiagram(m.tdwp,names=nm,cex=1,main=paste("\n\nTranscripts, p-value <= ",pvalue,sep=""),cex.main=1)
    vennDiagram(m.gdwp,names=nm,cex=1,main=paste("\n\nGenes, p-value <= ",pvalue,sep=""),cex.main=1)
    vennDiagram(m.tdwa,names=nm,cex=1,main=paste("\n\nTranscripts, FDR adjusted p-value <= ",adj.pvalue,sep=""),cex.main=1)
    vennDiagram(m.gdwa,names=nm,cex=1,main=paste("\n\nGenes, FDR adjusted p-value <= ",adj.pvalue,sep=""),cex.main=1)
    par(mfrow=c(1,1), mar=c(0,3,1,1))
    mtext("Down-regulated transcripts and genes.\n\n",side=2,cex=1.3)
    mtext(tx,side=1, cex=0.95)
    dev.off()
}










