exprDiff <- function(arrayData, contrasts) {
    name <- pval <- padj <- lgFC <- NA
    
    # Calculate p-values for each contrast:
    # see section 3.1 at page 6 http://bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf
    for (i in 1:length(contrasts)) {
        a <- unlist(strsplit(contrasts[i], split=" - "))
        
        print(paste("process",contrasts[i], Sys.time(),sep=" "))
        
        if (file.exists(paste("./Data/de_nbinom_",a[1],"_",a[2],".rds",sep=""))) {
            res<- readRDS(paste("./Data/de_nbinom_",a[1],"_",a[2],".rds",sep=""))
        } else {
            res <- nbinomTest(arrayData$dataCount, a[1], a[2])
            saveRDS(res,paste("./Data/de_nbinom_",a[1],"_",a[2],".rds",sep=""))
        }
        
        # P.Values, Adjusted P-values, and Fold changes for all contrasts
        pval <- cbind(pval,res$pval)
        padj <- cbind(padj,res$padj)
        lgFC <- cbind(lgFC,res$log2FoldChange)
    }
    
    name <- rownames(arrayData$dataNorm)
    pval <- as.data.frame(pval[,2:ncol(pval)],stringsAsFactors=FALSE)
    padj <- as.data.frame(padj[,2:ncol(padj)],stringsAsFactors=FALSE)
    lgFC <- as.data.frame(lgFC[,2:ncol(lgFC)],stringsAsFactors=FALSE)
    
    rownames(pval) <- rownames(padj) <- rownames(lgFC) <- res[,1]
    colnames(pval) <- colnames(padj) <- colnames(lgFC) <- contrasts
    
    return(list(name=name, pval=pval, padj=padj,lgFC=lgFC))
}