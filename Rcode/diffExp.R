diffExp <- function(arrayData, contrasts) {
    padj <- lgFC <- NA
	
    # Calculate p-values for each contrast:
	for (i in 1:length(contrasts)) {
	    a <- unlist(strsplit(contrasts[i], split=" - "))
	    
	    if (file.exists(paste("./Data/de_nbinom_",a[1],"_",a[2],".rds",sep=""))) {
	        res<- readRDS(paste("./Data/de_nbinom_",a[1],"_",a[2],".rds",sep=""))
	    } else {
	        res <- nbinomTest(arrayData$dataCount, a[1], a[2])
	        saveRDS(res,paste("./Data/de_nbinom_",a[1],"_",a[2],".rds",sep=""))
	    }
        
	    # Adjusted P-values and Fold changes for all contrasts
	    padj <- cbind(padj,res$padj)
	    lgFC <- cbind(lgFC,res$log2FoldChange)
	}

    padj <- as.data.frame(padj[,2:ncol(padj)],stringsAsFactors=FALSE)
    lgFC <- as.data.frame(lgFC[,2:ncol(lgFC)],stringsAsFactors=FALSE)
    
    rownames(padj) <- rownames(lgFC) <- res[,1]
	colnames(padj) <- colnames(lgFC) <- contrasts

	return(list(padj=padj,lgFC=lgFC))
}