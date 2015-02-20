# This function allows to check if experiments are properly defined. 
# Outliers look very suspicious here.

plotPCA <- function(arrayData) {
            
    colors=c("red","green","blue","yellow","orange", "purple","tan","cyan","gray60","black", "white")
    # Centralize data
    dataForPCA <- l_da$dataNorm - rowMeans(l_da$dataNorm)    
    
    # PCA
    dataPrcomp <- prcomp(dataForPCA)
    pc1 <- dataPrcomp$rotation[,1]
    pc2 <- dataPrcomp$rotation[,2]
    clr <- 1:length(pc1)
    
    # PCA plot
    
    png(paste("./Figures/pcaplot.png",sep=""),width=1500,height=1500,res=300,pointsize=8)
    par(mar=c(10, 4, 4, 2))
    par(mar=c(5,5,4,4))
    plot(cbind(pc1,pc2), 
         pch=15, 
         col=colors,
         cex=1.5,
         main="PCA plot, normalized dataset", 
         xlab=paste("PC1 (",round(summary(dataPrcomp)$importance[2,1]*100,digits=2),"%)",sep=""), 
         ylab=paste("PC2 (",round(summary(dataPrcomp)$importance[2,2]*100,digits=2),"%)",sep=""),
         xlim=c(min(pc1)-0.05,max(pc1)+0.05),
         ylim=c(min(pc2)-0.05,max(pc2)+0.05))
    text(cbind(pc1,pc2), labels=names(pc1), pos=4, cex=0.5)
    dev.off()
}