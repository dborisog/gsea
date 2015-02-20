plotBox <- function(arrayData) {
    
    # Plot file:
    png(paste("./Figures/boxplot.png",sep=""),width=1500,height=1500,res=300,pointsize=8)
    par(mar=c(10, 4, 4, 2))
    boxplot(arrayData$dataNorm, main="Boxplot of normalized data", ylab=expression(log[2] ~ intensity), las=2)
    dev.off()
    
}