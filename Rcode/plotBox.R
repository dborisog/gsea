plotBox <- function(arrayData) {
    
    # Plot file:
    dev.new()
    par(mar=c(10, 4, 4, 2))
    boxplot(arrayData$dataNorm, main="Boxplot of normalized data",
            ylab=expression(log[2] ~ intensity), las=2)
    
}