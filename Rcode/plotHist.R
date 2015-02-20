plotHist <- function(arrayData) {
    
    png(paste("./Figures/histplot.png",sep=""),width=1500,height=1500,res=300,pointsize=8)
    par(mar=c(10, 4, 4, 2))
    plot(density(arrayData$dataNorm[,1]), main="Histogram of normalized data",
         xlab=expression(log[2] ~ intensity),ylab="Density")
    for(i in 2:dim(arrayData$dataNorm)[2]){
        lines(density(arrayData$dataNorm[,i]))
    }
    dev.off()

    
}