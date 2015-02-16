plotHist <- function(arrayData) {
    
    dev.new()
    plot(density(arrayData$dataNorm[,1]), main="Histogram of normalized data",
         xlab=expression(log[2] ~ intensity),ylab="Density")
    for(i in 2:dim(arrayData$dataNorm)[2]){
        lines(density(arrayData$dataNorm[,i]))
    }
    
}