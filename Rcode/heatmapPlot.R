heatmapPlot <- function(data,contrasts,gsc,GSAsign) {    
    rm(rownm_non,gpvalues_non,rownm_up,gpvalues_up,rownm_down,gpvalues_down)
    rownm_non <- rownm_up <- rownm_down<- vector()
    for (contrast in contrasts) {
        if (file.exists(paste("./Data/GSA_hyper_",goSubgroup,"_non_",contrast,"_",GSAsign,".rds",sep=""))) {
            gsah_non <- readRDS(paste("./Data/GSA_hyper_",goSubgroup,"_non_",contrast,"_",GSAsign,".rds",sep=""))
            #         gsah_up <- readRDS(paste("./Data/GSA_hyper_",goSubgroup,"_up_",contrast,"_",GSAsign,".rds",sep=""))
            #         gsah_down <- readRDS(paste("./Data/GSA_hyper_",goSubgroup,"_down_",contrast,"_",GSAsign,".rds",sep=""))
        } else {
            # non-regulated
            dt <- dexp$resTable[[contrast]]
            dt[!complete.cases(dt),c(2,5,6)] <- c(0,1,1)
            gsah_non <- runGSAhyper(genes = dt[,1], 
                                    pvalues = dt[,6], 
                                    pcutoff = GSAsign, 
                                    gsc = myGsc, 
                                    gsSizeLim=c(1,Inf), 
                                    adjMethod="fdr")
            saveRDS(gsah_non,paste("./Data/GSA_hyper_",goSubgroup,"_non_",contrast,"_",GSAsign,".rds",sep="")) # fromfile
            
            # up-regulated
            #         dt_up <- dt[which(dt$adj.P.Val  <= GSAsign & dt$logFC > 0),]
            #         gsah_up <- runGSAhyper(genes = dt_up[,1], 
            #                                pvalues = dt_up[,6], 
            #                                pcutoff = GSAsign, 
            #                                gsc = myGsc, 
            #                                gsSizeLim=c(1,Inf), 
            #                                adjMethod="fdr")
            #         saveRDS(gsah_up,paste("./Data/GSA_hyper_",goSubgroup,"_up_",contrast,"_",GSAsign,".rds",sep=""))
            
            # down-regulated
            #         dt_down <- dt[which(dt$adj.P.Val  <= GSAsign & dt$logFC <= 0),]
            #         gsah_down <- runGSAhyper(genes = dt_down[,1], 
            #                                  pvalues = dt_down[,6], 
            #                                  pcutoff = GSAsign, 
            #                                  gsc = myGsc, 
            #                                  gsSizeLim=c(1,Inf), 
            #                                  adjMethod="fdr")
            #         saveRDS(gsah_down,paste("./Data/GSA_hyper_",goSubgroup,"_down_",contrast,"_",GSAsign,".rds",sep=""))
            
            
        } # fromfile
        
        rownm_non <- c(rownm_non,names(gsah_non$p.adj[gsah_non$p.adj <= GSAsign]))
        #     rownm_up <- c(rownm_up,names(gsah_up$p.adj[gsah_up$p.adj <= GSAsign]))
        #     rownm_down <- c(rownm_down,names(gsah_down$p.adj[gsah_down$p.adj <= GSAsign]))
        if (exists("gpvalues_non")==TRUE) {
            gpvalues_non[,contrast] <- gsah_non$p.adj
            #         gpvalues_up[,contrast] <- gsah_up$p.adj
            #         gpvalues_down[,contrast] <- gsah_down$p.adj
        } else {
            gpvalues_non <- matrix(nrow=length(gsah_non$p.adj),ncol=length(contrasts))
            colnames(gpvalues_non) <- contrasts
            rownames(gpvalues_non) <- names(gsah_non$p.adj)
            gpvalues_non[,contrast] <- gsah_non$p.adj
            
            #         gpvalues_up <- matrix(nrow=length(gsah_up$p.adj),ncol=length(contrasts))
            #         colnames(gpvalues_up) <- contrasts
            #         rownames(gpvalues_up) <- names(gsah_up$p.adj)
            #         gpvalues_up[,contrast] <- gsah_up$p.adj
            #         
            #         gpvalues_down <- matrix(nrow=length(gsah_down$p.adj),ncol=length(contrasts))
            #         colnames(gpvalues_down) <- contrasts
            #         rownames(gpvalues_down) <- names(gsah_down$p.adj)
            #         gpvalues_down[,contrast] <- gsah_down$p.adj
        }
    }
    
    rownm_non <- sort(unique(rownm_non))
    gpdata_non <- gpvalues_non[rownames(gpvalues_non) %in% rownm_non,]
    # rownm_up <- sort(unique(rownm_up))
    # gpdata_up <- gpvalues_up[rownames(gpvalues_up) %in% rownm_up,]
    # rownm_down <- sort(unique(rownm_down))
    # gpdata_down <- gpvalues_down[rownames(gpvalues_down) %in% rownm_down,]
    
    
    # 
    par(mar=c(15,15,10,10),pin=c(10,15))
    dev.new()
    heatmap.2(as.matrix(gpdata_non), Rowv = FALSE, Colv=FALSE,margins=c(20,20), key=T, keysize=1, density.info="none", trace="none")
    # 
    # dev.new()
    # heatmap.2(gpdata_up, Rowv = FALSE, Colv=FALSE,margins=c(20,20), key=T, keysize=1, density.info="none", trace="none")
    
    # dev.new()
    # heatmap.2(gpdata_down, Rowv = FALSE, Colv=FALSE,margins=c(20,20), key=T, keysize=1, density.info="none", trace="none")
}