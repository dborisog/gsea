enrichSets <- function(dexp, pvalues, sets, pcutoff) {
    pvalues[pvalues==0] <- 1e-10
    dt <- data.frame("FBtranscriptID"=dexp$name, "padj"=pvalues)
    dt$FBtranscriptID <- as.character(dt$FBtranscriptID)
    dt[is.na(dt$padj),2] <- 1; sets[is.na(sets$ENTREZID)==TRUE,2] <- "unknown"; sets[is.na(sets$GSET)==TRUE,3] <- "unknown"
    
    dtd <- merge(dt, sets,by="FBtranscriptID")

    
    dtd$strg <- dtd$padj <= pcutoff; dtd$weak <- dtd$padj > pcutoff
    strg <- dtd[,4:5]; weak <- dtd[,c(4,6)]

    strgs <- setkey(setDT(strg), GSET)[, lapply(.SD, sum), GSET] # calculate No of strong genes in each gene set
    weaks <- setkey(setDT(weak), GSET)[, lapply(.SD, sum), GSET] # calculate No of weak genes in each gene set
    
    #data for fisher.test
    ftdt <- data.frame("strg" = strgs$strg,                   # count of strong genes in gene set
                       "slft" = sum(strgs$strg) - strgs$strg, # count of strong genes left
                       "weak" = weaks$weak,                   # count of weak genes in gene set
                       "wlft" = sum(weaks$weak) - weaks$weak) # count of weak genes left
    
    
    pval <- padj <- NA # calculate p.value and adjusted p.value for gene sets
    for(i in 1:length(ftdt[,1])) {
        pval[i] <- fisher.test(matrix(as.integer(ftdt[i,]),2,2),alternative="greater")$p.value
    }
    padj <- p.adjust(pval, method="fdr")
    
    
    res <- data.frame("gset" = strgs$GSET, "pval" = pval, "padj" = padj); res$gset <- as.character(res$gset)
    return(res)
}
# http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html
# We are asking if diamond enriched or depleted in our 5 selected cards, comparing to the background.
#
#             Diamond    Non-Diamond
# selected          x            5-x  total 5 sampled cards
# left           13-x           34+x  total 47 left cards after sampling
#              13 Dia     39 Non-Dia  total 52 cards
# 
# phyper(x, 13, 39, 5, lower.tail=TRUE)
# Numerical parameters in order
# (success-in-sample, success-in-bkgd, failure-in-bkgd, sample-size)
# 
# fisher.test(matrix(c(x, 13-x, 5-x, 34+x), 2, 2), alternative='less');
# Numerical parameters in order:
# (success-in-sample, success-in-left-part, failure-in-sample, failure-in-left-part).