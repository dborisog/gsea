# ==================================
# ==================================
#  if (group=="rfbr_orthologs") {} at loadGSC
# ==================================


# read from manually prepared data files
# 'dm' for Drosophila Melanogaster, 'ce' for Caenorhabditis elegans, 'mm' for Mus musculus; 
# 'G' for gene, 'O' for ortholog, 'lg' for longevity gene;
# 'dmG' read as genes of 'dm';
# 'dmOce' othologs of/to 'dm' derived from 'ce'.
dmG_lg <- read.delim(file="./Data/dmG_lg.csv", header=TRUE)
dmOce_lg <- read.delim(file="./Data/dmOce_lg.csv", header=TRUE)
dmOmm_lg <- read.delim(file="./Data/dmOmm_lg.csv", header=TRUE)
ceG_lg <- read.delim(file="./Data/ceG_lg.csv", header=TRUE)
mmG_lg <- read.delim(file="./Data/mmG_lg.csv", header=TRUE)

# remove unnecessary strings and transform all to characters
dmG_lg <- data.frame(lapply(dmG_lg[,c(1,3,11)], as.character), stringsAsFactors=FALSE)
dmOce_lg <- data.frame(lapply(dmOce_lg[,c(1,4)], as.character), stringsAsFactors=FALSE)
dmOmm_lg <- data.frame(lapply(dmOmm_lg[,c(1,4)], as.character), stringsAsFactors=FALSE)
ceG_lg <- data.frame(lapply(ceG_lg[,c(1,3,5)], as.character), stringsAsFactors=FALSE)
mmG_lg <- data.frame(lapply(mmG_lg[,c(1,3,5)], as.character), stringsAsFactors=FALSE)

# transform colnames to the same
# c("ENTREZID","Name", "Influence") and c("ENTREZID","Name")
colnames(dmG_lg) <- c("ENTREZID","Name", "Influence")
colnames(ceG_lg) <- c("ENTREZID","Name", "Influence")
colnames(mmG_lg) <- c("ENTREZID","Name", "Influence")
colnames(dmOce_lg) <- c("ENTREZID","Name")
colnames(dmOmm_lg) <- c("ENTREZID","Name")

# clean ENTREZID from prefix
library("stringr")
dmOce_lg$ENTREZID <- str_trim(str_replace(dmOce_lg$ENTREZID, 'ENTREZGENE_ACC:', ''))
dmOmm_lg$ENTREZID <- str_trim(str_replace(dmOmm_lg$ENTREZID, 'ENTREZGENE_ACC:', ''))

# remove "N/A" from dmOce and dmOmm
dmOce_lg <- dmOce_lg[which(dmOce_lg$Name != "N/A"),]
dmOmm_lg <- dmOmm_lg[which(dmOmm_lg$Name != "N/A"),]

# add longevity influences to dmOce_lg and dmOnn_lg
dmOce_lg <- merge(dmOce_lg,ceG_lg[,c(1,3)],by="ENTREZID")
dmOmm_lg <- merge(dmOmm_lg,mmG_lg[,c(1,3)],by="ENTREZID")

# define classes
# > unique(dmG_lg$Influence)
# [1] "anti" "pro"  "none"
# > unique(dmOce_lg$Influence)
# [1] "Anti-Longevity" "Pro-Longevity" "Unclear" "Unannotated"   
# > unique(dmOmm_lg$Influence)
# [1] "Pro-Longevity"  "Anti-Longevity"
#### I should keep all these classes -- they might be relevant in some research. 
#### However, in other research some classes might be removed

# create classes
dmG_lg$Influence <- paste("dmG", dmG_lg$Influence, sep="_")

dmOce_lg$Influence <- str_trim(str_replace(dmOce_lg$Influence, 'Anti-Longevity', 'dmOce_anti'))
dmOce_lg$Influence <- str_trim(str_replace(dmOce_lg$Influence, 'Pro-Longevity', 'dmOce_pro'))
dmOce_lg$Influence <- str_trim(str_replace(dmOce_lg$Influence, 'Unclear', 'dmOce_unclear'))
dmOce_lg$Influence <- str_trim(str_replace(dmOce_lg$Influence, 'Unannotated', 'dmOce_unannotated'))

dmOmm_lg$Influence <- str_trim(str_replace(dmOmm_lg$Influence, 'Anti-Longevity', 'dmOmm_anti'))
dmOmm_lg$Influence <- str_trim(str_replace(dmOmm_lg$Influence, 'Pro-Longevity', 'dmOmm_pro'))

# merge rows
dm_lg <- rbind(dmG_lg[,c(1,3)], dmOce_lg[,c(1,3)],dmOmm_lg[,c(1,3)])

# Flybase transcript to ENTREZID
# at the end I should have F.Transcript in one column, and Influence in another
### There is a chance that this table might actually be more comlete if I would utilize multiple sources
# discover & select
keytypes(org.Dm.eg.db)
columns(org.Dm.eg.db)
fbids <- dm_lg$ENTREZID
cols <- c("ENSEMBLTRANS", "SYMBOL")
anno <- select(org.Dm.eg.db, keys = fbids, columns = cols, keytype = "ENTREZID")
dmT_lg <- merge(anno,dm_lg, by="ENTREZID")

saveRDS(dmT_lg,"./Data/dmT_lg.rds")