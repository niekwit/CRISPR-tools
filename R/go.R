library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(viridis)
#library(GOplot)

setwd("/run/user/1000/gvfs/smb-share:server=jcbc-store01.jcbc.private.cam.ac.uk,share=citiid_home/nw416/My Documents/Data/HDAC3 project/Results/Screen E data")

df <- read.csv("E6vs14.gene_summary copy.csv")

#get genes below FDR cut off
df.fdr <- df[df$neg.fdr < 0.25, ]
d <- df.fdr[,c("id","neg.lfc")]

#convert gene symbols to entrez gene IDs
d$id <- mapIds(org.Hs.eg.db, d$id,"ENTREZID","SYMBOL")

#remove rows with NA
d <- d[complete.cases(d$id),]

geneList <- d[ ,2]
names(geneList) <- as.character(d[,1])
geneList <- sort(geneList, decreasing = TRUE)


enrichPathway <- enrichPathway(d[ ,1],
                               organism = "human",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.2,
                               minGSSize = 10,
                               maxGSSize = 500,
                               readable = FALSE)
df <- enrichPathway@result
rownames(df) <- 1:nrow(df)
df <- df[1:10,]
df$log.pvalue <- -log10(df$p.adjust)

#plot enrichPathway results
p <- ggplot(df, aes(x=reorder(`Description`,`log.pvalue`),y=`log.pvalue`, fill=`Count`)) +
  geom_bar(stat='identity') + 
  coord_flip() + 
  labs(x='',y='-log(adj. p value)') +
  theme_bw(base_size = 20)+
  scale_fill_viridis()
  
ggsave(filename="enriched-pathways.pdf",
       p,
       width = 20,
       height= 5)





gsePathway <- gsePathway(geneList,
                         organism = "human",
                         exponent = 1,
                         minGSSize = 10,
                         maxGSSize = 500,
                         eps = 1e-10,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         verbose = TRUE,
                         seed = FALSE,
                         by = "fgsea")
