library(rJava)
library("RDAVIDWebService")
library("org.Hs.eg.db")
library(biomaRt)
library(dplyr)
library(GOplot)
library(stringr)
library(ggrepel)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
#email <- args[1]
#mageck <- args[2]
#fdr <- args[3]
#species <- args[4]
fdr <- 0.25
email <- "nw416@cam.ac.uk"

setwd("/run/user/1000/gvfs/smb-share:server=jcbc-store01.jcbc.private.cam.ac.uk,share=citiid_grp_nathan/Niek/LMB_data/Results/CRISPR SCREENING/NGS data and analysis/Screen E SHMT1+2 HDAC3 Ca9 HT Bassik/data-analysis-screen-e/mageck/E6vsE14")
df <- read.csv(file="E6vs14.gene_summary.txt", sep="\t")

#get ensemble gene IDs
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
df.temp <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"),
                 filters="hgnc_symbol",
                 values=df$id,
                 mart=ensembl)
df.temp <- df.temp %>%
  rename('id' = 'hgnc_symbol')
df <- merge(x=df,y=df.temp,by="id")

#select gene below fdr cut off
df.fdr <- df[df$neg.fdr < fdr, ]

gene.list <- df.fdr$ensembl_gene_id

david <- DAVIDWebService(email=email, url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

result <- addList(david, gene.list,idType="ENSEMBL_GENE_ID",listName="gene.list", listType="Gene")

setAnnotationCategories(david, "GOTERM_BP_ALL")

termCluster <- getClusterReport(david, type="Term")

getClusterReportFile(david, type="Term", fileName="termClusterReport1.tab")

davidGODag <- DAVIDGODag(members(termCluster) [[1]], pvalueCutoff=0.1,"BP")

pdf("GO-TermGraph.pdf")
plotGOTermGraph(g=goDag(davidGODag),r=davidGODag, max.nchar=30, node.shape="ellipse")
dev.off()

#extract most enriched GO from each cluster
cluster.total <- nrow(summary(termCluster))

for (i in 1:cluster.total){
  df.test <- termCluster@cluster[[i]]$Members
  df.test <- as.data.frame(do.call(cbind, df.test))
  cols.num <- c("Count","X.","PValue","Fold.Enrichment","Bonferroni","Benjamini","FDR")
  df.test[cols.num] <- sapply(df.test[cols.num], as.numeric)
  
  if(i == 1){df.go <- df.test[df.test$Fold.Enrichment == max(df.test$Fold.Enrichment),]}
  else {df.go <- rbind(df.go, df.test[df.test$Fold.Enrichment == max(df.test$Fold.Enrichment),])}
}

goplot.gene.list <- df[ ,c("ensembl_gene_id","neg.lfc")]
colnames(goplot.gene.list) <- c("ID","logFC")

df.go$Category <- "BP"
temp.list <- str_split_fixed(df.go$Term, "~", 2)
df.go <- cbind(df.go, temp.list)
df.go <- rename(df.go, "Term2"="Term")
df.go <- rename(df.go, c("ID"="1", "Term"="2"))
goplot.david <- df.go[ ,c("Category","ID","Term","Genes","Benjamini")]
goplot.david <- rename(goplot.david, "adj_pval"="Benjamini")

circ <- circle_dat(goplot.david, goplot.gene.list)
circ$log.pvalue <- -log(circ$adj_pval)
circ$term <- as.factor(circ$term)
#circ$log.pvalue <- as.factor(circ$log.pvalue)
#circ$count <- as.factor(circ$count)
circ.plot <- subset(circ, select=-c(genes,logFC))
circ.plot <- distinct(circ.plot)

circ.label <- circ.plot[circ.plot$log.pvalue > 3, ]
circ.table <- circ.label[ ,c("ID","term")]
table1 <- tableGrob(circ.table, rows=NULL, theme=ttheme_minimal())

p <- ggplot(circ.plot, aes(x = `zscore`, 
                    y = `log.pvalue`,
                    size= `count`)) +
  scale_size(range = c(5, 15)) +
  theme_bw(base_size = 22) + 
  theme(legend.position = "top") +
  geom_point(shape=21,
             alpha=0.5,
             fill="red") +
  guides(fill= FALSE) +
  xlab("z-score") +
  ylab("-log10(adj. p value)") +
  geom_text_repel(aes(x=`zscore`,y=`log.pvalue`,label=`ID`, size = 3),
                  data=circ.label,
                  nudge_y=c(-2,2)) +
  geom_hline(yintercept=3, linetype="dashed", color = "black") 

grid.arrange(p,table1,nrow=1, ncol=2, as.table=TRUE)


