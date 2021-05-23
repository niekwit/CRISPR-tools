#check packages
cran.packages <- c("BiocManager","ggplot2","yaml","ggrepel",
                   "viridis","dplyr","gridExtra","GOplot",
                   "stringr","cowplot","rJava")
installed.packages <- installed.packages()[,1]

#install CRAN packages
for (i in cran.packages){
  if(!i %in% installed.packages){
    cat(i,"package missing, installing now\n")
    install.packages(i,repos='http://cran.us.r-project.org')
  }
}

#install BiocManager packages
bioc.packages <- c("RDAVIDWebService", "biomaRt")

for (i in bioc.packages){
  if(!i %in% installed.packages){
    cat(i,"package missing, installing now\n")
    BiocManager::install(i, force = TRUE)
  }
}

library(rJava)
library("RDAVIDWebService")
library(biomaRt)
library(dplyr)
library(GOplot)
library(stringr)
library(ggrepel)
library(ggplot2)
library(gridExtra)
library(BiocManager)

args <- commandArgs(trailingOnly=TRUE)
email <- args[1]
mageck.file <- args[2]
fdr <- args[3]
species <- args[4]
save.path <- args[5]
go.test <- args[6]
go.term <- args[7]

df <- read.csv(file=mageck.file, sep="\t")

#get ensemble gene IDs (needed for DAVID)
if(species == "human"){ensembl.species <- "hsapiens_gene_ensembl"
} else if (species == "mouse"){ensembl.species <- "mmusculus_gene_ensembl"}

ensembl <- useMart("ensembl")
ensembl <- useDataset(ensembl.species,mart=ensembl)
df.temp <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"),
                 filters="hgnc_symbol",
                 values=df$id,
                 mart=ensembl)
df.temp <- df.temp %>%
  rename('id' = 'hgnc_symbol')
df <- merge(x=df,y=df.temp,by="id")

#function to perform GO analysis
go.analysis <- function(x){
  if(x == "neg"){
    z <- "neg.fdr"
    title <- "Drop out gene GO analysis"
    file.plot <- paste0(save.path,"/GO-analysis-dropout.pdf")
  } else if(x == "pos") {
    z <- "pos.fdr"
    title <- "Enrichment gene GO analysis"
    file.plot <- paste0(save.path,"/GO-analysis-enrichment.pdf")
  }
  
  if(go.term == "BP"){
    go.term <- "GOTERM_BP_ALL"
  } else if(go.term == "CC"){
      go.term <- "GOTERM_CC_ALL"
  } else if(go.term == "MF"){
        go.term <- "GOTERM_MF_ALL"}
  
  #select genes below FDR cutoff
  df.fdr <- df[df[z] < fdr, ]
  gene.list <- df.fdr$ensembl_gene_id
  
  #perform DAVID
  david <- DAVIDWebService(email=email, url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  result <- addList(david, gene.list,idType="ENSEMBL_GENE_ID",listName="gene.list", listType="Gene")
  setAnnotationCategories(david, go.term)
  termCluster <- getClusterReport(david, type="Term")
  getClusterReportFile(david, type="Term", fileName=paste0(save.path,"/GO-DAVID-cluster-output.tab"))
  
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
  
  #calculate z-scores with GOplot
  goplot.gene.list <- df[ ,c("ensembl_gene_id","neg.lfc")]
  colnames(goplot.gene.list) <- c("ID","logFC")
  df.go$Category <- "BP"
  temp.list <- str_split_fixed(df.go$Term, "~", 2)
  df.go <- cbind(df.go, temp.list)
  df.go <- rename(df.go, "Term2"="Term")
  df.go <- rename(df.go, c("ID"="1", "Term"="2"))
  goplot.david <- df.go[ ,c("Category","ID","Term","Genes")]
  goplot.david[go.test] <- df.go[go.test]
  goplot.david <- rename(goplot.david, "adj_pval"=`go.test`)
  circ <- circle_dat(goplot.david, goplot.gene.list)
  circ$log.pvalue <- -log(circ$adj_pval)
  circ$term <- as.factor(circ$term)
  
  #create df for plotting
  circ.plot <- subset(circ, select=-c(genes,logFC))
  circ.plot <- distinct(circ.plot)
  
  #create df for plot labels
  circ.label <- circ.plot[circ.plot$log.pvalue > 3, ]
  
  #create table for plotting
  circ.table <- circ.label[ ,c("ID","term")]
  table1 <- tableGrob(circ.table, rows=NULL, theme=ttheme_minimal(base_size=16))
  
  #create plot
  p <- ggplot(circ.plot, aes(x = `zscore`, 
                             y = `log.pvalue`,
                             size= `count`)) +
    scale_size(range = c(5, 15), name = "Number of genes in group:",
               guide = guide_legend(
                 direction = "horizontal",
                 title.position = "top")) +
    theme_bw(base_size = 22) + 
    theme(legend.justification= "top",legend.direction="horizontal") +
    geom_point(shape=21,
               alpha=0.5,
               fill="red") +
    guides(fill= FALSE) +
    xlab("z-score") +
    ylab("-log10(adj. p value)") +
    geom_text_repel(aes(x=`zscore`,y=`log.pvalue`,label=`ID`, size = 3),
                    data=circ.label,
                    nudge_y=c(-1,1),
                    nudge_x=c(0,0.5),
                    show.legend = FALSE) +
    geom_hline(yintercept=3, linetype="dashed", color = "black") +
    ggtitle(title)
  
  legend <- cowplot::get_legend(p)
  
  p <- grid.arrange(p + guides(size=FALSE),
                    legend,
                    table1,
                    layout_matrix = rbind(c(1,2),
                                          c(1,3),
                                          c(1,3),
                                          c(1,3),
                                          c(1,3)))
  
  #save plot
  ggsave(p,file=file.plot, width=15, height=7)
}

#perform GO analysis
input <- c("neg","pos")
lapply(input, go.analysis)