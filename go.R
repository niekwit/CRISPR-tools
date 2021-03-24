library(rJava)
library("RDAVIDWebService")
library(biomaRt)

args <- commandArgs(trailingOnly=TRUE)
email <- args[1]

david <- DAVIDWebService(email=email, url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

setwd("/home/niek/Documents/analyses/CRISPR-screens/tekle_2021_03_17-DUB-deplete/mageck-depletion/24H_vs_24N")
df <- read.csv(file="24H_vs_24N.gene_summary.txt", sep="\t")


df$ensid = mapIds(org.Hs.eg.db,
                  keys=df$id, 
                  column="ENSEMBL",
                  keytype="SYMBOL",
                  multiVals="first")

gene.list <- df$ensid

result <- addList(david, gene.list,idType="ENSEMBL_GENE_ID",listName="gene.list", listType="Gene")

setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

termCluster <- getClusterReport(david, type="Term")

getClusterReportFile(david, type="Term", fileName="termClusterReport1.tab")
