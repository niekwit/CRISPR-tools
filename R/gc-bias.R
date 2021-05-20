#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
#setwd(args[1])
#fasta <- args[2]

setwd("/mnt/sdb1/analyses/CRISPR-screens/guine/count")
fasta <- "/home/niek/Documents/references/fasta/Mouse/Mouse-SLC-Mito-2OGDD/mouse-2-OG-DD-SLC-mitoP-library.fasta"
df <- read.csv(file="counts-aggregated.tsv", 
               sep="\t")

df <- subset(x=df,
             select=c("sgRNA","gene","pre","post"))

#convert fasta to csv format
df.fasta <- read.csv(fasta,
                     header=FALSE)

sgRNA.index <- grep(">", df.fasta$V1, invert=FALSE)
sequence.index <- grep(">", df.fasta$V1, invert=TRUE)

df.csv <- as.data.frame(x=df.fasta[sgRNA.index, ])
df.csv$sequence <- df.fasta[sequence.index, ]
names(df.csv) <- c("sgRNA","sequence")
df.csv$sgRNA <- substring(df.csv$sgRNA, 2) #removes < from sgRNA names
df.csv$sgRNA <- gsub("^.*?_","",df.csv$sgRNA) #removes gene_ from gene_sg_gene_number
df.csv <- df.csv[order(df.csv$sgRNA), ]

#add sequences to main df
df.merge <- merge(x=df,
                  y=df.csv,
                  by="sgRNA")

#determine GC bias of guides with the lowest count
df.merge <- df.merge[order(df.merge$pre), ]
length <- nrow(df.merge)
fraction <- round(0.2 * length,0)
df.low <- df.merge[1:fraction,] #selects guides with lowest counts (bottom 5%)


#calculate GC content of bottom 5% count sgRNAs
df.low$gc.content <- NA
for (i in 1:fraction) {
  seq <- df.low$sequence[i]
  len.seq <- str_count(seq, pattern="")
  G <- str_count(seq, pattern="G")
  C <- str_count(seq, pattern="C")
  gc.content <- (G + C)/len.seq*100
  df.low$gc.content[i] <- gc.content
  }
low.gc.content <-mean(df.low$gc.content)

#calculate average GC content of all sgRNAs in library
df.merge$gc.content <- NA
for (i in 1:nrow(df.merge)) {
  seq <- df.merge$sequence[i]
  len.seq <- str_count(seq, pattern="")
  G <- str_count(seq, pattern="G")
  C <- str_count(seq, pattern="C")
  gc.content <- (G + C)/len.seq*100
  df.merge$gc.content[i] <- gc.content
}
library.gc.content <- mean(df.merge$gc.content)

#plot histogram
ggplot(df.low, aes(x=gc.content)) +
  geom_histogram(binwidth=3) +
  geom_vline(aes(xintercept=library.gc.content),
             colour="navy",
             linetype="dashed",
             size=1)
  
