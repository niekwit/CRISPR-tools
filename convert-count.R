#!/usr/bin/env Rscript
library(stringr)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)

setwd(args[1])
fasta.file <- args[2]
out.folder <- args[3]

df <- read.csv(file="count/counts-aggregated.tsv",sep="\t")
df$sgRNA <- paste0(df$gene,"_",df$sgRNA)
df <- df[ , ! names(df) %in% "gene"]

fasta <- read.csv(file=fasta.file, header=FALSE)

df.guides <- as.data.frame(fasta[grep(">", fasta$V1), ])
names(df.guides)[1] <- "sgRNA"

df.guides$sgRNA <- gsub(">","",df.guides$sgRNA)

df.guides$SEQUENCE <- fasta[grep(">", fasta$V1,invert=TRUE), ]
df.guides$GENE <- word(df.guides$sgRNA,1,sep="\\_")
df.guides$sgRNA <- sub("_*\\.","",df.guides$sgRNA)

df.output <- merge(df.guides,df,by="sgRNA")
df.output <- df.output[ , ! names(df.output) %in% "sgRNA"]

write.table(df.output,
            file=paste0(out.folder,"/counts-aggregated-",out.folder,".txt"),
            sep="\t",
            row.names=FALSE)