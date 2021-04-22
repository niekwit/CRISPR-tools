#!/usr/bin/env Rscript
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

setwd(paste0(args[1],"/count"))
fasta.file <- args[2]

df <- read.csv(file="counts-aggregated.tsv",sep="\t")
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
            file="counts-aggregated-bagel.txt",
            sep="\t",
            row.names=FALSE)
