#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)
library(GO.db)

args = commandArgs(trailingOnly=TRUE)

fdr <- args[1]
df <- read.csv(file=args[2]) #check seperator for gene_summary.txt
species <- args[3]

if (species == "human") {
	library(org.Hs.eg.db)
} else if (species == "mouse"){
	library(org.Mm.eg.db)
}

df$neg.log.negscore <- -log10(df$neg.score)

p <- ggplot(df, aes(x=`id`, y=`neg.log.negscore`)) +
        theme_bw() +
	theme(axis.text=element_text(size=16),
              axis.title=element_text(size=16),
              plot.title = element_text(hjust = 0.5),
              axis.line = element_line(colour = "black"),
              panel.border = element_blank(),
              panel.background = element_blank()) +
        xlab("Gene") +
        ylab("-log(MAGeCK score)") +
        guides(color = FALSE,
               shape =FALSE)

ggsave()
