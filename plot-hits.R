#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)
library(GO.db)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

fdr <- args[1]
df <- read.csv(file=args[2], sep="\t")
species <- args[3]
save.path <- args[4]

if (species == "human") {
	library(org.Hs.eg.db)
} else if (species == "mouse"){
	library(org.Mm.eg.db)
}

df$neg.log.score <- -log10(df$neg.score) #for plotting

#calculates score value for fdr cut off line in plot
df$fdr.cutoff <- fdr
df$fdr.diff <- as.numeric(df$fdr.cutoff) - as.numeric(df$neg.fdr)
df$fdr.diff.abs <- abs(df$fdr.diff)
fdr.min <- min(df$fdr.diff.abs)
df.temp <- df[df$fdr.diff.abs == fdr.min,]
fdr.cut.off.df <- df.temp[df.temp$neg.log.score == max(df.temp$neg.log.score), ]
fdr.cut.off <- fdr.cut.off.df$neg.log.score

#determines top 10 hits for ggrepel labels
df.label <- df[df$neg.rank %in% 1:10, ]

#parameters for plotting
plot.titles <- c("Dropout gene ranking","Enrichment gene ranking")
plot.df.neg <- df
plot.df.pos <- df[order(df$pos.rank), ]
plot.dfs <- as.list( c(plot.df.neg, plot.df.pos))

p <- ggplot(df, aes(x=`neg.rank`, y=`neg.log.score`)) +
        theme_bw() +
	theme(axis.text=element_text(size=16),
        axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.x = element_blank()) +
        xlab("Genes") +
        ylab("-log(MAGeCK score)") +
        guides(color = FALSE,
       		shape =FALSE) +
      	ggtitle("Dropout gene ranking") +
        geom_point(df, 
		        alpha = 0.6,
		        shape = 1,
		        size = 5,
		        colour = "black",
		        mapping = aes(x = `neg.rank`,
			                    y = `neg.log.score`)) + 
  geom_hline(yintercept= fdr.cut.off, 
             linetype="dashed", 
             color = "red") +
  annotate("text", 
        x = nrow(df)*0.95, 
        y = fdr.cut.off*1.05, 
        label = paste0("FDR < ",fdr),
        size = 5,
        colour="red")  +  
  geom_label_repel(size=4,
                   aes(x = `neg.rank`,
                       y = `neg.log.score`,
                       label = id), 
                   data = df.label,
                   nudge_x = nrow(df)*0.1,
                   nudge_y = max(df$neg.log.score)*0.1)

ggsave(plot=p,
       file=paste0(save.path,"drop-out.pdf"),
       width = 7,
       height = 5,
       useDingbats = FALSE)