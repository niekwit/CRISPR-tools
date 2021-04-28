#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("viridis")))

args = commandArgs(trailingOnly=TRUE)

df <- read.csv(file=args[1], sep="\t")
save.path <- args[2]
title <- args[3]

df.top.ten <- df[1:10, ]
#df.hits <- df[df$BF>5, ]
df <- arrange(df, df$Gene)

#plot hits
p <- ggplot(df, aes(x=`Gene`,y=`BF`)) +
  theme_bw() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("Genes") +
  ylab("log2 Bayes Factor") +
  geom_point(df, mapping=aes(size=`Precision`, fill=`Recall`)
             ,shape=21) +
  ggtitle(title) +
  scale_fill_viridis(discrete=FALSE,
                     guide = guide_colorbar(frame.colour = "black", 
                                                          ticks.colour = "black"), 
                     limits=c(0,1))+
  geom_hline(yintercept= 5, 
             linetype="dashed", 
             color = "red") +
  geom_label_repel(data = df.top.ten, 
                   aes(x = `Gene`, y = `BF`, label = `Gene`))

file.plot <- paste0("/bagel2-hits",title, ".pdf")

ggsave(plot=p,
     file=paste0(save.path,file.plot),
     width = 7,
     height = 5,
     useDingbats = FALSE)

df <- arrange(df, df$Recall)

#plot Precision-Recall plot
pp <- ggplot(df, aes(x=`Recall`,y=`Precision`)) +
  theme_bw() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank()) +
  xlab("Recall") +
  ylab("Precision (1-FDR)") +
  geom_line() +
  ggtitle(paste0("Precision-Recall plot ",title))

file.plot <- paste0("/precision-recall-",title, ".pdf")

ggsave(plot=pp,
       file=paste0(save.path,file.plot),
       width = 7,
       height = 5,
       useDingbats = FALSE)
