#!/usr/bin/env Rscript
library(yaml)

args = commandArgs(trailingOnly=TRUE)






plot.bagel2 <- function(){
  df <- read.csv(file=args[1], sep="\t")
  save.path <- args[2]
  title <- args[3]
  
  df.top.ten <- df[1:10, ]
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
}


plot.mageck <- function(){
  library(ggplot2)
  library(ggrepel)
  suppressWarnings(suppressMessages(library("dplyr")))
  
  args = commandArgs(trailingOnly=TRUE)
  fdr <- args[1]
  
  df <- read.csv(file=args[2], sep="\t")
  df$fdr.cutoff <- fdr
  save.path <- args[3]
  
  #function to plot top 10 hits
  plot.hits <- function(input){
    
    if(input == "neg"){
      x <- "neg.rank"
      y <- "neg.score"
      z <- "neg.fdr"
      title <- "Drop out gene ranking"
      file <- "dropout.pdf"
    } else if(input == "pos") {
      x <- "pos.rank"
      y <- "pos.score"
      z <- "pos.fdr"
      title <- "Enrichment gene ranking"
      file <- "enrichment.pdf"
    }
    
    df$log.score <- -log10(df[[y]]) #for plotting
    
    #calculates score value for fdr cut off line in plot
    df$fdr.diff <- as.numeric(df$fdr.cutoff) - as.numeric(df[[z]])
    df$fdr.diff.abs <- abs(df$fdr.diff)
    fdr.min <- min(df$fdr.diff.abs)
    df.temp <- df[df$fdr.diff.abs == fdr.min,]
    fdr.cut.off.df <- df.temp[df.temp$log.score == max(df.temp$log.score), ]
    fdr.cut.off <- fdr.cut.off.df$log.score
    
    #determines top 10 hits for ggrepel labels
    df <- arrange(df, x)
    df.label <- df[df[[x]] %in% 1:10, ]
    
    options(ggrepel.max.overlaps = Inf)
    
    #generates plot
    p <- ggplot(df, aes_string(x = x, y = df$log.score)) +
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
      ylab("-log(MAGeCK score)") +
      guides(color = FALSE,
             shape =FALSE) +
      ggtitle(title) +
      geom_point(df, 
                 alpha = 0.6,
                 shape = 1,
                 size = 5,
                 colour = "black",
                 mapping = aes_string(x = x,
                                      y = df$log.score)) + 
      geom_hline(yintercept= fdr.cut.off, 
                 linetype="dashed", 
                 color = "red") +
      annotate("text", 
               x = nrow(df)*0.95, 
               y = fdr.cut.off*1.05, 
               label = paste0("FDR < ",fdr),
               size = 5,
               colour="red") 
    
    if(input == "neg"){
      p <- p + geom_label_repel(data = df.label, 
                                aes(x = `neg.rank`, y = `log.score`, label = `id`))
    } else if(input == "pos") {
      p <- p + geom_label_repel(data = df.label, 
                                aes(x = `pos.rank`, y = `log.score`, label = `id`))
    }
    
    ggsave(plot=p,
           file=paste0(save.path,file),
           width = 7,
           height = 5,
           useDingbats = FALSE)
  }
  
  #Plots top 10 hits for dropout and enrichment
  input <- c("neg","pos")
  lapply(input, plot.hits)
}