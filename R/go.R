options(show.error.locations = TRUE)
options(error=traceback)

library(ggplot2)
library(viridis)
library(enrichR)
#library(GOplot)
library(stringr)
#library(yaml)

args <- commandArgs(trailingOnly=TRUE)

file <- args[1]
save.path <- args[2]
save.path <- file.path(save.path,"enrichR")
dir.create(save.path)
analysis <- args[3]

go.mageck <- function(file,save.path,analysis){
  df <- read.csv(file, sep="\t")
  
  #get genes below FDR cut off
  df.fdr <- df[df$neg.fdr < 0.25, ]
  
  #enrichR databases
  dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
  
  enriched <- enrichr(df.fdr$id, dbs)
  
  for (i in dbs){
    df <- enriched[[i]]
    df.plot <- df[1:10,]
    df.plot$log.pvalue <- -log10(df.plot$Adjusted.P.value)
    
    #calculate genes:group size ratio
    df.plot$count <- str_split_fixed(df.plot$Overlap,"/", n=2)[,1]
    df.plot$count <- as.numeric(df.plot$count)
    df.plot$group.size <- str_split_fixed(df.plot$Overlap,"/", n=2)[,2]
    df.plot$group.size <- as.numeric(df.plot$group.size)
    df.plot$ratio <- df.plot$count / df.plot$group.size
    
    #plot results
    p <- ggplot(df.plot, aes(x=reorder(`Term`,`log.pvalue`),y=`log.pvalue`, fill=`ratio`)) +
      geom_bar(stat='identity') + 
      coord_flip() + 
      labs(x='',y='-log(adj. p value)') +
      theme_bw(base_size = 16)+
      theme(legend.key.size= unit(1.5, "cm"))+
      scale_fill_viridis(name = "genes:group size") +
      ggtitle(paste0("Enriched in ",dbs[3]))
    
    ggsave(filename=file.path(save.path,paste0("enrichR-",i,".pdf")),
           p,
           width = 20,
           height= 5)
  }
  
  
}

go.bagel2 <- function(){}


if(analysis == "mageck"){
  go.mageck(file,save.path,analysis)
} else if (analysis == "bagel2"){
  go.bagel2()
}

traceback()