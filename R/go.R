library(ggplot2)
library(viridis)
library(enrichR)
library(GOplot)
library(stringr)
library(yaml)

args <- commandArgs(trailingOnly=TRUE)

file <- args[1]
save.path <- args[2]
analysis <- args[3]

go.mageck <- function(file,save.path,analysis){
  #setwd("/run/user/1000/gvfs/smb-share:server=jcbc-store01.jcbc.private.cam.ac.uk,share=citiid_home/nw416/My Documents/Data/HDAC3 project/Results/Screen E data")
  
  df <- read.csv(file)
  
  #get genes below FDR cut off
  df.fdr <- df[df$neg.fdr < 0.25, ]
  
  #enrichR databases
  dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
  
  enriched <- enrichr(d$id, dbs)
  
  
  df <- enriched[[dbs[3]]]
  
  df.plot <- df[1:10,]
  df.plot$log.pvalue <- -log10(df.plot$Adjusted.P.value)
  
  #calculate genes:group size ratio
  df.plot$count <- str_split_fixed(df.plot$Overlap,"/", n=2)[,1]
  df.plot$count <- as.numeric(df.plot$count)
  df.plot$group.size <- str_split_fixed(df.plot$Overlap,"/", n=2)[,2]
  df.plot$group.size <- as.numeric(df.plot$group.size)
  df.plot$ratio <- df.plot$count / df.plot$group.size
  
  #plot enrichPathway results
  p <- ggplot(df.plot, aes(x=reorder(`Term`,`log.pvalue`),y=`log.pvalue`, fill=`ratio`)) +
    geom_bar(stat='identity') + 
    coord_flip() + 
    labs(x='',y='-log(adj. p value)') +
    theme_bw(base_size = 16)+
    theme(legend.key.size= unit(1.5, "cm"))+
    scale_fill_viridis(name = "genes:group size") +
    ggtitle(paste0("Enriched in ",dbs[3]))
  
  ggsave(filename="enriched-pathways.pdf",
         p,
         width = 20,
         height= 5)
}

go.bagel2 <- function(){}


if(analysis == "mageck"){
  go.mageck(file,save.path,analysis)
} else if (analysis == "bagel2"){
  go.bagel2()
}



