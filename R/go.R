library(ggplot2)
library(viridis)
library(enrichR)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)

file <- args[1]
save.path <- args[2]
save.path <- file.path(save.path,"enrichR")
dir.create(save.path)
analysis <- args[3]
fdr.cutoff <- args[4]

if(analysis == "mageck"){
  go.mageck <- function(fdr.column,
                        rank,
                        title,
                        input=file,
                        savepath=save.path,
                        fdr=fdr.cutoff){
    
    dir.create(savepath, 
               showWarnings = FALSE)
    
    df <- read.csv(input, sep="\t")
    
    #enrichR databases
    dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
    
    
    #get genes below FDR cut off
    df.fdr <- df[df[fdr.column] < fdr, ]
    
    if(nrow(df.fdr) != 0){
      enriched <- enrichr(df.fdr$id, dbs)
      for (j in dbs){
        df <- enriched[[j]]
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
        
        ggsave(filename=file.path(savepath,paste0("enrichR-",title,j,".pdf")),
               p,
               width = 20,
               height= 5)
      } 
    } else {
      next
    }
  }
  
  neg <- c("neg.fdr","neg.rank","dropout-")
  pos <- c("pos.fdr","pos.rank","enrichment-")
  
  go.mageck(neg[1],neg[2],neg[3])
  go.mageck(pos[1],pos[2],pos[3])
} else if (analysis == "bagel2"){
  
}
