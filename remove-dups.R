args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

df <- read.csv("counts-aggregated.tsv", sep="\t")
df <- unique.data.frame(df)

write.table(df.unique, 
            file="counts-aggregated.tsv", 
            sep="\t", 
            quote=FALSE, 
            row.names = FALSE)