setwd("/home/niek/CRISPR-tools")
df <- read.csv(file="counts-aggregated.tsv", sep= '\t')

df.output <- df
for(x in 3:ncol(df)) {
  sum.column <- sum(df[x])
  df.output[x] <-ceiling(df[x]/sum.column*1E8)
}

write.csv(df.output, 
          file = "counts-aggregated-normalised.tsv", 
          row.names = FALSE)
