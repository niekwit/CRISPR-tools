#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  stop("Only one argument must be supplied (path to data folder)", call.=FALSE)
} else if (length(args)==1) {
  setwd(args[1])
}

df <- read.csv(file="counts-aggregated.tsv", sep= '\t')

df.output <- df
for(x in 3:ncol(df)) {
  sum.column <- sum(df[x],na.rm=TRUE)
  df.output[x] <-ceiling(df[x]/sum.column*1E8)
}

write.csv(df.output, 
          file = "counts-aggregated-normalised.csv", 
          row.names = FALSE)
