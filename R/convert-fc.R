library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)

fc.file <- args[1]
output.dir <- args[2]
out.file <- args[3]

df <- read.csv(file=fc.file, sep="\t")

#convert to CERES fc format (.gct)
df$GENE <- df$REAGENT_ID
names(df)[1:2] <- c("Name","Description")

write.table(df,
            file=paste0(output.dir,out.file,".ceres.foldchange"),
            sep="\t",
            row.names=FALSE)