#!/usr/bin/env Rscript

#check packages
cran.packages <- c("ggplot2","yaml","ggrepel","viridis","dplyr","enrichR")
installed.packages <- installed.packages()[,1]

for (i in cran.packages){
  if(!i %in% installed.packages){
    cat("R: ",i,"package missing, installing now\n")
    install.packages(i,repos='http://cran.us.r-project.org')
  }
}
