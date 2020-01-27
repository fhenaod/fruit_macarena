library(dplyr)
library(treeplyr)
library(bayou)
library(ggplot2)

d <- list.files("custom_models", all.files = T, include.dirs = T, recursive = T)
files <- d[grep("ss_", d)]
models <- sapply(strsplit(files, "/"), "[[", 1)
path <- ("custom_models/")
mar_Lik <- c()
for(i in 1:length(models)){
  mar_Lik[i] <-readRDS(paste0(path, files[i]))$lnr
}
sum_res <- data.frame(models, mar_Lik)
sum_res$BF <-round(abs(2*(-415.4994-sum_res$mar_Lik)), 2)
