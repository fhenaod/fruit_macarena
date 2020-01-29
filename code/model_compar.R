library(dplyr)
library(treeplyr)
library(bayou)
library(ggplot2)

path <- ("custom_models/")
d <- list.files(path, all.files = T, include.dirs = T, recursive = T)
files <- d[grep("ss_", d)]
models <- sapply(strsplit(files, "/"), "[[", 1)

mar_Lik <- c()
for(i in 1:length(models)){
  mar_Lik[i] <-readRDS(paste0(path, files[i]))$lnr
}

mod_comp <- data.frame(models, mar_Lik)
mod_comp$BF <-round(abs(2*(mod_comp$mar_Lik[which(max(mod_comp$mar_Lik)==mar_Lik)]-mod_comp$mar_Lik)), 2)
mod_comp[which(max(mod_comp$BF)==mod_comp$BF),] # Best model
