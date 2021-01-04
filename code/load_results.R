library(dplyr)
library(treeplyr)
library(bayou)
library(ggplot2)

# RJ
chain.RR000 <- readRDS("custom_models/RR000/RR000/chain_RR000.rds")
sum_RR000 <- readRDS("custom_models/RR000/RR000/sum_RR000.rds")
shift_sum <- readRDS("custom_models/RR000/RR000/shiftsum_RR000.rds")
ss_RR000 <- readRDS("custom_models/RR000/RR000/ss_RR000.rds")

sum_RR000$statistics
sum_RR000$branch.posteriors

shift_sum$descendents
sapply(shift_sum$descendents, length)

shift_sum$regressions

plotSimmap.mcmc(chain.RR000, burnin = 0.3, pp.cutoff = 0.3, cex = .01)
plotBranchHeatMap(tree_data$phy, chain.RR000, "theta", burnin = 0.3, pal = cm.colors, cex = .1)
plotBranchHeatMap(tree_data$phy, chain.RR000, "beta_ar_tot", burnin = 0.3, pal = cm.colors, cex = .1)

phenogram.density(tree_data$phy, getVector(tree_data, fr_len), burnin = 0.3, chain.RR000, pp.cutoff = 0.3, 
                  xlab = "Time (Myr)", ylab = "Fruit lenght", spread.labels=TRUE)
phenogram.density(tree_data$phy, getVector(tree_data, ar_tot), burnin = 0.3, chain.RR000, pp.cutoff = 0.3, 
                  xlab = "Time (Myr)", ylab = "Total area", spread.labels=TRUE)

# Phylogenetic half-life
round(data.frame(mean = log(2)/sum_RR000$statistics["alpha","Mean"], 
                 hpdL = log(2)/sum_RR000$statistics["alpha","HPD95Lower"], 
                 hpdU = log(2)/sum_RR000$statistics["alpha","HPD95Upper"]), 3)
