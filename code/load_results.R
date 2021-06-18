library(dplyr)
library(treeplyr)
library(bayou)
library(ggplot2)

# RR0 ####
chain.RR0 <- readRDS("custom_models/RR0/RR0/chain_RR000.rds")
sum_RR0 <- readRDS("custom_models/RR0/RR0/sum_RR000.rds")
shift_sum_RR0 <- readRDS("custom_models/RR0/RR0/shiftsum_RR000.rds")

sum_RR0$statistics
sum_RR0$branch.posteriors

shift_sum_RR0$descendents
sapply(shift_sum_RR0$descendents, length)

shift_sum_RR0$regressions

plotSimmap.mcmc(chain.RR0, burnin = 0.3, pp.cutoff = 0.3, cex = .01)
plotBranchHeatMap(tree_data$phy, chain.RR0, "theta", burnin = 0.3, pal = cm.colors, cex = .1)
plotBranchHeatMap(tree_data$phy, chain.RR0, "beta_ar_tot", burnin = 0.3, pal = cm.colors, cex = .1)

phenogram.density(tree_data$phy, getVector(tree_data, fr_len), burnin = 0.3, chain.RR0, pp.cutoff = 0.3, 
                  xlab = "Time (Myr)", ylab = "Fruit lenght", spread.labels=TRUE)
phenogram.density(tree_data$phy, getVector(tree_data, ar_tot), burnin = 0.3, chain.RR0, pp.cutoff = 0.3, 
                  xlab = "Time (Myr)", ylab = "Total area", spread.labels=TRUE)

# Phylogenetic half-life
round(data.frame(mean = log(2)/sum_RR0$statistics["alpha","Mean"], 
                 hpdL = log(2)/sum_RR0$statistics["alpha","HPD95Lower"], 
                 hpdU = log(2)/sum_RR0$statistics["alpha","HPD95Upper"]), 3)

# R0R ####
chain.R0R <- readRDS("custom_models/R0R/R0R/mcmc_R0R.rds")
mcmc.R0R <- readRDS("custom_models/R0R/R0R/mcmc.conf.R0R.rds")
sum_R0R <- readRDS("custom_models/R0R/R0R/sum_R0R.rds")
shift_sum_R0R <- readRDS("custom_models/R0R/R0R/shiftsum_R0R.rds")
ss_R0R <- readRDS("custom_models/R0R/")

sum_R0R$statistics
sum_R0R$branch.posteriors

shift_sum_R0R$descendents
sapply(shift_sum_R0R$descendents, length)

shift_sum_R0R$regressions

plotSimmap.mcmc(chain.R0R, burnin = 0.3, pp.cutoff = 0.3, cex = .01)

# Phylogenetic half-life
round(data.frame(mean = log(2)/sum_R0R$statistics["alpha","Mean"], 
                 hpdL = log(2)/sum_R0R$statistics["alpha","HPD95Lower"], 
                 hpdU = log(2)/sum_R0R$statistics["alpha","HPD95Upper"]), 3)

