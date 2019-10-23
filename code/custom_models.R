library(dplyr)
library(treeplyr)
library(bayou)

## priors ####
par.alpha <- list(scale = 1)
par.sig2 <- list(scale = 1)
par.beta_ar_tot <- list(mean = 0, sd = .1) 

par.k <- list(lambda = (2*Ntip(tree_data$phy)-2)*(2.5/100), kmax = (2*Ntip(tree_data$phy)-2)*(5/100))
par.sb <- list(bmax = 1, prob = 1)
par.theta <- list(mean = mean(getVector(tree_data, fr_len)), sd = sd2)

# Tuning pars & model making ####
D.XXX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_ar_tot=0.2, theta = 3,    slide = 1, missing.pred = 1)
D.1XX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_ar_tot=0.2, theta = 0.15, slide = 1, missing.pred = 1)
D.XNX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_ar_tot=0.2, theta = 1,    slide = 1, missing.pred = 1)

tree_data <- filter(tree_data, !is.na(ar_tot))

cache <- bayou:::.prepare.ou.univariate(tree_data$phy, tree_data[["fr_len"]], pred=select(tree_data$dat, ar_tot))
prior.RR000 <- make.prior(tree_data$phy, plot.prior = FALSE, 
                          dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                       dbeta_ar_tot = "dnorm",
                                       dsb = "dsb", dk = "cdpois", dtheta = "dnorm"), 
                          param = list(dalpha = par.alpha, dsig2 = par.sig2,
                                       dbeta_ar_tot = par.beta_ar_tot, sd = sd2, # regresion coeff. "dbeta_predictor"
                                       dk = par.k, dsb = par.sb, dtheta = par.theta)
)
model.RR000 <- makeBayouModel(fr_len ~ ar_tot, rjpars = c("theta", "ar_tot"), tree = tree_data$phy,  
                              dat = getVector(tree_data, fr_len),
                              pred = tree_data$dat, prior.RR000, D = D.XXX(2))
prior.RR000(model.RR000$startpar)
model.RR000$model$lik.fn(model.RR000$startpar, cache, cache$dat)$loglik

mod<-"RR000"
gens<-10000
closeAllConnections()
mymcmc<-bayou.makeMCMC(tree_data$phy, getVector(tree_data, fr_len), pred = select(tree_data$dat, ar_tot), 
                       model = model.RR000$model, prior = prior.RR000, startpar = model.RR000$startpar,
                       new.dir = TRUE, outname = paste(mod, "run1", sep = "_"), plot.freq = NULL, 
                       ticker.freq = 2000, samp = 200, perform.checks = FALSE)

saveRDS(mymcmc, file = "modelOU_rj/mcmcOU_rj.rds")
mymcmc$run(gens)
chain <- mymcmc$load()
chain_burn <- set.burnin(chain, 0.3)
sum_c <- summary(chain_burn)
sum_c$branch.posteriors$branch.loc <- seq(1:length(sum_c$branch.posteriors$pp))
saveRDS(sum_c, file = "sum_ou_rj.rds")

shiftsum <- shiftSummaries(chain, mymcmc)
pdf("./shiftsummaryplot.pdf")
plotShiftSummaries((shiftsum))
dev.off()
