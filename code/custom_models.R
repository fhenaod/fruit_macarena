library(tidyverse)
library(treeplyr)
library(bayou)

# data preparation ####
tree <- read.tree("data/tree_calib.nwk")
tree$tip.label <- tolower(tree$tip.label)
fruit_dat_raw <- read.csv("data/Base plantas Tinigua modificable_Agosto.csv", header = T )
fruit_dat_raw$ESPECIE <- tolower(fruit_dat_raw$ESPECIE)
fruit_dat_raw <- filter(fruit_dat_raw, !is.na(LARGO_FRUTO.cm.))

fruit_dat <- fruit_dat_raw %>% select(ESPECIE,
                                      LARGO_FRUTO.cm., Fr.width..cm.,
                                      Largo.hojas.total.cm., Largo.lamina.cm., 
                                      Area.total.cm2., Area.lamina.cm2.,
                                      Semillas.Fruto,
                                      Prom.Largo.semilla..mm., Prom.Ancho.sem..mm.,
                                      DBH.2017,
                                      Densidad..specific.gravity.,
                                      Tipo.de.hoja,
                                      Sistema_de_Dispersion,
                                      Color.Janson,
                                      FAMILIA) 

tree_data <- make.treedata(tree, fruit_dat)
tree_data <- tree_data %>%  mutate(fr_len = log(LARGO_FRUTO.cm.), fr_wd = log(Fr.width..cm.),
                                   lv_len = log(Largo.hojas.total.cm.), lam_len = log(Largo.lamina.cm.),
                                   ar_tot = log(Area.total.cm2.), ar_lam = log(Area.lamina.cm2.),
                                   sem_fr = log(Semillas.Fruto), 
                                   sem_len = log(Prom.Largo.semilla..mm.), sem_wd = log(Prom.Ancho.sem..mm.),
                                   dbh = log(DBH.2017),  
                                   dens = log(Densidad..specific.gravity.))

tree_data$dat %>% group_by(FAMILIA) %>% 
  summarize(Mean = mean(log(LARGO_FRUTO.cm.), na.rm = TRUE)) %>% arrange(desc(Mean))

tree_data <- filter(tree_data, !is.na(ar_tot)) %>% filter(!is.na(sem_fr)) %>% 
  filter(!is.na(dbh)) %>% filter(!is.na(sem_len))

cache <- bayou:::.prepare.ou.univariate(tree_data$phy, tree_data[["fr_len"]], 
                                        pred = select(tree_data$dat, ar_tot, sem_fr, dbh, sem_len))

sd2 <- sqrt(log(1+ (var(pull(tree_data$dat, fr_len), na.rm = T)) / (mean(pull(tree_data$dat, fr_len), na.rm = T))^2 ))

## priors ####
par.alpha <- list(scale = 1)
par.sig2 <- list(scale = 1)
par.beta_ar_tot <- list(mean = 0, sd = .1) 
par.beta_sem_fr <- list(mean = 0, sd = .1) 
par.beta_dbh <- list(mean = 0, sd = .1) 
par.beta_sem_len <- list(mean = 0, sd = .1) 

par.k <- list(lambda = (2*Ntip(tree_data$phy)-2)*(2.5/100), kmax = (2*Ntip(tree_data$phy)-2)*(5/100))
par.sb <- list(bmax = 1, prob = 1)
par.theta <- list(mean = mean(getVector(tree_data, fr_len)), sd = sd2)

# Tuning pars & model making ####
D.XXX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_ar_tot = 0.3, theta = 3,    slide = 1, missing.pred = 1)
D.1XX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_ar_tot = 0.3, theta = 0.15, slide = 1, missing.pred = 1)
D.XNX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_ar_tot = 0.3, theta = 1,    slide = 1, missing.pred = 1)

# Code explanation. r numbers given either: R - RJ; N - Fixed multiple; 1 - Fixed global; 0 - Absent
# 1:β0; 2: β_area_tot; 3: β_sem_fr ; 4: β_dbh; 5: β_sem_len
 
#####
mod <- "RR000"
prior.RR000 <- make.prior(tree_data$phy, plot.prior = FALSE, 
                          dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                       dbeta_ar_tot = "dnorm",
                                       #dbeta_sem_fr = "dnorm",
                                       #dbeta_dbh = "dnorm",
                                       #dbeta_sem_len = "dnorm",
                                       dsb = "dsb", dk = "cdpois", dtheta = "dnorm"), 
                          param = list(dalpha = par.alpha, dsig2 = par.sig2,
                                       dbeta_ar_tot = par.beta_ar_tot, sd = sd2,
                                       #dbeta_sem_fr = par.beta_sem_fr, 
                                       #dbeta_dbh = par.beta_dbh,
                                       #:dbeta_sem_len = par.beta_sem_len, 
                                       dk = par.k, dsb = par.sb, dtheta = par.theta)
)

model.RR000 <- makeBayouModel(fr_len ~ ar_tot, rjpars = c("theta", "ar_tot"), tree = tree_data$phy,  
                              dat = getVector(tree_data, fr_len),
                              pred = tree_data$dat, prior.RR000, D = D.XXX(2))

prior.RR000(model.RR000$startpar)
model.RR000$model$lik.fn(model.RR000$startpar, cache, cache$dat)$loglik

gens <- 10000
closeAllConnections()
mcmc.RR000 <- bayou.makeMCMC(tree_data$phy, getVector(tree_data, fr_len), pred = select(tree_data$dat, ar_tot), 
                             model = model.RR000$model, prior = prior.RR000, startpar = model.RR000$startpar,
                             new.dir = paste0(mod, "/"), outname = paste(mod, "run1", sep = "_"), plot.freq = NULL, 
                             ticker.freq = 2000, samp = 200, perform.checks = FALSE)

mcmc.RR000$run(gens)
chain.RR000 <- mcmc.RR000$load(saveRDS = T, file = paste0(mod, "/", "mcmc_", mod, ".rds"))
chain.RR000 <- set.burnin(chain.RR000, 0.3)
saveRDS(chain.RR000, file = paste0(mod, "/", "chain_", mod, ".rds"))
sum_c <- summary(chain.RR000)
saveRDS(sum_c, file = paste0(mod, "/", "sum_", mod, ".rds"))

shiftsum <- shiftSummaries(chain.RR000, mcmc.RR000)
saveRDS(shiftsum, file = paste0(mod, "/", "shiftsum_", mod, ".rds"))
pdf(paste0(mod, "/", "shiftsummaryplot.pdf"))
plotShiftSummaries((shiftsum))
dev.off()

require(foreach)
require(doParallel)
registerDoParallel(cores = 3)
Bk <- qbeta(seq(0, 1, length.out = 30), 0.3, 1)
ss.RR000 <- mcmc.RR000$steppingstone(gens, chain.RR000, Bk = Bk, burnin = 0.3, plot = F)
saveRDS(ss.RR000, paste0(mod, "/", "ss_", mod, ".rds"))
plot(ss.RR000)
print(ss.RR000$lnr)

#

fixed.k <- shift_sum$pars$k
fixed.sb <- shift_sum$pars$sb 
fixed.loc <- shift_sum$pars$loc

mod <- "NN000"
prior.NNN00 <- make.prior(tree_data$phy, plot.prior = FALSE, 
                          dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                       dbeta_ar_tot = "dnorm",
                                       #dbeta_sem_fr = "dnorm",
                                       #dbeta_dbh = "dnorm",
                                       #dbeta_sem_len = "dnorm",
                                       dsb = "dsb", dk = "cdpois", dtheta = "dnorm"), 
                          param = list(dalpha = par.alpha, dsig2 = par.sig2,
                                       dbeta_ar_tot = par.beta_ar_tot, sd = sd2,
                                       dbeta_sem_fr = par.beta_sem_fr, 
                                       #dbeta_dbh = par.beta_dbh,
                                       #dbeta_sem_len = par.beta_sem_len, 
                                       dk = par.k, dsb = par.sb, dtheta = par.theta)
)