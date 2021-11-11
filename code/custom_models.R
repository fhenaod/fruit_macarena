library(tidyverse)
library(treeplyr)
library(bayou)
library(ggpubr)

# data load and clean ####
tree <- read.tree("data/tree_mac_s3.tre")
tree$edge.length[tree$edge.length<=0] <- 1e-6
tree <- extract.clade(tree, 1056) # extrae Spermatophyta
tree$tip.label <- tolower(tree$tip.label)

fruit_dat_raw <- read.csv("data/Base plantas Tinigua modificable_Agosto.csv", header = T)
fruit_dat_raw$ESPECIE <- tolower(fruit_dat_raw$ESPECIE)
fruit_dat_raw$ESPECIE <- gsub(pattern = " ", replacement = "_", x =fruit_dat_raw$ESPECIE)
#fruit_dat_raw <- filter(fruit_dat_raw, !is.na(LARGO_FRUTO.cm.))

fruit_dat <- fruit_dat_raw %>% select(ESPECIE,
                                      LARGO_FRUTO.cm., Fr.width..cm.,
                                      Largo.hojas.total.cm., Largo.lamina.cm., 
                                      Area.total.cm2., Area.lámina.cm2.,
                                      Semillas.Fruto,
                                      Prom.Largo.semilla..mm., Prom.Ancho.sem..mm.,
                                      DBH.2017, Densidad..specific.gravity.,
                                      Tipo.de.hoja,
                                      Sistema_de_Dispersion = Sistema_de_Dispersión,
                                      #Color.Janson,
                                      FAMILIA) 

tree_data <- make.treedata(tree, fruit_dat)
tree_data <- tree_data %>%  mutate(fr_len = log(LARGO_FRUTO.cm.), fr_wd = log(Fr.width..cm.),
                                   lv_len = log(Largo.hojas.total.cm.), lam_len = log(Largo.lamina.cm.),
                                   ar_tot = log(Area.total.cm2.), ar_lam = log(Area.lámina.cm2.),
                                   sem_fr = log(Semillas.Fruto+1), 
                                   sem_len = log(Prom.Largo.semilla..mm.+1), sem_wd = log(Prom.Ancho.sem..mm.+1),
                                   dbh = log(DBH.2017),  
                                   dens = log(Densidad..specific.gravity.))

#tree_data <- filter(tree_data, !is.na(fr_len)) 
   #%>% filter(!is.na(sem_fr)) %>%  filter(!is.na(dbh)) %>% 
  #filter(!is.na(sem_len))

tree_data$dat %>% group_by(FAMILIA) %>% 
  summarize(Mean = mean(log(LARGO_FRUTO.cm.), na.rm = TRUE)) %>% arrange(desc(Mean))

tree_data$dat %>% group_by(FAMILIA) %>% 
  summarize(Mean = mean(log(Prom.Largo.semilla..mm.+1), na.rm = TRUE)) %>% arrange(desc(Mean))

tree_data$dat %>% group_by(Sistema_de_Dispersion) %>% 
  count() %>% arrange(desc(n))

# recode dispersal systems
library(nlme)
tree_data <- tree_data %>% filter(Sistema_de_Dispersion != "")

tree_data$dat <- tree_data$dat %>% 
  mutate(disp_lato = 
           factor(ifelse(Sistema_de_Dispersion == 
                           "Endozoocórica", "endozoocorica", "no_endozoocorica")))
tree_data$dat %>% pull(disp_lato) %>% summary()

tree_data$dat <- tree_data$dat %>% 
  mutate(disp_3cat = 
           factor(ifelse(Sistema_de_Dispersion == "Endozoocórica", "Endozoochoric",
                         ifelse(Sistema_de_Dispersion == "Anemocórica", "Anemochoric",
                         "Non endozoochoric"))))

tree_data$dat %>% pull(disp_3cat) %>% summary()

# phylogenetic signal per trait #####
phy_sig_fr <- filter(tree_data, !is.na(fr_len)) %>% treedply(list("K" = phylosig(phy, getVector(., fr_len), "K"), 
                         "lambda" = phylosig(phy, getVector(., fr_len),"lambda")))

phy_sig_ln <- filter(tree_data, !is.na(lv_len)) %>% treedply(list("K" = phylosig(phy, getVector(., lv_len), "K"), 
                         "lambda" = phylosig(phy, getVector(., lv_len),"lambda")))

phy_sig_ar <- filter(tree_data, !is.na(ar_tot)) %>% treedply(list("K" = phylosig(phy, getVector(., ar_tot), "K"), 
                         "lambda" = phylosig(phy, getVector(., ar_tot),"lambda")))


phy_sig_se_l <- filter(tree_data, !is.na(sem_len)) %>% treedply(list("K" = phylosig(phy, getVector(., sem_len), "K"), 
                                       "lambda" = phylosig(phy, getVector(., sem_len),"lambda")))

phy_sig_se_w <- filter(tree_data, !is.na(sem_wd)) %>% treedply(list("K" = phylosig(phy, getVector(., sem_wd), "K"), 
                                       "lambda" = phylosig(phy, getVector(., sem_wd),"lambda")))

phy_sig_dens <- filter(tree_data, !is.na(dens)) %>% treedply(list("K" = phylosig(phy, getVector(., dens), "K"), 
                                       "lambda" = phylosig(phy, getVector(., dens),"lambda")))

# phylolm regression ####

library(phylolm)
#phyloANCOVA fruit length vs leave feats
rownames(tree_data$dat) <- tree_data$phy$tip.label

fit_fr <- phylolm(formula = fr_len ~ lv_len, 
                   data = tree_data$dat, phy = tree_data$phy, 
                   model = "lambda")

fit1 <- phylolm(fr_len ~ lv_len*disp_3cat, 
                  data = tree_data$dat, phy = tree_data$phy, 
                  model = "lambda")

fit1.1 <- phylolm(fr_len ~ lv_len+disp_3cat, 
                data = tree_data$dat, phy = tree_data$phy, 
                model = "lambda")

fit2 <- phylolm(fr_len ~ ar_tot*disp_3cat, 
                data = tree_data$dat, phy = tree_data$phy, 
                model = "lambda")

fit2.1 <- phylolm(fr_len ~ ar_tot+disp_3cat, 
                data = tree_data$dat, phy = tree_data$phy, 
                model = "lambda")

data.frame(
  predc = c("fit_fr", "fit1","fit1.1" #, "ar*disp", "ar+disp"
            ),
  logLik = c(fit_fr$logLik, fit1$logLik, fit1.1$logLik#, fit2$logLik, fit2.1$logLik
             ),
  lambda = c(fit_fr$optpar, fit1$optpar, fit1.1$optpar#, fit2$optpar, fit2.1$optpar
             ),
  #sigma = c(fr_size$sigma2, fit1$sigma2, fit1.1$sigma2 #, fit2$sigma2, fit2.1$sigma2
  #          ),
  r2 = c(fit_fr$r.squared, fit1$r.squared, fit1.1$r.squared #, fit2$r.squared, fit2.1$r.squared
         ),
  aicw(c(fit_fr$aic, fit1$aic, fit1.1$aic#, fit2$aic, fit2.1$aic
  ))
  #AICw = c(phytools::aic.w(c(fr_size$aic, fit1$aic, fit1.1$aic#, fit2$aic, fit2.1$aic
  #                           )))
) %>% arrange(desc(w)) %>% mutate(across(where(is.numeric), round, 3))

fit1.1 %>% summary()
fit_fr %>% summary()
lm(fr_len ~ lv_len, data = tree_data$dat) %>% summary()
lm(fr_len ~ lv_len+disp_3cat, data = tree_data$dat) %>% summary()

#ml_step <- phylostep(fr_len ~ lv_len + ar_tot + dens,
#                     data = tree_data$dat, phy = tree_data$phy, 
#                     model = "lambda", direction = "both")

#phyloANCOVA seed isom
tree_data$dat <- tree_data$dat %>% 
  mutate(sem_len_new = ifelse(sem_len > sem_wd, sem_len, sem_wd),
         sem_wd_new = ifelse(sem_wd < sem_len, sem_wd, sem_len)) 

seed_iso <- phylolm(sem_len_new ~ sem_wd_new, 
                    data = tree_data$dat, phy = tree_data$phy, 
                    model = "lambda")

seed_iso_ln <- phylolm(sem_len_new ~ sem_wd_new, 
                    data = tree_data$dat, phy = tree_data$phy, 
                    model = "lambda")

seed_iso_disp <- phylolm(sem_len_new ~ sem_wd_new + disp_3cat, 
                    data = tree_data$dat, phy = tree_data$phy, 
                    model = "lambda")

seed_isoxdisp <- phylolm(sem_len_new ~ sem_wd_new * disp_3cat, 
                         data = tree_data$dat, phy = tree_data$phy, 
                         model = "lambda")

seed_isoxdisp_sat <- phylolm(exp(sem_len_new) ~ (sem_wd_new) * disp_3cat, 
                         data = tree_data$dat, phy = tree_data$phy, 
                         model = "lambda")

data.frame(
  predc = c("isom", "isom+disp", "seed_iso x disp"),
  logLik = c(seed_iso$logLik, seed_iso_disp$logLik, seed_isoxdisp$logLik),
  lambda = c(seed_iso$optpar, seed_iso_disp$optpar, seed_isoxdisp$optpar),
  sigma =  c(seed_iso$sigma2, seed_iso_disp$sigma2, seed_isoxdisp$sigma2),
  r2 = c(seed_iso$r.squared, seed_iso_disp$r.squared, seed_isoxdisp$r.squared),
  AIC = c(seed_iso$aic, seed_iso_disp$aic, seed_isoxdisp$aic), 
  AICw = c(phytools::aic.w(c(seed_iso$aic, seed_iso_disp$aic, seed_isoxdisp$aic)))
) %>% arrange(desc(AICw)) %>% mutate(across(where(is.numeric), round, 3))

summary(seed_isoxdisp)

# ou shifts in a trait
OUshifts(getVector(tree_data, fr_len), tree_data$phy,
         method = "mbic", nmax = 2)

# pgls fruit length vs leave length or area ####
tr_dt_fr_lv <- filter(tree_data, !is.na(fr_len)) %>% filter(!is.na(lv_len))
tr_dt_fr_ar <- filter(tree_data, !is.na(fr_len)) %>% filter(!is.na(ar_tot))
                                                            
bm_ln <- gls(fr_len ~ lv_len, correlation = corBrownian(phy = tr_dt_fr_lv$phy),
               data = tr_dt_fr_lv$dat, method = "ML")

bm_ar <- gls(fr_len ~ ar_tot, correlation = corBrownian(phy = tr_dt_fr_ar$phy), 
               data = tr_dt_fr_ar$dat, method = "ML")

bm_ln_disp <- gls(fr_len~lv_len*disp_lato, correlation = corBrownian(phy = tr_dt_fr_lv$phy), 
               data = tr_dt_fr_lv$dat, method = "ML")
bm_ar_disp <- gls(fr_len~ar_tot*disp_lato, correlation = corBrownian(phy = tr_dt_fr_ar$phy), 
               data = tr_dt_fr_ar$dat, method = "ML")

data.frame(
length = coef(bm_ln_disp),
area = coef(bm_ar_disp)
) %>% round(3)

tr_dt_fr_lv_tmp <- tr_dt_fr_lv$phy
tr_dt_fr_ar_tmp <- tr_dt_fr_ar$phy

tr_dt_fr_lv_tmp$edge.length <- tr_dt_fr_lv_tmp$edge.length*100
tr_dt_fr_ar_tmp$edge.length <- tr_dt_fr_ar_tmp$edge.length*100

ou_ln <- gls(fr_len~lv_len, correlation = corMartins(1, phy = tr_dt_fr_lv_tmp), 
             data = tr_dt_fr_lv$dat, method = "ML")

ou_ar <- gls(fr_len~ar_tot, correlation = corMartins(1, phy = tr_dt_fr_ar_tmp), 
             data = tr_dt_fr_ar$dat, method = "ML")


ou_ln_disp <- gls(fr_len~lv_len*disp_lato, correlation = corMartins(1, phy = tr_dt_fr_lv_tmp), 
                  data = tr_dt_fr_lv$dat, method = "ML")

ou_ar_disp <- gls(fr_len~ar_tot*disp_lato, correlation = corMartins(1, phy = tr_dt_fr_ar_tmp), 
                  data = tr_dt_fr_ar$dat, method = "ML")

s_bm_ln <- summary(bm_ln_disp)
s_ou_ln <- summary(ou_ln_disp)
s_bm_ar <- summary(bm_ar_disp)
s_ou_ar <- summary(ou_ar_disp)

data.frame(
model = c("BM", "OU", "BM", "OU"#, "BM", "OU"
          ),
predc = c("lv_len", "lv_len", "lv_len_disp", "lv_len_disp" #, "lv_ar", "lv_ar"
          ),
logLik = c(bm_ln$logLik ,ou_ln$logLik, s_bm_ln$logLik, s_ou_ln$logLik #, s_bm_ar$logLik, s_ou_ar$logLik
           ),
AIC = c(AIC(bm_ln),AIC(ou_ln), s_bm_ln$AIC, s_ou_ln$AIC #, s_bm_ar$AIC, s_ou_ar$AIC
        ), 
AICw = c(aic.w(c(AIC(bm_ln),AIC(ou_ln), s_bm_ln$AIC, s_ou_ln$AIC#, s_bm_ar$AIC, s_ou_ar$AIC
                 )))
) %>% arrange(desc(AICw)) %>% mutate(across(where(is.numeric), round, 3))

anova(s_bm_ln)
coef(s_bm_ln) %>% data.frame() %>% round(3)

# pgls seed isometry
tr_dt_seed <- filter(tree_data, !is.na(sem_len)) %>% filter(!is.na(sem_wd))
tr_dt_seed <- filter(tr_dt_seed, !is.infinite(sem_len)) %>% filter(!is.infinite(sem_wd))
tr_dt_seed$phy$edge.length <- tr_dt_seed$phy$edge.length*100

bm_seed_disp <- gls(sem_len~sem_wd*disp_lato, correlation = corBrownian(phy = tr_dt_seed$phy), 
                  data = tr_dt_seed$dat, method = "ML")
ou_seed_disp <- gls(sem_len~sem_wd*disp_lato, correlation = corMartins(1, phy = tr_dt_seed$phy), 
                  data = tr_dt_seed$dat, method = "ML")

s_bm_se <- bm_seed_disp %>% summary()
s_ou_se <- ou_seed_disp %>% summary()

data.frame(
  model = c("BM", "OU"),
  #predc = c("lv_len", "lv_len", "lv_ar", "lv_ar"),
  logLik = c(s_bm_se$logLik, s_ou_se$logLik),
  AIC = c(s_bm_se$AIC, s_ou_se$AIC), 
  AICw = c(aic.w(c(s_bm_se$AIC, s_ou_se$AIC)))
) %>% arrange(desc(AICw))

anova(s_ou_se)
coef(s_ou_se) %>% data.frame() %>% round(3)

tree_data$dat %>% ggplot(aes(x = sem_len, y = sem_wd, 
                             shape = disp_lato, color = disp_lato)) +
  geom_point() + geom_smooth(method = "lm", se = T, fullrange = F) + 
  theme_classic() +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = .5) +
  theme() + 
  labs(x = "Ln Seed length (mm)", y = "Ln Seed width (mm)")

###### bayou ######
cache <- bayou:::.prepare.ou.univariate(tree_data$phy, tree_data[["fr_len"]], 
                                        pred = select(tree_data$dat, ar_tot, lv_len, sem_fr, dbh, sem_len))

sd2 <- sqrt(log(1+ (var(pull(tree_data$dat, fr_len), na.rm = T)) / (mean(pull(tree_data$dat, fr_len), na.rm = T))^2 ))
sd2

## priors ####
par.alpha <- list(scale = 1)
par.sig2 <- list(scale = 1)
par.beta_ar_tot <- list(mean = 0, sd = .1)
par.beta_lv_len <- list(mean = 0, sd = .1)
par.beta_sem_fr <- list(mean = 0, sd = .1) 
par.beta_dbh <- list(mean = 0, sd = .1) 
par.beta_sem_len <- list(mean = 0, sd = .1) 

par.k <- list(lambda = (2*Ntip(tree_data$phy)-2)*(2.5/100), kmax = (2*Ntip(tree_data$phy)-2)*(5/100))
par.sb <- list(bmax = 1, prob = 1)
par.theta <- list(mean = mean(getVector(tree_data, fr_len)), sd = sd2)

## RR0 ####
# Tuning pars & model making
D.XXX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_ar_tot = 0.3, theta = 3,    slide = 1, missing.pred = 1)
D.1XX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_ar_tot = 0.3, theta = 0.15, slide = 1, missing.pred = 1)
D.XNX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_ar_tot = 0.3, theta = 1,    slide = 1, missing.pred = 1)

# Code explanation. 5 numbers given either: R - RJ; N - Fixed multiple; 1 - Fixed global; 0 - Absent
# 1:β0; 2: β_area_tot; 3: β_lv_len

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
                                       #dbeta_sem_len = par.beta_sem_len, 
                                       dk = par.k, dsb = par.sb, dtheta = par.theta)
)
prior.RR000

model.RR000 <- makeBayouModel(fr_len ~ ar_tot, rjpars = c("theta", "ar_tot"), tree = tree_data$phy,  
                              dat = getVector(tree_data, fr_len),
                              pred = tree_data$dat, prior.RR000, D = D.XXX(2))
model.RR000

prior.RR000(model.RR000$startpar)
model.RR000$model$lik.fn(model.RR000$startpar, cache, cache$dat)$loglik

# run model 
gens <- 5000000
closeAllConnections()
mcmc.RR000 <- bayou.makeMCMC(tree_data$phy, getVector(tree_data, fr_len), pred = select(tree_data$dat, ar_tot), 
                             model = model.RR000$model, prior = prior.RR000, startpar = model.RR000$startpar,
                             file.dir = paste0(mod, "/"), outname = paste(mod, "run1", sep = "_"), plot.freq = NULL, 
                             ticker.freq = 2000, samp = 200, perform.checks = FALSE)

saveRDS(mcmc.RR000, "RR000/mcmc.conf.RR000.rds")
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
registerDoParallel(cores = 15)
Bk <- qbeta(seq(0, 1, length.out = 30), 0.3, 1)
ss.RR000 <- mcmc.RR000$steppingstone(gens, chain.RR000, Bk = Bk, burnin = 0.3, plot = F)
saveRDS(ss.RR000, paste0(mod, "/", "ss_", mod, ".rds"))
#plot(ss.RR00)
print(ss.RR000$lnr)

# R0R ####
tree_data <- filter(tree_data, !is.na(lv_len))

# Code explanation. 5 numbers given either: R - RJ; N - Fixed multiple; 1 - Fixed global; 0 - Absent
# 1:β0; 2: β_area_tot; 3: β_lv_len
D.XXX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_lv_len = 0.3, theta = 3,    slide = 1, missing.pred = 1)
D.1XX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_lv_len = 0.3, theta = 0.15, slide = 1, missing.pred = 1)
D.XNX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_lv_len = 0.3, theta = 1,    slide = 1, missing.pred = 1)

mod <- "R0R"
prior.R0R <- make.prior(tree_data$phy, plot.prior = FALSE,  
                        dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy",
                                     dbeta_lv_len = "dnorm",
                                     #dbeta_sem_fr = "dnorm",
                                     #dbeta_dbh = "dnorm",
                                     #dbeta_sem_len = "dnorm",
                                     dsb = "dsb", dk = "cdpois", dtheta = "dnorm"),
                        param = list(dalpha = par.alpha, dsig2 = par.sig2,
                                     dbeta_lv_len = par.beta_lv_len, sd = sd2,
                                     #dbeta_sem_fr = par.beta_sem_fr, 
                                     #dbeta_dbh = par.beta_dbh,
                                     #dbeta_sem_len = par.beta_sem_len, 
                                     dk = par.k, dsb = par.sb, dtheta = par.theta)
)
prior.R0R

model.R0R <- makeBayouModel(fr_len ~ lv_len, rjpars = c("theta", "lv_len"), tree = tree_data$phy,
                              dat = getVector(tree_data, fr_len),
                              pred = tree_data$dat, prior.R0R, D = D.XXX(2))

prior.R0R(model.R0R$startpar)
model.R0R$model$lik.fn(model.R0R$startpar, cache, cache$dat)$loglik

# run model 
gens <- 5000000
closeAllConnections()
mcmc.R0R <- bayou.makeMCMC(tree_data$phy, getVector(tree_data, fr_len), pred = select(tree_data$dat, lv_len),
                           model = model.R0R$model, prior = prior.R0R, startpar = model.R0R$startpar,
                           file.dir = paste0(mod, "/"), outname = paste(mod, "run1", sep = "_"), plot.freq = NULL,
                           ticker.freq = 2000, samp = 200, perform.checks = FALSE)


saveRDS(mcmc.R0R, "R0R/mcmc.conf.R0R.rds")
mcmc.R0R$run(gens) 
chain.R0R <- mcmc.R0R$load(saveRDS = T, file = paste0(mod, "/", "mcmc_", mod, ".rds"))
chain.R0R <- set.burnin(chain.R0R, 0.3)
saveRDS(chain.R0R, file = paste0(mod, "/", "chain_", mod, ".rds"))
sum_c <- summary(chain.R0R)
saveRDS(sum_c, file = paste0(mod, "/", "sum_", mod, ".rds"))

shiftsum <- shiftSummaries(chain.R0R, mcmc.R0R)
saveRDS(shiftsum, file = paste0(mod, "/", "shiftsum_", mod, ".rds"))
pdf(paste0(mod, "/", "shiftsummaryplot.pdf"))
par(mfrow = c(2,2))
plotShiftSummaries(shiftsum, single.plot = F, label.pts = T)
dev.off()

require(foreach)
require(doParallel)
registerDoParallel(cores = 15)
Bk <- qbeta(seq(0, 1, length.out = 30), 0.3, 1)
ss.R0R <- mcmc.R0R$steppingstone(gens, chain.R0R, Bk = Bk, burnin = 0.3, plot = F)
saveRDS(ss.R0R, paste0(mod, "/", "ss_", mod, ".rds"))
#plot(ss.RR00)
print(ss.R0R$lnr)
