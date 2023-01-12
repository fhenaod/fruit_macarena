# pca ####
library(nlme)
library(vegan); library(factoextra); library(FactoMineR)

pca_trs <- 
  tree_data$dat %>% 
  dplyr::select(LARGO_FRUTO.cm., Largo.hojas.total.cm., Area.lámina.cm2.,
                Semillas.Fruto, Prom.Largo.semilla..mm., Prom.Ancho.sem..mm.,
                DBH.2017, Densidad..specific.gravity., dispersal = disp_3cat) %>% 
  drop_na() %>% PCA(., graph = F, quali.sup = 9)

# log-transformed
pca_trs <- 
  tree_data$dat %>% 
  dplyr::select(fruit_length = fr_len, leaf_length = lv_len, 
                laminar_area = ar_lam, seeds_fruit = sem_fr, 
                seed_length = sem_len, seed_width = sem_wd, 
                dbh, wood_density = dens, dispersal = disp_3cat) %>% 
  drop_na() %>% PCA(., graph = F, quali.sup = 9)

fviz_pca_biplot(pca_trs, col.var = "cos2", col.ind = "contrib", 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = T) + theme_classic()
plotellipses(pca_trs, 9)

fviz_eig(pca_trs, addlabels = T)
pca_vars <- get_pca_var(pca_trs)

# pca plot most important contributing vars
fviz_pca_var(pca_trs, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             #alpha.var = "contrib",
             repel = T) + theme_classic()
plotellipses(pca_trs, 9)

fviz_contrib(pca_trs, choice = "var", axes = 1, top = 10) # Vars contributions to PC1
fviz_contrib(pca_trs, choice = "var", axes = 2, top = 10) # Vars contributions to PC2

fviz_contrib(pca_trs, choice = "var", axes = 1:2, top = 10) # Vars contributions to PC1&2

# phylo PCA ####
library(phytools)

dat2phy_pca <- 
  tree_data %>% 
  dplyr::select(fr_len, fr_wd, lv_len, lam_len, 
                #ar_tot, ar_lam, 
                sem_fr, sem_len, sem_wd, 
                dbh, #dens, 
                dispersal = disp_3cat) 

dat2phy_pca %>% dplyr::select(-contains("NA"))
drop_na() %>% skimr::skim()

rownames(dat2phy_pca$dat) <- dat2phy_pca$phy$tip.label

phylo_pca <- phyl.pca(dat2phy_pca$phy, dat2phy_pca$dat, method = "lambda")

# LDA ####
library(caret)
library(MASS)

cat2filt <- tree_data$dat %>% 
  dplyr::select(LARGO_FRUTO.cm., Largo.hojas.total.cm., 
                Area.lámina.cm2., Prom.Largo.semilla..mm., FAMILIA) %>% 
  drop_na() %>% group_by(FAMILIA) %>% tally %>% 
  filter(n > 3) %>% dplyr::select(FAMILIA)

daf2lda <- tree_data$dat %>% 
  dplyr::select(LARGO_FRUTO.cm., Largo.hojas.total.cm., Area.lámina.cm2.,
                Prom.Largo.semilla..mm., FAMILIA) %>% 
  drop_na() %>% 
  filter(FAMILIA %in% c(cat2filt)$FAMILIA) %>% droplevels()

training.samples <- daf2lda %>% pull(FAMILIA) %>% 
  createDataPartition(p = .8, list = F)

train.data <- daf2lda[training.samples, ]
test.data <- daf2lda[-training.samples, ]

preproc.param <- train.data %>% 
  preProcess(method = c("center", "scale"))

train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)

# 
m_lda <- lda(FAMILIA~., data = train.transformed)
predictions <- m_lda %>% predict(test.transformed) # predictions
mean(predictions$class==test.transformed$FAMILIA) # model accuracy

lda.data <- cbind(train.transformed, predict(m_lda)$x)
ggplot(lda.data, aes(LD1, LD2)) + geom_point(aes(color = FAMILIA)) +
  theme_classic()

ggord::ggord(m_lda, train.transformed$FAMILIA) + 
  theme_classic() + theme(legend.position = "none")

ldahist(data = predictions$x[,1], 
        g = train.transformed$FAMILIA)

# LMM ####
library(lme4)

gen <- 
  fruit_dat %>% 
  mutate(genus = str_split(fruit_dat$ESPECIE, "_", n = 2, simplify = F) %>% 
           sapply("[", 1), genus = as.factor(genus)) %>% 
  lmer(log(Largo.hojas.total.cm.) ~ (1|genus), data = .)

fam <- 
  fruit_dat %>% mutate(FAMILIA = as.factor(FAMILIA)) %>% 
  lmer(log(Largo.hojas.total.cm.) ~ (1|FAMILIA), data = .)

gen %>% lattice::qqmath(id = 0.05)
fam %>% lattice::qqmath(id = 0.05)

VarCorr(gen, comp = c("Variance","Std.Dev."))
VarCorr(fam, comp = c("Variance","Std.Dev."))

# bayou ######
library(bayou)
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

## R0R ####
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
