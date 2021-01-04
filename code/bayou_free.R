library(tidyverse)
library(treeplyr)
library(bayou)

tree <- read.tree("data/tree_calib.nwk")
tree <- read.tree("data/tree_mac_s3.tre")

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

contMap(tree_data$phy, getVector(tree_data, fr_len), fsize = c(.001,1), lwd = 2, 
        outline = F, sig = 2, type = "fan", leg.txt = "Ln fruit length (cm)")

contMap(filter(tree_data, !is.na(ar_tot))$phy, getVector(filter(tree_data, !is.na(ar_tot)), ar_tot), fsize = c(.001,1), lwd = 2, 
        outline = F, sig = 2, type = "fan", leg.txt = "Ln total area (cm2)")

treedply(tree_data, list("K" = phytools::phylosig(phy, getVector(tree_data, fr_len), "K"),
                         "lambda" = phytools::phylosig(phy, getVector(tree_data, fr_len), "lambda")))

getVector(tree_data, ar_tot)[!is.na()]

# Free OU
sd2<-sqrt(log(1+ (var(pull(tree_data$dat, fr_len), na.rm = T)) / (mean(pull(tree_data$dat, fr_len), na.rm = T))^2 ))

priorOU<-make.prior(tree_data$phy, 
                    dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", dk = "cdpois", dtheta = "dnorm"),
                    param = list(dalpha = list(scale = 0.1), dsig2 = list(scale = 0.1),
                                 dsb = list(bmax = 1, prob = 1), 
                                 dk = list(lambda = (2*Ntip(tree_data$phy)-2)*(2.5/100), kmax = (2*Ntip(tree_data$phy)-2)*(5/100)),
                                 dtheta = list(mean = mean(getVector(tree_data, fr_len)), 
                                               sd = sd2))
)

mcmcOU<-bayou.makeMCMC(tree_data$phy, getVector(tree_data, fr_len),
                       prior = priorOU, new.dir = "modelOU/", outname = "modelOU_r001", plot.freq = NULL,
                       ticker.freq = 2000, samp = 200) 

saveRDS(mcmcOU, file = "modelOU/mcmcOU.rds")
gens<-10000
mcmcOU$run(gens)
setwd("~/")
mcmcOU<-readRDS("modelOU/mcmcOU.rds")
chainOU<-mcmcOU$load()
#chainOU<-mcmcOU$load(saveRDS = T, file = 'modelOU/mcmcOU_chain.rds')
chainOU<-readRDS("mac_bayou/modelOU/mcmcOU_chain.rds")
chainOU<-set.burnin(chainOU, 0.3)
sum_c<-summary(chainOU)
sum_c$branch.posteriors$branch.loc<-seq(1:length(sum_c$branch.posteriors$pp))
saveRDS(sum_c, file = "sum_free_ou.rds")

sum_c$branch.posteriors %>% filter(pp >= 0.3) %>% arrange(desc(pp)) # which/where shifts
sum_c$branch.posteriors %>% filter(pp >= 0.3) %>%  tally() %>% 
  pull()/length(tree_data$phy$edge.length)*100 %>% as.numeric() # %_branch shifts

plotSimmap.mcmc(chainOU, burnin = 0.3, pp.cutoff = 0.3, cex = .01, no.margin = T)
plotBranchHeatMap(tree_data$phy, chainOU, "theta", burnin = 0.3, pal = cm.colors, cex = .1, type = "fan")
phenogram.density(tree_data$phy, getVector(tree_data, fr_len), burnin = 0.3, chainOU, pp.cutoff = 0.3, 
                  xlab = "Time (Myr)", ylab = "Phenotype", spread.labels = TRUE)

# RJ
setwd("Documents/UBC/Projects/fruit_macarena/ou_rj/")
mcmc_rj <- readRDS("model_ou_rj/mcmcOU_rj.rds")
chain_rj <- mcmc_rj$load()
chain_rj <- set.burnin(chain_rj, 0.3)
sum_rj <- readRDS("sum_ou_rj.rds")

shift_sum <- shiftSummaries(chain_rj, mcmc_rj, pp.cutoff = .3)
shift_sum$descendents

plotSimmap.mcmc(chain_rj, burnin = 0.3, pp.cutoff = 0.3, cex = .01, no.margin = T)
plotBranchHeatMap(tree_data$phy, chain_rj, "theta", burnin = 0.3, pal = cm.colors, cex = .1, type = "fan")
phenogram.density(tree_data$phy, getVector(tree_data, fr_len), burnin = 0.3, chain_rj, pp.cutoff = 0.3, 
                  xlab = "Time (Myr)", ylab = "Phenotype", spread.labels = TRUE)

# RJ V.PhyloTreeS3
tree_data$phy <- read.tree("RR000/data/tree_mac_s3.tre")
tree_data$phy$edge.length[tree$edge.length<=0] <- 1e-6
tree_data$phy <- extract.clade(tree, 1056) # extrae Spermatophyta

mcmc_rj <- readRDS("RR000/RR000/mcmc_RR000.rds")
chain_rj <- readRDS("RR000/RR000/chain_RR000.rds")
sum_rj <- readRDS("RR000/RR000/sum_RR000.rds")

shift_sum <- readRDS("RR000/RR000/shiftsum_RR000.rds")
shift_sum$tree
shift_sum$descendents
shift_sum$descendents[sapply(shift_sum$descendents, length)>3]
shift_sum$cladesummaries
shift_sum$regressions[sapply(shift_sum$descendents, length)>3,] %>% data.frame() %>% 
  mutate(rl_theta = exp(theta)) %>% round(3)

# phylo-half life
as.table(log(2)/sum_rj$statistics['alpha',c('Mean', 'HPD95Lower', 'HPD95Upper')]) %>% round(2)

plotSimmap.mcmc(chain_rj, burnin = 0.3, pp.cutoff = 0.3, cex = .01, no.margin = T)
plotBranchHeatMap(shift_sum$tree, chain_rj, "theta", burnin = 0, pal = cm.colors, cex = .1, type = "fan")
phenogram.density(shift_sum$tree, getVector(tree_data, fr_len), burnin = 0.3, chain_rj, pp.cutoff = 0.3, 
                  xlab = "Time (Myr)", ylab = "Phenotype", spread.labels = TRUE)

pdf(paste0("shiftsummaryplot.pdf"))
par(mfrow = c(2,2))
plotShiftSummaries(shift_sum, single.plot = F, label.pts = T) 
dev.off()

