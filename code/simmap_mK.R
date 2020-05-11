library(tidyverse)
library(treeplyr)
library(bayou)

# Data housekeeping ####
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
                                      Color.Janson) 

tree_data <- make.treedata(tree, fruit_dat)
tree_data <- tree_data %>%  mutate(fr_len = log(LARGO_FRUTO.cm.), fr_wd = log(Fr.width..cm.),
                                   lv_len = log(Largo.hojas.total.cm.), lam_len = log(Largo.lamina.cm.),
                                   ar_tot = log(Area.total.cm2.), ar_lam = log(Area.lamina.cm2.),
                                   sem_fr = log(Semillas.Fruto), 
                                   sem_len = log(Prom.Largo.semilla..mm.), sem_wd = log(Prom.Ancho.sem..mm.),
                                   dbh = log(DBH.2017),  
                                   dens = log(Densidad..specific.gravity.))

tree_data <- tree_data %>% filter(Sistema_de_Dispersion %in% c("Anemocorica", "Dehiscencia_explosiva", "Endozoocorica", 
                                                  "Epizoocorica", "Hidrocorica", "Inasistida", "Sinzoocorica"))

tree_data$dat <- droplevels(tree_data$dat)
summary(getVector(tree_data, Sistema_de_Dispersion))

# Simulate stochastic character maps ####

disp_simp_er <- make.simmap(tree_data$phy, getVector(tree_data, Sistema_de_Dispersion), model = "ER", nsim = 1000)
saveRDS(disp_simp_er, "output/disp_simp_er.rds")
disp_simp_ard <- make.simmap(tree_data$phy, getVector(tree_data, Sistema_de_Dispersion), model = "ARD", nsim = 1000)
saveRDS(disp_simp_ard, "ou:tput/disp_simp_ard.rds")
disp_simp_sym <- make.simmap(tree_data$phy, getVector(tree_data, Sistema_de_Dispersion), model = "SYM", nsim = 1000)
saveRDS(disp_simp_sym, "output/disp_simp_sym.rds")

disp_simp_er <- readRDS("output/disp_simp_er.rds")
disp_simp_ard <- readRDS("output/disp_simp_ard.rds")
disp_simp_sym <- readRDS("output/disp_simp_sym.rds")

disp_simp_er_sum <- summary(disp_simp_er, plot = F)
disp_simp_ard_sum <- summary(disp_simp_ard, plot = F)
disp_simp_sym_sum <- summary(disp_simp_sym, plot = F)

cols <- setNames(palette()[1:length(levels(getVector(tree_data, Sistema_de_Dispersion)))],levels(getVector(tree_data, Sistema_de_Dispersion)))

# Equal Rates
plot(disp_simp_er_sum, cols, type = "fan", fsize = .01, lwd = .5, cex = c(.4,.2))
add.simmap.legend(colors = cols, x = 0.95*par()$usr[1],y = 0.9*par()$usr[4],prompt = FALSE, fsize = 0.6)

disp_simp_er_dd <- density(disp_simp_er)
plot(disp_simp_er_dd)

# All-Rates-Different
plot(disp_simp_ard_sum, cols, type = "fan", fsize = .01, lwd = .5, cex = c(.4,.2))
add.simmap.legend(colors = cols, x = 0.95*par()$usr[1], y = 0.9*par()$usr[4], prompt = FALSE, fsize = 0.6)

disp_simp_ard_dd <- density(disp_simp_ard)
plot(disp_simp_ard_dd)

# Symmetrical model
plot(disp_simp_sym_sum, cols, type = "fan", fsize = .01, lwd = .5, cex = c(.4,.2))
add.simmap.legend(colors = cols, x = 0.95*par()$usr[1], y = 0.9*par()$usr[4], prompt = FALSE, fsize = 0.6)

disp_simp_sym_dd <- density(disp_simp_sym)
plot(disp_simp_sym_dd)

# Mk models ####
disp_er <- fitMk(tree_data$phy, getVector(tree_data, Sistema_de_Dispersion), model = "ER")
saveRDS(disp_er, "output/disp_er.rds")
disp_ard <- fitMk(tree_data$phy, getVector(tree_data, Sistema_de_Dispersion), model = "ARD")
saveRDS(disp_ard, "output/disp_ard.rds")
disp_sym <- fitMk(tree_data$phy, getVector(tree_data, Sistema_de_Dispersion), model = "SYM")
saveRDS(disp_sym, "output/disp_sym.rds")

disp_er <- readRDS("output/disp_er.rds")
disp_ard <- readRDS("output/disp_ard.rds") 
disp_sym <- readRDS("output/disp_sym.rds")

mod_sum <- data.frame(Mk_model = c("ARD","SYM","ER"),
           logLik = c(logLik(disp_ard), logLik(disp_sym), logLik(disp_er)),
           k = c(attr(AIC(disp_ard),"df"), attr(AIC(disp_sym),"df"), attr(AIC(disp_er),"df")),
           AIC = c(AIC(disp_ard), AIC(disp_sym), AIC(disp_er)),
           ΔAIC = round(qpcR::akaike.weights(c(AIC(disp_ard), AIC(disp_sym), AIC(disp_er)))$deltaAIC, 2),
           AICw = round(qpcR::akaike.weights(c(AIC(disp_ard), AIC(disp_sym), AIC(disp_er)))$weights, 2)
           )

plot(disp_er, show.zeros = F, main = "Equal rates")
plot(disp_ard, show.zeros = F, main = "All-rates different")
plot(disp_sym, show.zeros = F, main = "Symmetrical rates")
