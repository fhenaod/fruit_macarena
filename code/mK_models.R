library(tidyverse)
library(treeplyr)
library(bayou)

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

disp_simp_er <- make.simmap(tree_data$phy, getVector(tree_data, Sistema_de_Dispersion), model = "ER", nsim = 1000)
saveRDS(disp_simp_er, "output/disp_simp_er.rds")
disp_simp_ard <- make.simmap(tree_data$phy, getVector(tree_data, Sistema_de_Dispersion), model = "ARD", nsim = 1000)
saveRDS(disp_simp_ard, "ou:tput/disp_simp_ard.rds")
disp_simp_sym <- make.simmap(tree_data$phy, getVector(tree_data, Sistema_de_Dispersion), model = "SYM", nsim = 1000)
saveRDS(disp_simp_sym, "output/disp_simp_sym.rds")

plot(disp_simp_er, setNames(palette()[1:length(levels(getVector(tree_data, Sistema_de_Dispersion)))],
                            levels(getVector(tree_data, Sistema_de_Dispersi0n))), type = "fan", fsize = .1, lwd = .9, ftype = "i")

disp_er <- fitMk(tree_data$phy, getVector(tree_data, Sistema_de_Dispersion), model = "ER")
saveRDS(disp_er, "output/disp_er.rds")
disp_ard <- fitMk(tree_data$phy, getVector(tree_data, Sistema_de_Dispersion), model = "ARD")
saveRDS(disp_ard, "output/disp_ard.rds")
disp_sym <- fitMk(tree_data$phy, getVector(tree_data, Sistema_de_Dispersion), model = "SYM")
saveRDS(disp_sym, "output/disp_sym.rds")

rbind(disp_er$logLik, disp_ard$logLik, disp_sym$logLik)

plot(disp_er, show.zeros = F, main = "Equal rates")
plot(disp_ard, show.zeros = F, main = "Equal rates")
plot(disp_sym, show.zeros = F, main = "Equal rates")
