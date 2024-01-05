library(tidyverse)
library(treeplyr)

# load and prepare data  ####
tree_s3 <- read.tree("data/tree_mac_s3.tre")
#tree <- read.tree("data/Genus_phylogeny.tre")
#tree <- read.tree("data/neves_species_phylo.tre")
tree <- tree_s3

fruit_dat_raw <- read.csv("data/Base plantas Tinigua 2023.csv", header = T )
fruit_dat_raw$ESPECIE <- (fruit_dat_raw$ESPECIE)

fruit_dat_raw$Genus <- (sapply(strsplit(as.character(fruit_dat_raw$ESPECIE), split = "_", fixed = T), "[", 1))

#fruit_dat_raw <- filter(fruit_dat_raw, !is.na(LARGO_FRUTO.cm.))
#fruit_dat_raw <- filter(fruit_dat_raw, !is.na(Sistema_de_Dispersion))

fruit_dat <- fruit_dat_raw %>% dplyr::select(Genus, ESPECIE,
                                      LARGO_FRUTO.cm., Fr.width..cm.,
                                      Largo.hojas.total.cm., Largo.lamina.cm., 
                                      Area.total.cm2., Area.lamina.cm2.,
                                      Prom.Largo.semilla..mm., Prom.Ancho.sem..mm.,
                                      DBH.2017, Habito, 
                                      Densidad..specific.gravity.,
                                      Tipo.de.hoja, TIPODEINFL, 
                                      TIPODEFRUT_General, TIPODEFRUT_Especifico, FORMAFRUT, 
                                      Sistema_de_Dispersion,
                                      Semillas.Fruto) 

tree_data <- make.treedata(tree, fruit_dat)
tree_data <- tree_data %>%  mutate(fr_len = log(LARGO_FRUTO.cm.), fr_wd = log(Fr.width..cm.),
                                   lv_len = log(Largo.hojas.total.cm.), lam_len = log(Largo.lamina.cm.),
                                   ar_tot = log(Area.total.cm2.), ar_lam = log(Area.lamina.cm2.),
                                   sem_len = log(Prom.Largo.semilla..mm.), sem_wd = log(Prom.Ancho.sem..mm.),
                                   dbh = log(DBH.2017),  
                                   dens = log(Densidad..specific.gravity.),
                                   sem_fr = log(Semillas.Fruto)
                                   )

# recode dispersal systems
tree_data$dat$Sistema_de_Dispersion <- recode(tree_data$dat$Sistema_de_Dispersion, 
                                              "Inasistida/Hidrocorica/Anemocorica" = "Anemocorica", 
                                              "Hidrocorica/Dehiscencia_explosiva" = "Dehiscencia_explosiva", 
                                              "Inasistida/Anemocorica" = "Anemocorica", 
                                              "Inasistida/Hidrocorica" = "Hidrocorica",
                                              "myrmecocorica" = "Myrmecocorica")

tree_data$dat$Sistema_de_Dispersion <- recode(tree_data$dat$Sistema_de_Dispersion, 
                                              'Endozoocorica' = "Endozoochory","Anemocorica" = "Anemochory", 
                                              'Inasistida' = "Unassisted", 'Dehiscencia_explosiva' = "Explosive dehiscence", 
                                              'Hidrocorica' = "Hydrochory", 'Sinzoocorica' = "Synzoochory", 
                                              "Epizoocorica" = "Epizoochory", "Myrmecocorica" = "Myrmecochory") 

# dispersal three-categories
tree_data$dat <- tree_data$dat %>% 
  mutate(disp_3cat = factor(ifelse(Sistema_de_Dispersion == "Endozoochory", "Endozoochorous", 
                                   ifelse(Sistema_de_Dispersion == "Anemochory", "Anemochorous", "Non-endozoochorous"))))

# recode fruit types
tree_data$dat$TIPODEFRUT_General <- recode(tree_data$dat$TIPODEFRUT_General, 
                                              "Baya_en_amento" = "Baya", "Baya_protegida" = "Baya", 
                                              "Capsula_indehiscente" = "Capsula", 
                                              "Legumbre_indehiscente" = "Legumbre")
tree_data$dat <- droplevels(tree_data$dat)

# Mk models ####
library(phytools)
dir.create("mk_models/output/")
n_cores <- 10

## fruit type ####
tree_data_fr <- tree_data %>% filter(TIPODEFRUT_General %in% 
                                         c("Capsula", "Baya", "Drupa", "Legumbre", "Samara", 
                                           "Amento", "Grano_de_pasto", "Aquenio", "Sicono", "Foliculo", "Lomento"))
runan <- T
if(runan){
  fr_type_er <- fitMk.parallel(tree_data_fr$phy, getVector(tree_data_fr, TIPODEFRUT_General), model = "ER", ncores = n_cores)
  saveRDS(fr_type_er, "mk_models/output/fr_type_er.rds")
  fr_type_ard <- fitMk.parallel(tree_data_fr$phy, getVector(tree_data_fr, TIPODEFRUT_General), model = "ARD", ncores = n_cores)
  saveRDS(fr_type_ard, "mk_models/output/fr_type_ard.rds")
  fr_type_sym <- fitMk.parallel(tree_data_fr$phy, getVector(tree_data_fr, TIPODEFRUT_General), model = "SYM", ncores = n_cores)
  saveRDS(fr_type_sym, "mk_models/output/fr_type_sym.rds") 
} else {
  fr_type_er  <- readRDS("mk_models/output/fr_type_er.rds")
  fr_type_ard <- readRDS("mk_models/output/fr_type_ard.rds") 
  fr_type_sym <- readRDS("mk_models/output/fr_type_sym.rds")
}

mod_sum <- data.frame(Mk_model = c("ARD","SYM","ER"),
                      logLik = c(logLik(fr_type_ard), logLik(fr_type_sym), logLik(fr_type_er)),
                      k = c(length(fr_type_ard$rates), length(fr_type_sym$rates), length(fr_type_er$rates)),
                      AIC = c(AIC(fr_type_ard), AIC(fr_type_sym), AIC(fr_type_er)),
                      deltaAIC = round(qpcR::akaike.weights(c(AIC(fr_type_ard), AIC(fr_type_sym), AIC(fr_type_er)))$deltaAIC, 2),
                      AICw = round(qpcR::akaike.weights(c(AIC(fr_type_ard), AIC(fr_type_sym), AIC(fr_type_er)))$weights, 3)
)
mod_sum %>% arrange(desc(AICw))

## dispersal system ####
tree_data_disp <- tree_data %>% filter(Sistema_de_Dispersion %in% 
                       c("Endozoochory", "Anemochory", "Unassisted", "Hydrochory", 
                         "Explosive dehiscence", "Synzoochory", "Epizoochory", "Myrmecochory"))
runan <- T
if(runan){
  disp_er <- fitMk.parallel(tree_data_disp$phy, getVector(tree_data_disp, Sistema_de_Dispersion), model = "ER", ncores = n_cores)
  saveRDS(disp_er, "mk_models/output/disp_er.rds")
  disp_ard <- fitMk.parallel(tree_data_disp$phy, getVector(tree_data_disp, Sistema_de_Dispersion), model = "ARD", ncores = n_cores)
  saveRDS(disp_ard, "mk_models/output/disp_ard.rds")
  disp_sym <- fitMk.parallel(tree_data_disp$phy, getVector(tree_data_disp, Sistema_de_Dispersion), model = "SYM", ncores = n_cores)
  saveRDS(disp_sym, "mk_models/output/disp_sym.rds") 
  } else {
  disp_er  <- readRDS("mk_models/output/disp_er.rds")
  disp_ard <- readRDS("mk_models/output/disp_ard.rds") 
  disp_sym <- readRDS("mk_models/output/disp_sym.rds")
}

mod_sum <- data.frame(Mk_model = c("ARD","SYM","ER"),
                      logLik = c(logLik(disp_ard), logLik(disp_sym), logLik(disp_er)),
                      k = c(length(disp_ard$rates), length(disp_sym$rates), length(disp_er$rates)),
                      AIC = c(AIC(disp_ard), AIC(disp_sym), AIC(disp_er)),
                      ΔAIC = round(qpcR::akaike.weights(c(AIC(disp_ard), AIC(disp_sym), AIC(disp_er)))$deltaAIC, 2),
                      AICw = round(qpcR::akaike.weights(c(AIC(disp_ard), AIC(disp_sym), AIC(disp_er)))$weights, 3)
)
mod_sum

# transition rates estimation and plots ####
library(circlize)

## fruit type ####
qm_ard <- fr_type_sym %>% phytools::as.Qmatrix() %>% as.table()
dimnames(qm_ard) = list(source = colnames(qm_ard), 
                        sink = rownames(qm_ard))

png("fruit_type_qm.png", width = 20, height = 20, units = "cm", res = 300)
cols <- hcl.colors(length(dimnames(qm_ard)$sink), "Geyser")
circos.clear()
circos.par(start.degree = 90, gap.degree = 5)
chordDiagram(x = qm_ard, #grid.col = cols, 
             grid.border = "black", 
             transparency = 0.5, order = rownames(qm_ard), 
             directional = T, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.1),
             annotationTrack = "grid", annotationTrackHeight = c(0.1, 0.1),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.15, link.arr.lty = 3, link.arr.col = "#252525", 
             link.largest.ontop = F)

# Add labels to chord diagram
circos.trackPlotRegion(track.index = 1, bg.border = NA,
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         sector.index = get.cell.meta.data("sector.index")
                         # Text direction
                         theta = circlize(mean(xlim), 1)[1, 1] %% 360
                         dd = ifelse(theta < 180 || theta > 360, "bending.inside", "bending.outside")
                         circos.text(x = mean(xlim), y = 0.8, labels = sector.index, facing = dd,
                                     niceFacing = TRUE, cex = 1, font = 4)})

# Add title
mtext("Transitions: Fruit type", outer = FALSE, cex = 2, font = 2, line = -1)
dev.off()

plot(fr_type_er, show.zeros = F, main = "Equal rates")
plot(fr_type_ard, show.zeros = F, main = "All-rates different")
plot(fr_type_sym, show.zeros = F, main = "Symmetrical rates")

## dispersal system #### 
qm_ard <- disp_ard %>% phytools::as.Qmatrix() %>% as.table()
dimnames(qm_ard) = list(source = colnames(qm_ard), 
                        sink = rownames(qm_ard))

png("disp_syst_qm.png", width = 20, height = 20, units = "cm", res = 300)
cols <- c("#FF410D", "#6EE2FF", "#F7C530", "#95CC5E", "#D0DFE6")
circos.clear()
circos.par(start.degree = 90, gap.degree = 6)
chordDiagram(x = qm_ard, grid.col = cols, grid.border = "black", 
             transparency = 0.5, order = rownames(qm_ard), 
             directional = T, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.1),
             annotationTrack = "grid", annotationTrackHeight = c(0.1, 0.1),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.15, link.arr.lty = 3, link.arr.col = "#252525", 
             link.largest.ontop = F)

# Add labels to chord diagram
circos.trackPlotRegion(track.index = 1, bg.border = NA,
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         sector.index = get.cell.meta.data("sector.index")
                         # Text direction
                         theta = circlize(mean(xlim), 1)[1, 1] %% 360
                         dd = ifelse(theta < 180 || theta > 360, "bending.inside", "bending.outside")
                         circos.text(x = mean(xlim), y = 0.8, labels = sector.index, facing = dd,
                                     niceFacing = TRUE, cex = 1, font = 4)})

# Add title
mtext("Transitions: Dispersal Systems ", outer = FALSE, cex = 2, font = 2, line = -1)
dev.off()

plot(disp_er, show.zeros = F, main = "Equal rates")
plot(disp_ard, show.zeros = F, main = "All-rates different")
plot(disp_sym, show.zeros = F, main = "Symmetrical rates")

# Simulate stochastic character maps ####
## fruit type ####
runan <- T
if(runan){
  fr_type_simp_er <- make.simmap(tree_data_fr$phy, getVector(tree_data_fr, TIPODEFRUT_General), model = "ER", nsim = 1000)
  saveRDS(fr_type_simp_er, "mk_models/output/fr_type_simp_er.rds")
  fr_type_simp_ard <- make.simmap(tree_data_fr$phy, getVector(tree_data_fr, TIPODEFRUT_General), model = "ARD", nsim = 1000)
  saveRDS(fr_type_simp_ard, "mk_models/output/fr_type_simp_ard.rds")
  fr_type_simp_sym <- make.simmap(tree_data_fr$phy, getVector(tree_data_fr, TIPODEFRUT_General), model = "SYM", nsim = 1000)
  saveRDS(fr_type_simp_sym, "mk_models/output/fr_type_simp_sym.rds")
} else {
  fr_type_simp_er  <- readRDS("mk_models/output/fr_type_simp_er.rds")
  fr_type_simp_ard <- readRDS("mk_models/output/fr_type_simp_ard.rds")
  fr_type_simp_sym <- readRDS("mk_models/output/fr_type_simp_sym.rds")
}

# simmaps summaries
runan <- T
if(runan){
  fr_type_simp_er_sum  <- summary(fr_type_simp_er, plot = T)
  saveRDS("mk_models/output/fr_type_er_sum.rds")
  fr_type_simp_ard_sum <- summary(fr_type_simp_ard, plot = F)
  saveRDS("mk_models/output/fr_type_ard_sum.rds")
  fr_type_simp_sym_sum <- summary(fr_type_simp_sym, plot = F)
  saveRDS("mk_models/output/fr_type_sym_sum.rds")
} else{ 
  fr_type_simp_er_sum  <- readRDS("mk_models/output/fr_type_er_sum.rds")
  fr_type_simp_ard_sum <- readRDS("mk_models/output/fr_type_ard_sum.rds")
  fr_type_simp_sym_sum <- readRDS("mk_models/output/fr_type_sym_sum.rds")
}

# paint it on the tree
#cols <- setNames(palette()[1:length(levels(getVector(tree_data, Sistema_de_Dispersion)))], 
#                 levels(getVector(tree_data, Sistema_de_Dispersion)))
cols <- setNames(palette()[1:dim(fr_type_simp_er_sum$ace)[2]],
                 colnames(fr_type_er_sum$ace))
# Equal Rates
plot(fr_type_simp_er_sum, cols, type = "fan", fsize = .01, lwd = .5, cex = c(.4,.2))
add.simmap.legend(colors = cols, x = 0.95*par()$usr[1],y = 0.9*par()$usr[4], prompt = FALSE, fsize = 0.6)

fr_type_simp_er_dd <- density(fr_type_simp_er)
plot(fr_type_simp_er_dd)

# All-Rates-Different
plot(fr_type_simp_ard_sum, cols, type = "fan", fsize = .01, lwd = .5, cex = c(.4,.2))
add.simmap.legend(colors = cols, x = 0.95*par()$usr[1], y = 0.9*par()$usr[4], prompt = FALSE, fsize = 0.6)

fr_type_simp_ard_dd <- density(fr_type_simp_ard)
plot(fr_type_simp_ard_dd)

# Symmetrical model
plot(fr_type_simp_sym_sum, cols, type = "fan", fsize = .01, lwd = .5, cex = c(.4,.2))
add.simmap.legend(colors = cols, x = 0.95*par()$usr[1], y = 0.9*par()$usr[4], prompt = FALSE, fsize = 0.6)

fr_type_simp_sym_dd <- density(fr_type_simp_sym)
plot(fr_type_simp_sym_dd)

## dispersal system ####
runan <- T
if(runan){
  disp_simp_er <- make.simmap(tree_data_disp$phy, getVector(tree_data_disp, Sistema_de_Dispersion), model = "ER", nsim = 1000)
  saveRDS(disp_simp_er, "mk_models/output/disp_simp_er.rds")
  disp_simp_ard <- make.simmap(tree_data_disp$phy, getVector(tree_data_disp, Sistema_de_Dispersion), model = "ARD", nsim = 1000)
  saveRDS(disp_simp_ard, "mk_models/output/disp_simp_ard.rds")
  disp_simp_sym <- make.simmap(tree_data_disp$phy, getVector(tree_data_disp, Sistema_de_Dispersion), model = "SYM", nsim = 1000)
  saveRDS(disp_simp_sym, "mk_models/output/disp_simp_sym.rds")
} else {
  disp_simp_er  <- readRDS("mk_models/output/disp_simp_er.rds")
  disp_simp_ard <- readRDS("mk_models/output/disp_simp_ard.rds")
  disp_simp_sym <- readRDS("mk_models/output/disp_simp_sym.rds")
}

# simmaps summaries
runan <- T
if(runan){
  disp_simp_er_sum  <- summary(disp_simp_er, plot = T)
  saveRDS("mk_models/output/disp_er_sum.rds")
  disp_simp_ard_sum <- summary(disp_simp_ard, plot = F)
  saveRDS("mk_models/output/disp_ard_sum.rds")
  disp_simp_sym_sum <- summary(disp_simp_sym, plot = F)
  saveRDS("mk_models/output/disp_sym_sum.rds")
} else{ 
  disp_simp_er_sum  <- readRDS("mk_models/output/disp_er_sum.rds")
  disp_simp_ard_sum <- readRDS("mk_models/output/disp_ard_sum.rds")
  disp_simp_sym_sum <- readRDS("mk_models/output/disp_sym_sum.rds")
  }

# paint it on the tree
#cols <- setNames(palette()[1:length(levels(getVector(tree_data, Sistema_de_Dispersion)))], 
#                 levels(getVector(tree_data, Sistema_de_Dispersion)))
cols <- setNames(palette()[1:dim(disp_simp_er_sum$ace)[2]],
                 colnames(disp_simp_er_sum$ace))
# Equal Rates
plot(disp_simp_er_sum, cols, type = "fan", fsize = .01, lwd = .5, cex = c(.4,.2))
add.simmap.legend(colors = cols, x = 0.95*par()$usr[1],y = 0.9*par()$usr[4], prompt = FALSE, fsize = 0.6)

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
