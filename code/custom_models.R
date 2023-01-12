library(treeplyr); library(tidyverse)
library(tidytree)
library(ggpubr)

# data load and clean ####
tree <- read.tree("data/tree_mac_s3.tre")
tree$edge.length[tree$edge.length<=0] <- 1e-6
#tree <- extract.clade(tree, 1056) # extrae Spermatophyta
tree$tip.label <- tolower(tree$tip.label)

fruit_dat_raw <- read.csv("data/Base plantas Tinigua modificable_Agosto.csv", header = T)
fruit_dat_raw$ESPECIE <- tolower(fruit_dat_raw$ESPECIE)
fruit_dat_raw$ESPECIE <- gsub(pattern = " ", replacement = "_", x = fruit_dat_raw$ESPECIE)
#fruit_dat_raw <- filter(fruit_dat_raw, !is.na(LARGO_FRUTO.cm.))

fruit_dat <- fruit_dat_raw %>% 
  dplyr::select(ESPECIE,
                LARGO_FRUTO.cm., Fr.width..cm.,
                Largo.hojas.total.cm., Largo.lamina.cm., 
                Area.total.cm2., Area.lámina.cm2.,
                Semillas.Fruto,
                Prom.Largo.semilla..mm., Prom.Ancho.sem..mm.,
                DBH.2017, Densidad..specific.gravity.,
                Tipo.de.hoja,
                Sistema_de_Dispersion = Sistema_de_Dispersión,
                #Color.Janson,
                Habito = Hábito,
                Fruit.Dry.Weight,
                FAMILIA) 

tree_data <- make.treedata(tree, fruit_dat)
tree_data <- tree_data %>% 
  mutate(fr_len = log1p(LARGO_FRUTO.cm.), 
         fr_wd = log1p(Fr.width..cm.),
         lv_len = log1p(Largo.hojas.total.cm.), 
         lam_len = log1p(Largo.lamina.cm.),
         ar_tot = log1p(Area.total.cm2.), 
         ar_lam = log1p(Area.lámina.cm2.),
         sem_fr = log1p(Semillas.Fruto), 
         sem_len = log1p(Prom.Largo.semilla..mm.), 
         sem_wd = log1p(Prom.Ancho.sem..mm.),
         dbh = log1p(DBH.2017), 
         fr_dry_w = log1p(Fruit.Dry.Weight),
         dens = log1p(Densidad..specific.gravity.))

#tree_data <- filter(tree_data, !is.na(fr_len)) 
   #%>% filter(!is.na(sem_fr)) %>%  filter(!is.na(dbh)) %>% 
  #filter(!is.na(sem_len))

tree_data$dat %>% group_by(FAMILIA) %>% 
  summarize(Mean = mean(log(LARGO_FRUTO.cm.), na.rm = TRUE)) %>% 
  arrange(desc(Mean))

tree_data$dat %>% group_by(FAMILIA) %>% 
  summarize(Mean = mean(log1p(Prom.Largo.semilla..mm.), na.rm = TRUE)) %>% 
  arrange(desc(Mean))

tree_data <- tree_data %>% 
  filter(!Sistema_de_Dispersion %in% c("Autocórica", 
                                       "Myrmecocórica",
                                       "myrmecocórica",
                                       "Sinzoocórica/Endozoocórico",
                                       "",
                                       "Anemocórica/Hidrocórica")) %>% 
  filter(!Habito %in% c(
    "Hemiepífito"       ,  "Hierba/Arbusto" , "Palma_de_dosel"     ,
    "Hemiepífito-Arbol" ,  "Hierba_bejucosa", "Palma_de_sotobosque",
    "Epífito"           ,  "Hemiparásita"   , "Árbol_bejucoso"     ,
    "Arbolito/Árbol"    ,  "Arbusto/Hierba" , "Bejuco_epífito"     ,
    "Guadua"            ,  "Hierba/Caña"    , "Palma_acaule"       ,
    "Palma_bejucosa"    ,  "Platanillo"     , "Arbusto/Arbolito"   ,
    "Hierba_parásita"   ,  ""
  )) 

tree_data$dat$Sistema_de_Dispersion <- 
  recode(tree_data$dat$Sistema_de_Dispersion, 
         "Inasistida/Hidrocórica/Anemocórica" = "Anemocórica",
         "Hidrocórica/Dehiscencia_explosiva" = "Dehiscencia_explosiva",
         "Inasistida/Anemocórica" = "Anemocórica",
         "Inasistida/Hidrocórica" = "Hidrocórica")

tree_data$dat$Sistema_de_Dispersion <- 
  recode(tree_data$dat$Sistema_de_Dispersion,  
         'Endozoocórica' = "Endozoochory", 
         "Anemocórica" = "Anemochory", 
         'Inasistida' = "Unassisted", 
         'Dehiscencia_explosiva' = "Explosive dehiscence", 
         'Hidrocórica' = "Hydrochory", 
         'Sinzoocórica' = "Synzoochory",
         "Epizoocórica" = "Epizoochory"
  ) #%>% droplevels()

#tree_data$dat$Habito <- droplevels(tree_data$dat$Habito)
  #recode(tree_data$dat$Habito, 
  #       "Inasistida/Hidrocórica/Anemocórica" = "Anemocórica",
  #       "Hidrocórica/Dehiscencia_explosiva" = "Dehiscencia_explosiva",
  #       "Inasistida/Anemocórica" = "Anemocórica",
  #       "Inasistida/Hidrocórica" = "Hidrocórica")

tree_data$dat %>% 
  group_by(Sistema_de_Dispersion) %>% 
  count() %>% arrange(desc(n))

tree_data$dat %>% 
  group_by(Habito) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n))

# recode dispersal systems #
# dispersal two-categories
tree_data$dat <- tree_data$dat %>% 
  mutate(disp_lato = 
           factor(ifelse(Sistema_de_Dispersion == 
                           "Endozoochory", "Endozoochory", "Non endozoochoric")))

tree_data$dat %>% 
  pull(disp_lato) %>% summary()

# dispersal three-categories
tree_data$dat <- tree_data$dat %>% 
  mutate(disp_3cat = 
           factor(ifelse(Sistema_de_Dispersion == "Endozoochory", "Endozoochoric",
                         ifelse(Sistema_de_Dispersion == "Anemochory", "Anemochoric",
                         "Non endozoochoric"))))
tree_data$dat %>% 
  pull(disp_3cat) %>% summary()

# phylogenetic signal per trait #####
library(phytools)

trait_sig <- list(
  phy_sig_fr_dry_w = filter(tree_data, !is.na(fr_dry_w)) %>% 
    treedply(list("K" = phylosig(phy, getVector(., fr_dry_w), "K"),
                  "lambda" = phylosig(phy, getVector(., fr_dry_w),"lambda"))),
  
  phy_sig_fr = filter(tree_data, !is.na(fr_len)) %>% 
    treedply(list("K" = phylosig(phy, getVector(., fr_len), "K"), 
                  "lambda" = phylosig(phy, getVector(., fr_len),"lambda"))),
  
  phy_sig_ar = filter(tree_data, !is.na(ar_tot)) %>% 
    treedply(list("K" = phylosig(phy, getVector(., ar_tot), "K"), 
                  "lambda" = phylosig(phy, getVector(., ar_tot),"lambda"))),
  
  phy_sig_ln = filter(tree_data, !is.na(lv_len)) %>% 
    treedply(list("K" = phylosig(phy, getVector(., lv_len), "K"), 
                  "lambda" = phylosig(phy, getVector(., lv_len),"lambda"))),
  
  phy_sig_se_l = filter(tree_data, !is.na(sem_len)) %>% 
    treedply(list("K" = phylosig(phy, getVector(., sem_len), "K"),
                  "lambda" = phylosig(phy, getVector(., sem_len),"lambda"))),
  
  phy_sig_se_fr = filter(tree_data, !is.na(sem_fr)) %>% 
    treedply(list("K" = phylosig(phy, getVector(., sem_fr), "K"),
                  "lambda" = phylosig(phy, getVector(., sem_fr),"lambda"))),
  
  phy_sig_se_w = filter(tree_data, !is.na(sem_wd)) %>% 
    treedply(list("K" = phylosig(phy, getVector(., sem_wd), "K"),
                  "lambda" = phylosig(phy, getVector(., sem_wd),"lambda"))),
  
  phy_sig_dens = filter(tree_data, !is.na(dens)) %>% 
    treedply(list("K" = phylosig(phy, getVector(., dens), "K"),
                  "lambda" = phylosig(phy, getVector(., dens),"lambda")))
)

sapply(trait_sig, "[[", 2) 

## stepwise phyloANCOVA: seed isometrics ######
library(phylolm)
library(kableExtra)
# here put the longest dimension in seed length !!
tree_data$dat <- tree_data$dat %>% 
  mutate(sem_len_new = ifelse(sem_len > sem_wd, sem_len, sem_wd),
         sem_wd_new = ifelse(sem_wd < sem_len, sem_wd, sem_len)) 

rownames(tree_data$dat) <- tree_data$phy$tip.label
step_mod_seed <- 
  phylostep(sem_wd_new ~ sem_len_new * disp_3cat + dbh + lv_len + sem_fr , #+ Habito, 
            data = tree_data$dat, phy = tree_data$phy, 
            model = "lambda", 
            direction = "both")

sum_mod_seed <- step_mod_seed %>% summary() 
rownames(sum_mod_seed$coefficients) <- 
  c("Intercept", 
    "Seed\nlength", 
    "Endozoochoric", 
    "Non-endozoochoric",
    "DBH", 
    "Leaf\nlength", 
    "Seeds\nper fruit",
    "Seed length x\nEndozoochoric", 
    "Seed length x\nNon-endozoochoric")

isom_plot <- 
data.frame(
x = rownames(sum_mod_seed$coefficients) %>% as.factor(),  
beta = sum_mod_seed$coefficients[,"Estimate"],
conf = confint(step_mod_seed)
) %>% filter(x != "Intercept") %>% 
  rename(low_ci = "conf.2.5..", high_ci = "conf.97.5..") %>% 
  ggplot(aes(x = stats::reorder(x, -beta), 
             y = beta)) +
  geom_point() +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = .2) + 
  theme_classic(base_size = 12) + 
  theme(legend.position = "none",
        axis.text.x = 
          element_text(angle = 0, vjust = 0.9, hjust = .5, colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  labs(x = "", y = "β coefficient", subtitle = "Seed width")

sum_mod_seed$coefficients %>% 
  round(3) %>% data.frame() %>% #select(-StdErr) %>% 
  rename(Beta = Estimate, SE = StdErr, 
         "t-value" = t.value, "p-value" = p.value) %>%  
  kbl(caption = paste0("Table 1S. Seed isometry model coefficients") ) %>%
  kable_classic(html_font = "Cambria", #font_size = 20, 
              full_width = F, position = "float_right") %>% 
  save_kable(file = "figures/seed_step_fit.png"
             , self_contained = T)

### seed elongation
tree_data <- tree_data %>% 
  mutate(sem_elong = (sem_wd_new/sem_len_new))

rownames(tree_data$dat) <- tree_data$phy$tip.label

fit_rd <- phylolm(sem_elong ~ sem_fr + disp_3cat, #+ Habito, 
                  data = tree_data$dat, phy = tree_data$phy, 
                  model = "lambda")

sum_mod_elong <- summary(fit_rd)

sum_mod_elong$coefficients %>% 
  round(3) %>% data.frame() %>% #select(-StdErr) %>% 
  rename(Beta = Estimate, SE = StdErr, 
         "t-value" = t.value, "p-value" = p.value) %>%  
  kbl(caption = paste0("Table 3S. Seed elongation model coeficients")) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>% 
  save_kable(file = "figures/seed_elong_fit.png"
             , self_contained = T)

## stepwise phyloANCOVA: fruit length vs leaf feats #####
rownames(tree_data$dat) <- tree_data$phy$tip.label
step_mod_fr <- 
phylostep(fr_len ~ lv_len + sem_len + sem_wd + dbh + sem_fr + disp_3cat, #+ Habito, 
             data = tree_data$dat, phy = tree_data$phy, 
             model = "lambda", 
             direction = "both")

sum_mod_fr <- step_mod_fr %>% summary() 
rownames(sum_mod_fr$coefficients) <- 
  c("Intercept", "Leaf\nlength",
    "Seed\nlength", "Seed\nwidth",
    "DBH", 
    "Seeds\nper fruit",
    "Endozoochoric", 
    "Non-endozoochoric")

fruit_plot <- 
data.frame(
  x = rownames(sum_mod_fr$coefficients) %>% as.factor(),  
  beta = sum_mod_fr$coefficients[,"Estimate"],
  conf = confint(step_mod_fr)
) %>% filter(x != "Intercept") %>% 
  rename(low_ci = "conf.2.5..", high_ci = "conf.97.5..") %>% 
  ggplot(aes(x = stats::reorder(x, -beta), 
             y = beta)) +
  geom_point() +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = .2) + 
  theme_classic(base_size = 12) + 
  theme(legend.position = "none",
        axis.text.x = 
          element_text(angle = 0, vjust = 0.9, hjust = .5, colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  labs(x = "", y = "β coefficient", subtitle = "Fruit length")

sum_mod_fr$coefficients %>% 
  round(3) %>% data.frame() %>% #select(-StdErr) %>% 
  rename(Beta = Estimate, SE = StdErr, 
         "t-value" = t.value, "p-value" = p.value) %>%  
  kbl(caption = paste0("Table 2S. Fruit length model coefficients")) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>% 
  save_kable(file = "figures/fru_step_fit.png", density = 600,
             self_contained = T)

library(patchwork)
isom_plot + fruit_plot +
  plot_layout(design = "A
              B", guides = 'collect') +
  plot_annotation(tag_levels = 'a', tag_sep = " ")
ggsave("figures/eff_ci.png", 
       width = 10, height = 12, dpi = 600)

# simple anovas: phylom does not work !! ?? 
# Largo.hojas per family
lm_f <- phylolm(log1p(Largo.hojas.total.cm.) ~ as.factor(FAMILIA), 
                data = tree_data$dat, phy = tree_data$phy, 
                model = "lambda")

tab <- lm(log1p(Largo.hojas.total.cm.) ~ FAMILIA, 
          data = tree_data$dat) %>% summary()

tab$coefficients %>% data.frame() %>% # select(Pr...t..) %>% 
  filter(Pr...t.. <= .05) %>% arrange((Pr...t..)) %>% 
  mutate(fams = rownames(.) %>% gsub("FAMILIA", "", .)) %>% 
  dplyr::select(fams, Estimate, SE = Std..Error, t.value, p.val = Pr...t..) %>% 
  write.csv(file = "largo.hojas.anova.csv", row.names = F)

tab2 <- lm(fr_len ~ FAMILIA, data = tree_data$dat)

summary(tab2)
anova(tab2)

# laminar area vs DBH
lm_m <- phylolm(log1p(Area.lámina.cm2.) ~ log1p(DBH.2017), 
                data = tree_data$dat, phy = tree_data$phy, 
                model = "lambda")
lm_m %>% summary()
lm_m$n
coefficients(lm_m)[2] %>% exp() -1 # coeff in real space

# individual models #####
# seed isometry
seed_isom_3disp <- phylolm(sem_wd_new ~ sem_len_new * disp_3cat, 
                           data = tree_data$dat, phy = tree_data$phy, 
                           model = "lambda")
summary(seed_isom_3disp)
emmeans::emmeans(seed_isom_3disp, ~ disp_3cat)

seed_isom_3disp_lm <- lm(sem_wd_new ~ sem_len_new * disp_3cat, 
                           data = tree_data$dat)


seed_iso <- phylolm(sem_len_new ~ sem_wd_new, 
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

seed_disp_dbh <- phylolm(sem_len_new ~ dbh + sem_wd_new * disp_3cat, 
                         data = tree_data$dat, phy = tree_data$phy, 
                         model = "lambda")

models_seed <- list(seed_iso, seed_iso_disp, 
                    seed_isoxdisp, seed_isoxdisp_sat
                    , seed_disp_dbh
)

data.frame(
  "Model_formula" = sapply(models_seed, '[[', 16, simplify = T) %>% as.character(),
  logLik = sapply(models_seed, '[[', 5),
  lambda = sapply(models_seed, '[[', 3),
  #sigma = sapply(models_seed, '[[', 2),
  r2 = sapply(models_seed, '[[', 20),
  geiger::aicw(sapply(models_seed, '[[', 7))
  #,AICw = c(phytools::aic.w(sapply(models_ls, '[[', 3)))
) %>% rename(AIC = fit, ΔAIC = delta, AICw = w) %>% 
  arrange(desc(AICw)) %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  kbl(caption = "Seed isometry models") %>%
  kable_classic(full_width = F, html_font = "Cambria") #%>% 
save_kable(file = "figures/seed_mod_fits.pdf"
           , self_contained = T)

summary(seed_disp_dbh)

# fruit size
rownames(tree_data$dat) <- tree_data$phy$tip.label

fr_len_m1 <- phylolm(formula = fr_len ~ lv_len, 
                     data = tree_data$dat, phy = tree_data$phy, 
                     model = "lambda")

fr_len_m2 <- phylolm(fr_len ~ lv_len * disp_3cat, 
                     data = tree_data$dat, phy = tree_data$phy, 
                     model = "lambda")

fr_len_m3 <- phylolm(fr_len ~ lv_len * disp_3cat + sem_len_new, 
                     data = tree_data$dat, phy = tree_data$phy, 
                     model = "lambda")

fr_len_m4 <- phylolm(fr_len ~ lv_len * disp_3cat + dbh + sem_len_new, 
                     data = tree_data$dat, phy = tree_data$phy, 
                     model = "lambda")

#fr_len_m5 <- phylolm(fr_len ~ ar_tot * disp_3cat, 
#                data = tree_data$dat, phy = tree_data$phy, 
#                model = "lambda")

#fr_len_m6 <- phylolm(fr_len ~ ar_tot + disp_3cat, 
#                data = tree_data$dat, phy = tree_data$phy, 
#                model = "lambda")

fr_len_m7 <- phylolm(fr_len ~ lv_len * disp_3cat + dbh, 
                     data = tree_data$dat, phy = tree_data$phy, 
                     model = "lambda")

fr_len_m8 <- phylolm(fr_len ~ lv_len + sem_len_new + dbh, 
                     data = tree_data$dat, phy = tree_data$phy, 
                     model = "lambda")

fr_len_m9 <- phylolm(fr_len ~ lv_len + disp_3cat + dbh, 
                     data = tree_data$dat, phy = tree_data$phy, 
                     model = "lambda")

models_ls <- list(#fr_len_m1, 
  fr_len_m2, fr_len_m3, 
  fr_len_m4
  #, fr_len_m5, fr_len_m6,
  , fr_len_m7 , fr_len_m8, fr_len_m9
)

# models summary table
data.frame(
  "Model_formula" = sapply(models_ls, '[[', 16, simplify = T) %>% 
    as.character(),
  logLik = sapply(models_ls, '[[', 5),
  lambda = sapply(models_ls, '[[', 3),
  #sigma = sapply(models_ls, '[[', 2),
  r2 = sapply(models_ls, '[[', 20),
  geiger::aicw(sapply(models_ls, '[[', 7))
  #,AICw = c(phytools::aic.w(sapply(models_ls, '[[', 3)))
) %>% 
  rename(AIC = fit, ΔAIC = delta, AICw = w) %>%  
  arrange(desc(AICw)) %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  kbl(caption = "Fruit length models") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>% 
  save_kable(file = "figures/fru_mod_fits.pdf"
             , self_contained = T)

fr_len_m4 %>% summary()
fr_len_m4 %>% coefficients()

fr_len_m2 %>% summary()
fr_len_m2 %>% coefficients()

# plot predictions and residuals
library(modelr)
grid <- tree_data$dat %>% 
  filter(!is.na(lv_len), !is.na(disp_3cat), !is.na(dbh)) %>%  
  data_grid(lv_len, disp_3cat, dbh) %>% 
  gather_predictions(step_mod_fr, lm_reg, lm_reg_disp)

tree_data$dat %>% 
  ggplot(aes(lv_len, fr_len, shape = disp_3cat, colour = disp_3cat)) + 
  geom_point() + 
  geom_line(data = grid, aes(y = pred)) + 
  theme_classic() +
  scale_color_viridis_d(option = "viridis", direction = -1, begin = 0, end = .65) +
  facet_wrap(~ model)

tree_data$dat %>% filter(!is.na(lv_len), !is.na(disp_3cat), !is.na(dbh)) %>% 
  gather_residuals(fr_len_m2, lm_reg, lm_reg_disp) %>% 
  ggplot(aes(lv_len, resid, colour = disp_3cat)) + 
  geom_point() +
  scale_color_viridis_d(option = "viridis", direction = -1, begin = 0, end = .65) +
  facet_grid(model ~ disp_3cat)
  
# stepwise with leaf area
#ml_step <- phylostep(fr_len ~ lv_len + ar_tot + dens,
#                     data = tree_data$dat, phy = tree_data$phy, 
#                     model = "lambda", direction = "both")


# OU shifts in a traits
#OUshifts(getVector(tree_data, fr_len), tree_data$phy,
#         method = "mbic", nmax = 2)


# laminar area vs plant habit
# habits summary
tree_data$dat %>% 
  group_by(Habito) %>% count() %>% 
  write_delim(file = "habitos_mac.csv", delim = ";")

# habits load
hab_df <- read.delim("data/habitos_mac.csv", sep = ";")
tree_data$dat <- 
left_join(
  tree_data$dat,
  hab_df %>% dplyr::select(Habito, new_habito = X), 
  by = c("Habito")
)

tree_data$dat$new_habito <- recode(tree_data$dat$new_habito, "Palma" = "Arbol") 
tree_data$dat %>% 
  group_by(new_habito) %>% count()

tree_data$dat %>% filter(new_habito != "") %>%  
  ggplot(aes(x = new_habito, y = log(Area.lámina.cm2.))) +
  geom_boxplot() + theme_classic() +
  labs(x = "Habit", y = "Ln Laminar area (cm2)")

lm_h <- phylolm(log(Area.lámina.cm2.) ~ new_habito, 
                data = tree_data$dat, phy = tree_data$phy, model = "lambda")
lm_h %>% summary()
lm_h$n
coefficients(lm_h)[9] %>% exp() -1 # coeff in real space

# pgls fruit length vs leaf length or area #
tr_dt_fr_lv <- filter(tree_data, !is.na(fr_len)) %>% 
  filter(!is.na(lv_len))

tr_dt_fr_ar <- filter(tree_data, !is.na(fr_len)) %>% 
  filter(!is.na(ar_tot))
                                                            
bm_ln <- gls(fr_len ~ lv_len, 
             correlation = corBrownian(phy = tr_dt_fr_lv$phy),
               data = tr_dt_fr_lv$dat, method = "ML")

bm_ar <- gls(fr_len ~ ar_tot, 
             correlation = corBrownian(phy = tr_dt_fr_ar$phy), 
               data = tr_dt_fr_ar$dat, method = "ML")

bm_ln_disp <- gls(fr_len ~ lv_len*disp_lato, 
                  correlation = corBrownian(phy = tr_dt_fr_lv$phy), 
               data = tr_dt_fr_lv$dat, method = "ML")
bm_ar_disp <- gls(fr_len ~ ar_tot*disp_lato, 
                  correlation = corBrownian(phy = tr_dt_fr_ar$phy), 
               data = tr_dt_fr_ar$dat, method = "ML")

data.frame(
length = coef(bm_ln_disp),
area = coef(bm_ar_disp)
) %>% round(3)

tr_dt_fr_lv_tmp <- tr_dt_fr_lv$phy
tr_dt_fr_ar_tmp <- tr_dt_fr_ar$phy

tr_dt_fr_lv_tmp$edge.length <- tr_dt_fr_lv_tmp$edge.length*100
tr_dt_fr_ar_tmp$edge.length <- tr_dt_fr_ar_tmp$edge.length*100

ou_ln <- gls(fr_len ~ lv_len, 
             correlation = corMartins(1, phy = tr_dt_fr_lv_tmp), 
             data = tr_dt_fr_lv$dat, method = "ML")

ou_ar <- gls(fr_len ~ ar_tot, 
             correlation = corMartins(1, phy = tr_dt_fr_ar_tmp), 
             data = tr_dt_fr_ar$dat, method = "ML")


ou_ln_disp <- gls(fr_len ~ lv_len*disp_lato, 
                  correlation = corMartins(1, phy = tr_dt_fr_lv_tmp), 
                  data = tr_dt_fr_lv$dat, method = "ML")

ou_ar_disp <- gls(fr_len ~ ar_tot*disp_lato, 
                  correlation = corMartins(1, phy = tr_dt_fr_ar_tmp), 
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


