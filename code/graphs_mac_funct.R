library(tidyverse) ; library(treeplyr)
library(tidytree)
library(ggpubr)
library(wesanderson)

# data load and clean ####
tree <- read.tree("data/tree_mac_s3.tre")
tree$edge.length[tree$edge.length<=0] <- 1e-6
#tree <- extract.clade(tree, 1056) # extrae Spermatophyta
tree$tip.label <- tolower(tree$tip.label)

fruit_dat_raw <- read.csv("data/Base plantas Tinigua modificable_Agosto.csv", header = T)
fruit_dat_raw$ESPECIE <- tolower(fruit_dat_raw$ESPECIE)
fruit_dat_raw$ESPECIE <- gsub(pattern = " ", replacement = "_", x = fruit_dat_raw$ESPECIE)
#fruit_dat_raw <- filter(fruit_dat_raw, !is.na(LARGO_FRUTO.cm.))

fruit_dat <- fruit_dat_raw %>% dplyr::select(ESPECIE,
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
tree_data <- tree_data %>%  mutate(fr_len = log1p(LARGO_FRUTO.cm.), fr_wd = log1p(Fr.width..cm.),
                                   lv_len = log1p(Largo.hojas.total.cm.), lam_len = log1p(Largo.lamina.cm.),
                                   ar_tot = log1p(Area.total.cm2.), ar_lam = log1p(Area.lámina.cm2.),
                                   sem_fr = log1p(Semillas.Fruto), 
                                   sem_len = log1p(Prom.Largo.semilla..mm.), sem_wd = log1p(Prom.Ancho.sem..mm.),
                                   dbh = log1p(DBH.2017), 
                                   fr_dry_w = log1p(Fruit.Dry.Weight),
                                   dens = log1p(Densidad..specific.gravity.))

tree_data <- tree_data %>% 
  filter(!Sistema_de_Dispersion %in% c("Autocórica", 
                                       "Myrmecocórica",
                                       "myrmecocórica",
                                       "Sinzoocórica/Endozoocórico",
                                       "Anemocórica/Hidrocórica",
                                       ""))

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

tree_data$dat %>% 
  group_by(Sistema_de_Dispersion) %>% 
  count() %>% arrange(desc(n))

# Figure 1, plot phylo with data and panel hist vars####
library(ggtree)
library(ggnewscale)

#library(TDbook) ; library(ggimage)

tree_data$dat %>% 
  group_by(FAMILIA) %>% 
  summarize(n = n()) %>% arrange(desc(n)) %>% 
  mutate(tog = paste0(FAMILIA," (", n,")")) %>% 
  slice(1:10) %>% pull(tog)

# dispersal three-categories
tree_data$dat <- tree_data$dat %>% 
  mutate(disp_3cat = 
           factor(ifelse(Sistema_de_Dispersion == "Endozoochory", "Endozoochoric",
                         ifelse(Sistema_de_Dispersion == "Anemochory", "Anemochoric",
                                "Non endozoochoric"))))
tree_data$dat %>% 
  pull(disp_3cat) %>% summary()

# basic tree
circ <- ggtree(tree_data$phy, 
               layout = "circular",
               size = 0.15) #+ 
  #geom_text(aes(label = node), color = "red", cex = 2.5, hjust = -.3) # red node numbers
# %>% open_tree(angle = 45) # fish eye

#add dispersal systems
p1 <- gheatmap(circ, 
               tree_data$dat %>% 
                 dplyr::select(disp_3cat) %>% 
                 data.frame(row.names = tree_data$phy$tip.label), 
               offset = 90, width = .1,
               colnames = F) +
  scale_fill_manual(name = "Dispersal\nsystem", 
                    values = wes_palette("FantasticFox1",  3, type = "discrete"))
p2 <- p1 + new_scale_fill()

p3 <- 
  gheatmap(p2, 
           tree_data$dat %>% 
             dplyr::select(fr_len, lv_len, sem_len) %>% 
             data.frame(row.names = tree_data$phy$tip.label), 
           offset = 5, width = .25,
           colnames = F
           #, custom_column_labels = c("A", "B", "C"),
           #colnames_angle = 90, colnames_offset_y = .25
  ) +
  scale_fill_gradientn(name = "Functional\ntraits", 
                       colours = wes_palette("Zissou1",  100, type = "continuous")) + 
  scale_fill_viridis_c(option = "D", name = "Functional\ntraits")

# families to paint
fam2paint <- c("Fabaceae", "Rubiaceae", "Moraceae", 
               "Araceae",  "Arecaceae", "Melastomataceae", 
               "Apocynaceae", "Bignoniaceae", "Malvaceae")

fams_nodes <- 
lapply(fam2paint, 
       function(x) castor::get_mrca_of_set(tree_data$phy, 
                        filter(tree_data, FAMILIA == x)$phy$tip.label)
) %>% unlist()

# shading clades
phylo_data <- 
p3 + #theme(legend.position = "top") +
  geom_hilight(node = fams_nodes[1], "min", fill = 'darkgreen', alpha = .5) +
  geom_hilight(node = fams_nodes[2], 'min', fill = 'firebrick', alpha = .5) +
  geom_hilight(node = fams_nodes[3], 'min', fill = 'salmon',    alpha = .65) +
  geom_hilight(node = fams_nodes[4], 'min', fill = 'goldenrod', alpha = .5) +
  geom_hilight(node = fams_nodes[5], 'min', fill = 'steelblue', alpha = .5) +
  geom_hilight(node = fams_nodes[6], 'min', fill = 'pink',      alpha = .5) +
  geom_hilight(node = fams_nodes[7], 'min', fill = 'blue',      alpha = .5) +
  geom_hilight(node = fams_nodes[8], 'min', fill = 'cyan2',      alpha = .5) +
  geom_hilight(node = fams_nodes[9], 'min', fill = 'violetred',  alpha = .5)
ggsave("figures/tree_data.png", 
       width = 25, height = 25, units = "cm", dpi = 600)

# other ways to show clades in trees
# collapsing clades
p3 %>% 
  ggtree::collapse(node = fams_nodes[1], "mixed", fill = 'darkgreen', alpha = .5) %>% 
  ggtree::collapse(node = fams_nodes[2], 'mixed', fill = 'firebrick', alpha = .5) %>%
  ggtree::collapse(node = fams_nodes[7], 'mixed', fill = 'steelblue', alpha = .5)

# bar and letters on clades
#p3 +
circ + 
  geom_cladelabel(node = fams_nodes[1], label = fam2paint[1], 
                  fontsize = 2, angle = 90, hjust = "center",
                  geom = "label", fill = "white",
                  #barsize = .001, 
                  color = "darkgreen", offset.text = 5, align = T) +
  geom_cladelabel(node = fams_nodes[2], label = fam2paint[2], 
                  fontsize = 2, angle = 90, hjust = "center",
                  geom = "label", fill = "white",
                  #barsize = .001, 
                  color = "salmon", offset.text = 5, align = T) +
  geom_cladelabel(node = fams_nodes[3], label = fam2paint[3], 
                  fontsize = 2, angle = 90, hjust = "center",
                  geom = "label", fill = "white",
                  #barsize = .001, 
                  color = "firebrick", offset.text = 5, align = T) +
  geom_cladelabel(node = fams_nodes[4], label = fam2paint[4], 
                  fontsize = 2, angle = 90, hjust = "center",
                  geom = "label", fill = "white",
                  #barsize = .001, 
                  color = "goldenrod", offset.text = 5, align = T) +
  geom_cladelabel(node = fams_nodes[5], label = fam2paint[5], 
                  fontsize = 2, angle = 90, hjust = "center",
                  geom = "label", fill = "white",
                  #barsize = .001, 
                  color = "steelblue", offset.text = 5, align = T) +
  geom_cladelabel(node = fams_nodes[6], label = fam2paint[6], 
                  fontsize = 2, angle = 90, hjust = "center",
                  geom = "label", fill = "white",
                  #barsize = .001, 
                  color = "pink", offset.text = 5, align = T) +
  geom_cladelabel(node = fams_nodes[7], label = fam2paint[7], 
                  fontsize = 2, angle = 90, hjust = "center",
                  geom = "label", fill = "white",
                  #barsize = .001, 
                  color = "blue", offset.text = 5, align = T) +
  geom_cladelabel(node = fams_nodes[8], label = fam2paint[8], 
                  fontsize = 2, angle = 90, hjust = "center",
                  geom = "label", fill = "white",
                  #barsize = .001, 
                  color = "cyan2", offset.text = 5, align = T)+
  geom_cladelabel(node = fams_nodes[9], label = fam2paint[9], 
                  fontsize = 2, angle = 90, hjust = "center",
                  geom = "label", fill = "white",
                  #barsize = .001, 
                  color = "violetred", offset.text = 5, align = T) +
  geom_hilight(node = fams_nodes[1], "min", fill = 'darkgreen', alpha = .5) +
  geom_hilight(node = fams_nodes[2], 'min', fill = 'firebrick', alpha = .5) +
  geom_hilight(node = fams_nodes[3], 'min', fill = 'salmon',    alpha = .65) +
  geom_hilight(node = fams_nodes[4], 'min', fill = 'goldenrod', alpha = .5) +
  geom_hilight(node = fams_nodes[5], 'min', fill = 'steelblue', alpha = .5) +
  geom_hilight(node = fams_nodes[6], 'min', fill = 'pink',      alpha = .5) +
  geom_hilight(node = fams_nodes[7], 'min', fill = 'blue',      alpha = .5) +
  geom_hilight(node = fams_nodes[8], 'min', fill = 'cyan2',      alpha = .5) +
  geom_hilight(node = fams_nodes[9], 'min', fill = 'violetred',  alpha = .5)

# check
tree_data$dat %>%
  select(FAMILIA,Sistema_de_Dispersion, Habito, disp_3cat) %>% 
  data.frame(row.names = tree_data$phy$tip.label) %>% tail()

library(ggthemes)

mu <- 
  tree_data$dat %>% 
  select("Log Fruit length (cm)" = LARGO_FRUTO.cm., 
         "Log Leaf length (cm)" = Largo.hojas.total.cm., 
         "Log Leaf area (cm2)" = Area.total.cm2.,
         "Log Seeds per fruit" = Semillas.Fruto, 
         "Log Seed length (mm)" = Prom.Largo.semilla..mm.,
         "Log Fruit dry weight (gr)" = "Fruit.Dry.Weight") %>% 
  mutate_all(log1p) %>% 
  gather(key = trait, value) %>%
  plyr::ddply("trait", summarise, 
              grp.mean = mean(value, na.rm = T),
              sd_d = sd(value, na.rm = T),
              cv = (sd_d/grp.mean)*100,
              l_sd = grp.mean-sd_d,
              u_sd = grp.mean+sd_d) %>% 
  mutate_if(is.numeric, round, 2)

paste0(
  "λ = ", 
  traits_phylo_sig_df %>%
    slice(1:6) %>% 
    pull(lambda) %>% 
    unlist() %>% round(2)
)

paste0(
  "p = ", 
  traits_phylo_sig_df %>%
    slice(1:6) %>% 
    pull(`p-value`) %>% 
    unlist() #%>% round(3)
)

lbd <- 
  data.frame(
    trait = c("Log Fruit dry weight (gr)",
              "Log Fruit length (cm)", 
              "Log Leaf area (cm2)", 
              "Log Leaf length (cm)", 
              "Log Seed length (mm)",
              "Log Seeds per fruit"),
    lambda_txt = paste0(
      "λ = ", 
      traits_phylo_sig_df %>%
        slice(1:6) %>% 
        pull(lambda) %>% 
        unlist() %>% round(2),
      " ***"
    ),
    x_p = c(5, 4, 9, 5.5, 3.7, 7),
    y_p  = c(42, 90, 72, 148, 77, 300)
  )

hist_traits_plot <-   
  tree_data$dat %>% 
  select("Log Fruit length (cm)" = LARGO_FRUTO.cm., 
         "Log Leaf length (cm)" = Largo.hojas.total.cm.,
         "Log Leaf area (cm2)" = Area.total.cm2.,
         "Log Seeds per fruit" = Semillas.Fruto, 
         "Log Seed length (mm)" = Prom.Largo.semilla..mm.,
         "Log Fruit dry weight (gr)" = "Fruit.Dry.Weight") %>% 
  gather(key = trait, value) %>%
  ggplot(aes(x = log1p(value))) +
  geom_histogram(color = "gray", fill = "lightgray",
                 position = "identity", na.rm = T) +
  theme_classic(base_size = 12) + 
  geom_vline(data = mu, aes(xintercept = grp.mean),
             color = "red", linetype = "dashed", size = .65) +
  geom_vline(data = mu, aes(xintercept = l_sd),
             color = "black", linetype = "longdash", size = .5) +
  geom_vline(data = mu, aes(xintercept = u_sd),
             color = "black", linetype = "longdash", size = .5) +
  geom_segment(data = data.frame(trait = "Log Fruit dry weight (gr)", 
                                 xvalue = 4.63, yvalue = 10, 
                                 xend = 4.63, yend = 5,
                                 lineend = "mitre",
                                 linejoin = "mitre"),
               aes(x = xvalue, y = yvalue, xend = xend, yend = yend),
               arrow = arrow(length = unit(.25, "cm"), type = "open"),
               size = 1, colour = "steelblue") +
  geom_text(data = lbd, 
            mapping = aes(x = x_p, y = y_p, label = lambda_txt)) +
  facet_wrap(~trait, scales = "free", 
             #switch = "x"
             #strip.position = c("bottom")
  ) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(y = "Count", x = "")
ggsave("figures/vars_hists_log1p.png", bg = "white",
       width = 25, height = 25, units = "cm", dpi = 600)

# Fig. 1 ab verical
library(patchwork)
phylo_data + hist_traits_plot +
  plot_layout(design = "AA
                        AA
                        BB") +
  plot_annotation(tag_levels = 'a', tag_sep = " ")
ggsave("figures/fig_1ab.png", 
       width = 28, height = 35, units = "cm", dpi = 600)

# Figure 2, seed width vs length by dispersal system ####
tree_data$dat %>% 
  ggplot(aes(x = exp(sem_len_new), y = exp(sem_wd_new), 
             shape = disp_3cat,
             color = disp_3cat)) + 
  geom_point() + theme_classic(base_size = 12) + 
  geom_smooth( # "formula = log(y) ~ (x)"
    method = "lm", 
    se = T, show.legend = F) +
  theme(legend.position = "top") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "Seed length (mm)", y = "Seed width (mm)", 
       color = "Dispersal system", shape = "Dispersal system") +
  scale_color_manual(values = wes_palette("FantasticFox1",  3, type = "continuous"))
  #scale_color_viridis_d(option = "magma", direction = 1, begin = 0, end = .65) 
#geom_abline(slope = .97, intercept = -0.49, lwd = 1, col = "#488f31") +
#geom_abline(slope = .86, intercept = -.76, lwd = 1, col = "#e08861") 
ggsave("figures/seed_isom_3cat.png", 
       width = 20, height = 20, units = "cm", dpi = 600)

# regression formula per category
paste0("Y_endo = ", round(coefficients(seed_isoxdisp)[1], 2), " + ", round(coefficients(seed_isoxdisp)[2],2), "x")
paste0("y_no = ", round(coefficients(seed_isoxdisp)[1]+coefficients(seed_isoxdisp)[3], 2), " + ", round(coefficients(seed_isoxdisp)[2]+coefficients(seed_isoxdisp)[4],2), "x")

# with density distributions
p <- tree_data$dat %>% 
  ggplot(aes(x = exp(sem_len_new), 
             y = exp(sem_wd_new),
             #shape = Sistema_de_Dispersion,
             color = Sistema_de_Dispersion)) +
  geom_point() + theme_classic() +
  theme(legend.position = "top") +
  labs(x = "Seed length (mm)", y = "Seed width (mm)") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_smooth(#formula = (y) ~ log1p(x)
    method = "lm"
    )
ggExtra::ggMarginal(p, type = "density", 
                    groupColour = T, groupFill = T)

# isometry # other figs
# seed width ~ seed length
fruit_dat_raw %>% 
  ggplot(aes(x = log(Prom.Largo.semilla..mm.), y = log(Prom.Ancho.sem..mm.))) +
  geom_point() + theme_classic(base_size = 10) + 
  geom_smooth(method = "lm", se = T, show.legend = F, col = "red") +
  theme(legend.position = "none") +
  #geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "Ln Mean Seed length (mm)", y = "Ln Mean Seed width (mm)") +
  scale_color_viridis_d(option = "viridis", direction = -1, begin = 0, end = .65)
ggsave("figures/seed_len-seed_widt.eps", width = 8.9, height = 8.9, units = "cm")

# scatter for log-log with regression info
ggscatter(tree_data$dat, x = "sem_len", y = "sem_wd", 
          color = "disp_3cat", add = "reg.line", palette = "Set1",
          xlab = "Ln Seed length (mm)", ylab = "Ln Seed width (mm)", 
          shape = "disp_3cat", fullrange = F) +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), 
        color = disp_3cat)) + 
  geom_abline(slope = 1, intercept = 0, lty = 2)

# seed data transposition 
tree_data$dat <- tree_data$dat %>% 
  mutate(sem_len_new = ifelse(sem_len > sem_wd, sem_len, sem_wd),
         sem_wd_new = ifelse(sem_wd < sem_len, sem_wd, sem_len)) 

ggscatter(tree_data$dat, x = "sem_len_new", y = "sem_wd_new", 
          color = "disp_3cat", add = "reg.line", palette = "Set1",
          xlab = "Seed length (mm)", ylab = "Seed width (mm)", 
          shape = "disp_3cat", fullrange = F) +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~"), 
        color = disp_3cat)) + 
  geom_abline(slope = 1, intercept = 0, lty = 2)

# y ~ logp1(x) with regression pars
ggplot(tree_data$dat, aes(x = exp(sem_len_new), y = exp(sem_wd_new), 
                          color = disp_3cat, shape = disp_3cat)) +
  geom_point() + theme_classic() +
  stat_smooth(aes(fill = disp_3cat, color = disp_3cat),
              formula = y ~ log1p(x), fullrange = F) +
  stat_regline_equation(formula = y ~ log1p(x), 
                        aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), 
                            color = disp_3cat)) + 
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "Seed length (mm)", y = "Seed width (mm)") +
  theme(legend.position = "top", panel.background = element_blank())

# Figure 3, fruit length ~ leaf length by dispersal system ####
# log-transf datavars
tree_data$dat %>% 
  ggplot(aes(x = (lv_len), y = (fr_len), 
             #shape = Sistema_de_Dispersion, 
             color = disp_3cat)) +
  geom_point() + 
  geom_smooth(method = "lm", 
              se = T, fullrange = F, alpha = .2) + 
  theme_classic(base_size = 12) +
  theme(legend.position = "top") + 
  labs(x = "Ln Leaf length (cm)", y = "Ln Fruit length (cm)",
       color = "Dispersal system", shape = "Dispersal system") +
  scale_color_manual(values = wes_palette("FantasticFox1",  3, type = "continuous"))
ggsave("figures/fr_lf_len_3cat.png", 
       width = 20, height = 20, units = "cm", dpi = 600)

# other figs
# scatter with regression info per dispersal system
ggscatter(tree_data$dat, x = "Largo.hojas.total.cm.", y = "LARGO_FRUTO.cm.", 
          color = "disp_3cat", add = "reg.line", palette = "Set1",
          xlab = "Leaf length (cm)", ylab = "Fruit length (cm)", 
          shape = "disp_3cat", fullrange = F) +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), 
        color = disp_3cat)) 
#+ geom_abline(slope = 1, intercept = 0, lty = 2)

# scatter with regression info per category, log-transformed
ggscatter(tree_data$dat, x = "lv_len", y = "fr_len", 
          color = "disp_3cat", add = "reg.line", palette = "Set1",
          xlab = "Ln Leaf length (cm)", ylab = "Ln Fruit length (cm)", 
          shape = "disp_3cat", fullrange = F) + 
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), 
        color = disp_3cat))

# Figure 4. average fruit and laminar length ####
# biggest 10 families
fams_av <- fruit_dat_raw %>% group_by(FAMILIA) %>% 
  summarise(ave_fr = mean(LARGO_FRUTO.cm., na.rm = T), n = n()) %>% 
  arrange(desc(ave_fr)) %>% filter(n>=3) %>% pull(FAMILIA)

g1 <- 
  fruit_dat_raw %>% 
  ggplot(aes(x = (FAMILIA), y = (LARGO_FRUTO.cm.))) +
  geom_boxplot() + 
  theme_classic(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  labs(x = "", y = "Fruit length (cm)") + 
  scale_x_discrete(limits = fams_av[1:10])

# average leave length by family
fams_lav <- fruit_dat_raw %>% group_by(FAMILIA) %>% 
  summarise(ave_lv = mean(Largo.lamina.cm., na.rm = T), n = n()) %>% 
  arrange(desc(ave_lv)) %>% filter(n>=3) %>% pull(FAMILIA)

g2 <- 
  fruit_dat_raw %>% 
  ggplot(aes(x = (FAMILIA), y = (Largo.lamina.cm.))) +
  geom_boxplot() + 
  theme_classic(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  labs(x = "", y = "Leaf lamina length (cm)") + 
  scale_x_discrete(limits = fams_lav[1:10])

ggarrange(g1, g2, 
          ncol = 1, nrow = 2, 
          labels = c("a", "b"),
          legend = "top")
ggsave("figures/big_aver_panel.png", 
       width = 21, height = 20, units = "cm", dpi = 600)

# Figure X Supp. panel with leaf length and area by dispersal system ####
# fruit length ~ leaf length 
pp1 <- tree_data$dat %>% 
  ggplot(aes(x = log(Largo.hojas.total.cm.), 
             y = log(LARGO_FRUTO.cm.))) + 
  geom_point(cex = .8) + theme_classic(base_size = 12) + 
  #geom_smooth(method = "lm", se = T, show.legend = F, col = "#67017E") +
  #geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme(legend.position = "none") +
  labs(x = "Ln Leaf length (cm)", y = "Ln Fruit length (cm)")

# fruit length ~ total leaf area 
pp2 <- tree_data$dat %>% 
  ggplot(aes(x = (ar_tot), y = (fr_len), color = disp_3cat, shape = disp_3cat)) + 
  geom_point() + theme_classic(base_size = 12) + 
  geom_smooth(method = "lm", se = T, show.legend = F) +
  theme(legend.position = "top") +
  #geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "Ln Leaf area (mm2)", y = "", 
       color = "Dispersal system", shape = "Dispersal system") +
  scale_color_manual(values = wes_palette("FantasticFox1",  3, type = "continuous"))

ggarrange(pp1, pp2,
          ncol = 2, nrow = 1,
          labels = c("a", "b"),
          legend = "top", common.legend = TRUE)
ggsave("figures/phylo_lm_leaf_seed.png", width = 18.4, height = 8.9, units = "cm")

# fruit length ~ seed length
pp3 <- 
  fruit_dat_raw %>% ggplot(aes(x = log(Prom.Largo.semilla..mm.), 
                               y = log(LARGO_FRUTO.cm.))) +
  geom_point(cex = .8) + theme_classic(base_size = 12) + 
  geom_smooth(method = "lm", se = T, show.legend = F, col = "#67017E") +
  theme(legend.position = "none") +
  #geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "Ln Seed length (mm)", y = "") 
ggsave("figures/fr_len-seed_len.eps", width = 8.9, height = 8.9, units = "cm")

ggarrange(pp1, pp3,
          ncol = 2, nrow = 1,
          labels = c("a", "b"),
          legend = "top")
ggsave("figures/fr_len-seed_len.png", 
       width = 20, height = 11, units = "cm", dpi = 600)

# fruit vs seed scatter with regression info per category
ggscatter(tree_data$dat, x = "Prom.Largo.semilla..mm.", y = "LARGO_FRUTO.cm.", 
          color = "disp_3cat", add = "reg.line", palette = "Set1",
          xlab = "Seed length (mm)", ylab = "Fruit length (cm)", 
          shape = "disp_3cat", fullrange = F) +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), 
        color = disp_3cat)) + 
  geom_abline(slope = 1, intercept = 0, lty = 2)

# Figure 3 Supp. fruit length ~ leaf length all dispersal system ####
# log-transf datavars
tree_data$dat %>% 
  ggplot(aes(x = (lv_len), y = (fr_len), 
             #shape = Sistema_de_Dispersion, 
             color = Sistema_de_Dispersion)) +
  geom_point() + 
  geom_smooth(method = "lm", 
              se = T, fullrange = F, alpha = .2) + 
  theme_classic(base_size = 12) +
  theme(legend.position = "top") + 
  labs(x = "Ln Leaf length (cm)", y = "Ln Fruit length (cm)",
       color = "Dispersal system", shape = "Dispersal system") +
  scale_color_manual(values = wes_palette("FantasticFox1",  7, type = "continuous"))
ggsave("figures/fr_lf_len_alldisp.png", 
       width = 20, height = 20, units = "cm", dpi = 600)
# seed area/roundness ####
# log seed area by dispersal system
tree_data$dat %>% 
  mutate(largo = round(Prom.Largo.semilla..mm., 2), 
         ancho = round(Prom.Ancho.sem..mm., 2)) %>% #select(s_l_t:cond) %>% tail()
  mutate(seed_area = ifelse(largo==ancho, 
                         ((largo/2))^2*pi, (largo/2)*(ancho/2)*pi)) %>% 
  #select(FAMILIA, Prom.Largo.semilla..mm., Prom.Ancho.sem..mm., sem_red) %>% 
  #filter(!is.na(sem_red)) %>% 
  ggplot(aes(x = disp_3cat, y = log1p(seed_area))) + 
  geom_boxplot(notch = T) +
  theme_classic()

# seed area vs. leaf length 
tree_data$dat %>% 
  mutate(largo = round(Prom.Largo.semilla..mm., 2), 
         ancho = round(Prom.Ancho.sem..mm., 2)) %>% #select(s_l_t:cond) %>% tail()
  mutate(seed_area = ifelse(largo==ancho, 
                            ((largo/2))^2*pi, (largo/2)*(ancho/2)*pi)) %>% 
  #select(FAMILIA, Prom.Largo.semilla..mm., Prom.Ancho.sem..mm., sem_red) %>% 
  #filter(!is.na(sem_red)) %>% 
  ggplot(aes(x = (lv_len), y = log1p(seed_area))) + 
  geom_point() + theme_classic() +
  labs(x = "Leaf length (cm)", y = "Seed area (mm2)") +
  geom_smooth(method = "lm")
  
