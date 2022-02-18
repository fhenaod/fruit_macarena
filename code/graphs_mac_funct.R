library(tidyverse)
library(ggpubr)
library(wesanderson)

fruit_dat_raw <- read.csv("data/Base plantas Tinigua modificable_Agosto.csv", header = T )
fruit_dat_raw %>% #names() %>% 
  select(where(is.numeric)) %>% names()

fruit_dat_raw %>% select(
                         LARGO_FRUTO.cm., Fr.width..cm.,
                         Largo.hojas.total.cm., Largo.lamina.cm., 
                         Area.total.cm2., Area.lámina.cm2.,
                         Largo.peciolo.cm.,
                         Prom.Largo.semilla..mm., Prom.Ancho.sem..mm., Redondez..L.A.
                         , Fruit.Dry.Weight, Semillas.Fruto,
                         DBH.2017, Densidad..specific.gravity.,
                         ) %>% 
  gather("variable", "value", -LARGO_FRUTO.cm.) %>% 
  ggplot(aes(x = (value), y = log1p(LARGO_FRUTO.cm.))) + geom_point() + 
  geom_smooth(method = "lm", col = "red", se = T, fullrange = F) + 
  #theme_minimal() + 
  ggthemes::theme_clean() + theme(panel.grid.major = element_blank()) +
  facet_wrap(~variable, scales = "free")

# Figure 1, panel with leaf length and area by dispersal system ####
# fruit length ~ leaf length 
pp1 <- tree_data$dat %>% 
  ggplot(aes(x = log(Largo.hojas.total.cm.), 
             y = log(LARGO_FRUTO.cm.))) + 
  geom_point(cex = .8) + theme_classic(base_size = 12) + 
  geom_smooth(method = "lm", se = T, show.legend = F, col = "#67017E") +
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

# fruit lenght ~ seed lenght
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

# isometry ####
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

# seed data transpostition 
tree_data$dat <- tree_data$dat %>% 
  mutate(sem_len_new = ifelse(sem_len > sem_wd, sem_len, sem_wd),
         sem_wd_new = ifelse(sem_wd < sem_len, sem_wd, sem_len)) 

ggscatter(tree_data$dat, x = "sem_len_new", y = "sem_wd_new", 
          color = "disp_3cat", add = "reg.line", palette = "Set1",
          xlab = "Seed length (mm)", ylab = "Seed width (mm)", 
          shape = "disp_3cat", fullrange = F) +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), 
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

# Figure 2, y ~ logp1(x) by dispersal system ####
tree_data$dat %>% 
  ggplot(aes(x = exp(sem_len_new), y = exp(sem_wd_new), 
             color = disp_3cat, shape = disp_3cat)) + 
  geom_point() + theme_classic(base_size = 12) + 
  geom_smooth(formula = y ~ log1p(x), se = T, show.legend = F) +
  theme(legend.position = "top") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "Seed length (mm)", y = "Seed width (mm)", 
       color = "Dispersal system", shape = "Dispersal system") +
  scale_color_manual(values = wes_palette("FantasticFox1",  3, type = "continuous"))
  #scale_color_viridis_d(option = "magma", direction = 1, begin = 0, end = .65) 
#geom_abline(slope = .97, intercept = -0.49, lwd = 1, col = "#488f31") +
#geom_abline(slope = .86, intercept = -.76, lwd = 1, col = "#e08861") 
ggsave("figures/seed_isom_3cat.png", 
       width = 15, height = 15, units = "cm", dpi = 600)

# regression formula per category
paste0("Y_endo = ", round(coefficients(seed_isoxdisp)[1], 2), " + ", round(coefficients(seed_isoxdisp)[2],2), "x")
paste0("y_no = ", round(coefficients(seed_isoxdisp)[1]+coefficients(seed_isoxdisp)[3], 2), " + ", round(coefficients(seed_isoxdisp)[2]+coefficients(seed_isoxdisp)[4],2), "x")

# with density distributions
p <- tree_data$dat %>% 
  ggplot(aes(x = exp(sem_len_new), 
             y = exp(sem_wd_new), 
             color = disp_3cat, shape = disp_3cat)) +
  geom_point() + theme_classic() +
  theme(legend.position = "top") +
  labs(x = "Seed length (mm)", y = "Seed width (mm)") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_smooth(formula = (y) ~ log1p(x))
ggExtra::ggMarginal(p, type = "density", 
                    groupColour = T, groupFill = T)

# Figure 3 average fruit and laminar lenght ####
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

# average leave lenght by family
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

# Figure S1, fruit length ~ leaf length by dispersal system ####
# log-transf datavars
tree_data$dat %>% ggplot(aes(x = lv_len, y = fr_len, 
                             shape = disp_3cat, color = disp_3cat)) +
  geom_point() + geom_smooth(method = "lm", se = T, fullrange = F) + 
  theme_classic(base_size = 12) +
  theme(legend.position = "top") + 
  labs(x = "Ln Leaf length (cm)", y = "Ln Fruit length (cm)",
       color = "Dispersal system", shape = "Dispersal system") +
  scale_color_manual(values = wes_palette("FantasticFox1",  3, type = "continuous"))
ggsave("figures/fr_lf_len_disp.png", 
       width = 15, height = 15, units = "cm", dpi = 600)

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
  



