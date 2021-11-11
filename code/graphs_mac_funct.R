library(tidyverse)

fruit_dat_raw <- read.csv("data/Base plantas Tinigua modificable_Agosto.csv", header = T )
fruit_dat_raw %>% #names() %>% 
  select(where(is.numeric)) %>% names()

# fruit length ~ leaf length ####
fruit_dat_raw %>% ggplot(aes(x = log(Largo.hojas.total.cm.), y = log(LARGO_FRUTO.cm.))) +
  geom_point() + theme_classic(base_size = 10) + 
  geom_smooth(method = "lm", se = T, show.legend = F, col = "red") +
  theme(legend.position = "none") + 
  #geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "Ln Leaf length (cm)", y = "Ln Fruit length (cm)") +
  scale_color_viridis_d(option = "viridis", direction = -1, begin = 0, end = .65)
ggsave("fr_len-lv_len.eps", width = 8.9, height = 8.9, units = "cm")

# fruit length ~ leaf length by dispersal system
tree_data$dat %>% ggplot(aes(x = lv_len, y = fr_len, 
                             shape = disp_3cat, color = disp_3cat)) +
  geom_point() + geom_smooth(method = "lm", se = T, fullrange = F) + 
  theme_classic() +
  scale_color_viridis_d(begin = 0.2, end = .5) +
  theme() + 
  labs(x = "Ln Leaf length (cm)", y = "Ln Fruit length (cm)")

# panel with leaf length and area by dispersal system
pp1 <- tree_data$dat %>% 
  ggplot(aes(x = (lv_len), y = (fr_len), color = disp_3cat, shape = disp_3cat)) + 
  geom_point() + theme_classic(base_size = 12) + 
  geom_smooth(method = "lm", se = T, show.legend = F) +
  theme(legend.position = "none") +
  #geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "Ln Leaf length (cm)", y = "Ln Fruit length (cm)", 
       color = "Dispersal system", shape = "Dispersal system") +
  scale_color_viridis_d(option = "viridis", direction = -1, begin = 0, end = .65)

pp2 <- tree_data$dat %>% 
  ggplot(aes(x = (ar_tot), y = (fr_len), color = disp_3cat, shape = disp_3cat)) + 
  geom_point() + theme_classic(base_size = 12) + 
  geom_smooth(method = "lm", se = T, show.legend = F) +
  theme(legend.position = "none") +
  #geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "Ln Leaf area (mm2)", y = "Ln Fruit length (mm)", 
       color = "Dispersal system", shape = "Dispersal system") +
  scale_color_viridis_d(option = "viridis", direction = -1, begin = 0, end = .65)

ggarrange(pp1, pp2,
          ncol = 2, nrow = 1,
          labels = c("A", "B"),
          legend = "top")
ggsave("phylo_lm_leav_seed.eps", width = 18.4, height = 8.9, units = "cm")

# fruit lenght ~ seed lenght
fruit_dat_raw %>% ggplot(aes(x = log(Prom.Largo.semilla..mm.), y = log(LARGO_FRUTO.cm.))) +
  geom_point() + theme_classic(base_size = 10) + 
  geom_smooth(method = "lm", se = T, show.legend = F, col = "red") +
  theme(legend.position = "none") +
  #geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "Ln Seed length (mm)", y = "Ln Fruit length (cm)") +
  scale_color_viridis_d(option = "viridis", direction = -1, begin = 0, end = .65)
ggsave("fr_len-seed_len.eps", width = 8.9, height = 8.9, units = "cm")


# isometry ####
# seed width ~ seed length
fruit_dat_raw %>% ggplot(aes(x = log(Prom.Largo.semilla..mm.), y = log(Prom.Ancho.sem..mm.))) +
  geom_point() + theme_classic(base_size = 10) + 
  geom_smooth(method = "lm", se = T, show.legend = F, col = "red") +
  theme(legend.position = "none") +
  #geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "Ln Mean Seed length (mm)", y = "Ln Mean Seed width (mm)") +
  scale_color_viridis_d(option = "viridis", direction = -1, begin = 0, end = .65)
ggsave("seed_len-seed_widt.eps", width = 8.9, height = 8.9, units = "cm")

# scatter for log-log with regression info
ggscatter(tree_data$dat, x = "sem_len", y = "sem_wd", 
          color = "disp_3cat", add = "reg.line", palette = "Set1",
          xlab = "Ln Seed length (mm)", ylab = "Ln Seed width (mm)", 
          shape = "disp_lato", fullrange = F) +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), 
        color = disp_3cat)) + 
  geom_abline(slope = 1, intercept = 0, lty = 2)

# color by dispersal system
tree_data$dat %>% 
  ggplot(aes(x = exp(sem_len_new), y = exp(sem_wd_new), 
             color = disp_3cat, shape = disp_3cat)) + 
  geom_point() + theme_classic(base_size = 12) + 
  geom_smooth(formula = y ~ log1p(x), se = T, show.legend = F) +
  theme(legend.position = "top") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "Seed length (mm)", y = "Seed width (mm)", 
       color = "Dispersal system", shape = "Dispersal system") +
  #scale_color_manual(values = c("#488f31", "#e08861")) +
  scale_color_viridis_d(option = "magma", direction = 1, begin = 0, end = .65) 
#geom_abline(slope = .97, intercept = -0.49, lwd = 1, col = "#488f31") +
#geom_abline(slope = .86, intercept = -.76, lwd = 1, col = "#e08861") 
ggsave("seed_isom_3cat_.eps", width = 8.9, height = 8.9, units = "cm")

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

