library(tidyverse)
library(vegan)

diet <- read.delim("DietaPrimatesOrinoq.csv", sep = ",")
diet_data <- diet %>% select(ESPECIES, Sapajus, Lagothrix, Ateles, Alouatta, Callicebus, Saimiri) %>% 
  replace(., is.na(.), 0) %>% t() %>% data.frame()
colnames(diet_data) <- diet_data[1,]
diet_data <- diet_data[-1, -395]
diet_data <- sapply(diet_data, as.numeric)
row.names(diet_data) <- c("Sapajus", "Lagothrix", "Ateles", "Alouatta", "Callicebus", "Saimiri")
  
nmds <- metaMDS(diet_data, distance = "bray", trymax=1000, binary = T, 
                autotransform=FALSE)
nmds1<-plot(nmds, type = "t", cex= .05)
nmds2<-ordiplot(nmds1$sites, type="t", cex=.7, cex.axis=1, 
                cex.lab=1.5, xlim = c(-0.6,1), ylim = c(-0.4,0.6), main = "Bray-Curtis")
nmds2<-ordiplot(nmds1$species, type="t", cex=.7,)
