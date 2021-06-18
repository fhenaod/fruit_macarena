library(V.PhyloMaker)
library(tidyverse)

df2tree <- fruit_dat_raw %>% select(species = ESPECIE, genus = Genus, family = FAMILIA)
df2tree$species <- sub("_", " ", df2tree$species)

head(df2tree)
dim(df2tree)

tree_s3 <- phylo.maker(sp.list = df2tree, tree = GBOTB.extended, 
            nodes = nodes.info.1, scenarios = "S3")

plot(tree_s3$scenario.3, #type = "fan", 
     show.tip.label = F, no.margin = T)
nodelabels(col = "red", frame = "none", cex = .5)

tree_s3$scenario.3 <- tree_s3$scenario.3 %>% multi2di()
tree_s3$scenario.3 %>% is.binary()
tree_s3$scenario.3 %>% is.ultrametric()

write.tree(tree_s3$scenario.3, file = "data/tree_mac_s3.tre")

plot(extract.clade(tree_s3$scenario.3, 1056), type = "fan", 
     show.tip.label = F, no.margin = T) # extrae Spermatophyta
