library(treeplyr)
tree <- read.tree("data/phylo_spp_tinigua_noded_nwk.txt")
tree <- multi2di(tree)
is.binary(tree)
class<-tree$node.label[grep("sida",tree$node.label)]
orders<-tree$node.label[grep("ales",tree$node.label)]
family<-tree$node.label[grep("ceae",tree$node.label)]
no_rank <- c(tree$node.label[grep("permae",tree$node.label)], "Spermatophyta")

n_node <- function(tree, taxa){
  n_node<-c()
  for (i in 1:length(taxa)){
    n_node[i] <- Ntip(tree)+which(tree$node.label==paste0(taxa[i]))
  }
  n_node
}
n_class_node <- n_node(tree, class)
n_ord_node <- n_node(tree, orders)
n_fam_node <- n_node(tree, family)

dt <- data.frame(taxon = c(class, orders, no_rank),node = c(n_class_node, n_ord_node, n_no_rank_node)) 

taxon <- c("Poales", "Zingiberales", "Ericales", "Apiales", "Magnoliales", "Piperales", 
          "Commelinales", "Myrtales", "Cucurbitales", "Asparagales", "Gentianales", 
          "Malpighiales", "Celastrales", "Dioscoreales", "Sapindales", "Ranunculales", 
          "Lamiales", "Liliopsida", "Rosidae", "Spermatophyta", "Rosales", "Aquifoliales", "Caryophyllales", "Fabales")
age.min <- c(101,  30, 100, 44, 108, 111,  80,  68,  89,  97, 58, 92,   81, 101,  98, 103,  68, 123, 105, 289, 42,   68, 78, 75)
age.max <- c(115, 109, 118, 85, 125, 122, 110, 110, 103, 113, 80, 92, 81.2, 124, 125, 122, 103, 142, 115, 337, 161, 100, 92, 87)
soft.bounds <- rep(FALSE, length(node))
cal_ages <- data.frame(node, age.min, age.max, soft.bounds)

plot(tree, show.tip.label = F, show.node.label = F, no.margin = T, type = "fan")
nodelabels(tree$node.label[node-Ntip(tree)], node = node, frame = "none", col = "red", cex=.7, adj = c(1,1.2))
nodelabels(node = node, frame = "none", col = "blue", cex=.5, adj = c(3,3))

tree_calib <- chronos(tree, calibration = cal_ages)
is.ultrametric(tree_calib)
write.tree(tree_calib,"data/tree_calib.nwk")
