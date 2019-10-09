library(treeplyr)
tree <- read.tree("data/phylo_spp_tinigua_noded_nwk.txt")
tree <- multi2di(tree)
is.binary(tree)
class<-tree$node.label[grep("sida",tree$node.label)]
orders<-tree$node.label[grep("ales",tree$node.label)]

n_node <- function(tree, taxa){
  n_node<-c()
  for (i in 1:length(taxa)){
    n_node[i] <- Ntip(tree)+which(tree$node.label==paste0(taxa[i]))
  }
  n_node
}
n_class_node <- n_node(tree, class)
n_ord_node <- n_node(tree, orders)

dt<-data.frame(taxon = c(class, orders),node = c(n_class_node, n_ord_node) ) 

node <- c(dt[which(dt$taxon=="Poales"),2], dt[which(dt$taxon=="Zingiberales"),2], dt[which(dt$taxon=="Ericales"),2], 
          dt[which(dt$taxon=="Apiales"),2], dt[which(dt$taxon=="Magnoliales"),2], dt[which(dt$taxon=="Piperales"),2], 
          dt[which(dt$taxon=="Commelinales"),2], dt[which(dt$taxon=="Myrtales"),2], dt[which(dt$taxon=="Cucurbitales"),2],
          dt[which(dt$taxon=="Asparagales"),2], dt[which(dt$taxon=="Gentianales"),2], dt[which(dt$taxon=="Malpighiales"),2],
          dt[which(dt$taxon=="Celastrales"),2], dt[which(dt$taxon=="Dioscoreales"),2], dt[which(dt$taxon=="Sapindales"),2],
          dt[which(dt$taxon=="Ranunculales"),2], dt[which(dt$taxon=="Lamiales"),2], dt[which(dt$taxon=="Liliopsida"),2],
          dt[which(dt$taxon=="Rosidae"),2], Ntip(tree)+which(tree$node.label==paste0("Spermatophyta")))

age.min <- c(101,  30, 100, 44, 108, 111,  80,  68,  89,  97, 58, 92,   81, 101,  98, 103,  68, 123, 105, 289)
age.max <- c(115, 109, 118, 85, 125, 122, 110, 110, 103, 113, 80, 92, 81.2, 124, 125, 122, 103, 142, 115, 337)
soft.bounds <- rep(FALSE, length(node))
cal_ages <- data.frame(node, age.min, age.max, soft.bounds)

plot(tree, show.tip.label = F, show.node.label = F, no.margin = T, type = "fan")
nodelabels(tree$node.label[node-Ntip(tree)], node = node, frame = "none", col = "red", cex=.7, adj = c(1,1.2))
nodelabels(node = node, frame = "none", col = "blue", cex=.5, adj = c(3,3))

tree_calib <- chronos(tree, calibration = cal_ages)
is.ultrametric(tree_calib)
write.tree(tree_calib,"data/tree_calib.nwk")
