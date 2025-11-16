###############################################################################################
## IDENTIFYING WHICH CLADES SIGNIFICANTLY CONTRIBUTE TO THE PHYLOGENETIC STRUCTURE
## For each node in the phylogeny, it is tested whether it has significantly
## more descendent taxa in a sample than would be expected by chance by
## means of a randomization test (equivalent to NODESIG function in Phylocom software)
##
##  Pedro Abellan
##  December 2014
##
###############################################################################################

library(geiger)
library(picante)

## Load data
tree0 <- read.nexus("tree.nex")  # load file with phylogenetic tree
comm0 <- read.table("comm.txt", h=T, row.names=1) # load file with presence/absence. Species in columns and asemblages in rows

dat <- match.phylo.comm(tree0, comm0)  # match assemblage data and phylogenies
phy <- dat$phy
comm <- dat$comm[2,]  # select relevant assemblage for analyses. In my case, the second row

n.internal.nodes <- phy$Nnode # number of nodes in tree
n.tips <- length(phy$tip.label) # number of tips in tree
spp.pool <- names(comm) # get species names in species pool
spp.obs <- colnames(comm)[which(comm==1)] # taxa in sample
sr.sample <- length(spp.obs) # number of species in sample

## funtion to compute random sr in a node
sr.node.rand <- function(spp.pool, sr.sample, spp.node){
spp.rand <- sample(spp.pool, sr.sample) # sample a equivalent number of spp from total spp pool
sr.node.rand  <- length(spp.rand[spp.rand%in%spp.node])
return(sr.node.rand)
}

out <- data.frame() # prepare dataframe to store reults
## calculating mpdd.obs for node
for(i in (n.tips+1):(n.tips+n.internal.nodes)){
  node <- i
  spp.node <- tips(phy, node) # taxa in node
  sr.node.obs <- length(spp.obs[spp.obs%in%spp.node])  # number taxa of sample in node

  ## randomization test for node
  runs <- 999
  null.sr.node.rand <- replicate(runs, sr.node.rand(spp.pool, sr.sample, spp.node))  # null distribution of sr in node
  ma.sr.node <- c(sr.node.obs,null.sr.node.rand) # add sr.node.obs to null distribution
  r1<-rank(ma.sr.node)[1] # rank of mpd.obs in null distribution
  p.value <- 1-(r1/(runs+1)) # compute  one-tailed p-value
  test.sr.node <- data.frame(sr.node=length(spp.node),sr.node.obs, mean.sr.rand=mean(null.sr.node.rand), rank=r1, p.value, runs, row.names=node)
  out <- rbind(out, test.sr.node)
}

out2 <- out[which(out$p.value<0.05),]
write.table(out, "out.nodesig.txt")
write.table(out2, "out.nodesig_sig.txt")

## Export significant clades with species in each clade
sig.clades <- as.numeric(row.names(out2))
sink("clades.nodesig.txt")
for (i in sig.clades){
clades.i <- as.vector(tips(phy, i))
print(i)
print(clades.i)
}
sink()


## Identify significant clades on tree
pdf("tree.nodesig.pdf", width=20, height=200)
reproLabel <- character(length(phy$tip.label))
names(reproLabel) <- names(comm)
reproLabel[comm==0] <- NA
reproLabel[comm==1] <- "red"
plot(phy, cex=0.9, label.offset=0.2, tip.color= c(reproLabel[match(phy$tip.label, names(reproLabel))]), lwd=0.02)
add.scale.bar(length=10, x=0.7, y=-8)
nodelabels("*", sig.clades, cex=2, frame="none", col = "red")
dev.off()


