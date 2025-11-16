###############################################################################################
## PHYLOGENETIC BIAS ANALISIS                                                                ##
## Script to assess phylogenetic clustering in an assemblage in relation to a species pool   ##
##                                                                                           ##
##  Pedro Abellan                                                                            ##
##  December 2014                                                                            ##
###############################################################################################

###############################################################################################
## FUNCTIONS FOR PHYLOGENETIC BIAS ANALYSIS
###############################################################################################

## Load libraries
library(ape)
library(picante)

## function to compute random MPD given a subsample of n.obs species 
mpd.random <- function(comm, dists, n.obs){
  ids <- sample(names(comm), n.obs)
  mpd.r <- mean(as.dist(dists[ids, ids]))
  return(mpd.r)
}

## function to compute random MNTD given a subsample of n.obs species
mntd.random <- function(comm, dists, n.obs){
  ids <- sample(names(comm), n.obs)
  dis <- dists[ids, ids]
  dis[which(dis==0)] <- NA
  mntd.r <- mean(apply(dis, 1, FUN=min, na.rm=T))
  return(mntd.r)
}

## function to compute mpd of an assemblage for a tree
mpd.tree <- function(tree.i, comm, runs ){
  id.obs <- colnames(comm)[which(comm[1,]==1)]
  n.obs <- length(id.obs)
  samp <- as.character(row.names(comm[1,]))
  dists.i = cophenetic(tree.i)      # the results from a call to cophenetic() on a phylo object

  ## calculating mntd.obs
  mpd.obs <- mean(as.dist(dists.i[id.obs, id.obs]))

  ## randomization test
  null.random.mpd <- replicate(runs, mpd.random(comm, dists.i, n.obs))

  ma.mpd<-c(mpd.obs,null.random.mpd) # add mpd.obs to null distribution
  r1<-rank(ma.mpd)[1] # rank of mpd.obs in null distribution
  p.value <- r1/(runs+1) # compute  one-tailed p-value
  test.mpd <- data.frame(sr=n.obs, mpd.obs, mean.mpd.rand=mean(null.random.mpd), rank=r1, p.value, runs, row.names=samp)
  return(test.mpd)
}

## function to compute mntd of an assemblage for a tree
mntd.tree <- function(tree.i, comm, runs ){

  id.obs <- colnames(comm)[which(comm[1,]==1)]
  n.obs <- length(id.obs)
  samp <- as.character(row.names(comm[1,]))
  dists.i = cophenetic(tree.i)      # the results from a call to cophenetic() on a phylo object
  dis <- dists.i[id.obs, id.obs]
  dis[which(dis==0)] <- NA

  ## calculating mntd.obs
  mntd.obs <- mean(apply(dis, 1, FUN=min, na.rm=T))

  ## randomization test
  null.random.mntd <- replicate(runs, mntd.random(comm, dists.i, n.obs))

  ma.mntd<-c(mntd.obs,null.random.mntd) # add mntd.obs to null distribution
  r1<-rank(ma.mntd)[1] # rank of mntd.obs in null distribution
  p.value <- r1/(runs+1) # compute  one-tailed p-value
  test.mntd <- data.frame(sr=n.obs, mntd.obs, mean.mntd.rand=mean(null.random.mntd), rank=r1, p.value, runs, row.names=samp)
  return(test.mntd)
}


###############################################################################################
## ANALYSIS OF PHYLOGENETIC BIAS FOR ONE TREE
###############################################################################################

## Load data
tree0 <- read.nexus("tree.nex")  # load file with one tree
comm0 <- read.table("comm.txt", h=T, row.names=1) # load file with presence/absence. Species in columns and asemblages in rows

dat <- match.phylo.comm(tree0, comm0)  # match assemblage data and phylogenies
tree <- dat$phy
comm <- dat$comm[2,]  # select relevant assemblage for analyses. In my case, the second row

mpd.comm <- mpd.tree(tree,comm,999)
mntd.comm <- mntd.tree(tree,comm,999)


###############################################################################################
## ANALYSIS OF PHYLOGENETIC BIAS FOR MULTIPLE TREES
###############################################################################################

## Load data
trees <- read.nexus("trees.nex")  # load file with multiple trees
comm0 <- read.table("comm.txt", h=T, row.names=1) # load file with presence/absence. Species in columns and asemblages in rows

dat <- match.phylo.comm(trees[[1]], comm0)  # match assemblage data and phylogenies (uses just the first tree)
comm <- dat$comm[2,]  # select relevant assemblage for analyses. In my case, the second row

## MPD calculation for multiple trees
out.mpd.0 <- lapply(FUN=mpd.tree, trees, com, 999)
out.mpd <- do.call("rbind", out.mpd.0)
write.table(out.mpd, "results_mpd.txt", row.names=F)

## MNTD calculation for multiple trees
out.mntd.0 <- lapply(FUN=mntd.tree, trees, com, 999)
out.mntd <- do.call("rbind", out.mntd.0)
write.table(out.mntd, "results_mntd.txt", row.names=F)












































