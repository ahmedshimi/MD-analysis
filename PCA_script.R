setwd("~/Desktop/Papers-thesis/Data analysis/PCA analysis/")

library(bio3d)
#demo("md")

#trajectory file in dcd format instead of xtc
dcdfile <- "PKD2_interface_150.dcd"
#starting PDB structure
pdbfile <- "firstframe_PKD2.pdb"
dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)

dim(dcd)
length(pdb$xyz)
ncol(dcd) == length(pdb$xyz)

# Trajectory Frame Superposition on all protein atoms
ca.inds <- atom.select(pdb, "protein")
xyz <- fit.xyz(fixed = pdb$xyz, mobile = dcd, 
                  fixed.inds = ca.inds$xyz, 
                  mobile.inds = ca.inds$xyz)

# Trajectory Frame Superposition on Calpha atoms (Residues)
ca.inds2 <- atom.select(pdb,elety = "CA")
xyz2 <- fit.xyz(fixed = pdb$xyz, mobile = dcd, 
                fixed.inds = ca.inds2$xyz, 
                mobile.inds = ca.inds2$xyz)

# Root Mean Square Deviation (RMSD)

rd <- rmsd(xyz[1, ca.inds2$xyz], xyz[, ca.inds2$xyz])
plot(rd, typ = "l", ylab = "RMSD", xlab = "Frame No.", main="PKD2 RMSD")
points(lowess(rd), typ = "l", col = "red", lty = 2, lwd = 2)
hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram", xlab="RMSD")
lines(density(rd), col="gray", lwd=3)
summary(rd)

# Root Mean Squared Fluctuations (RMSF)
rf <- rmsf(xyz[, ca.inds$xyz])
plot(rf, ylab = "RMSF", xlab = "Atom Position", typ="l")

rf2 <- rmsf(xyz2[, ca.inds2$xyz])
plot(rf2, ylab = "RMSF", xlab = "Residue Position", typ="l", main="PKD2 RMSF")

# Principal Component Analysis
set.seed(1)
pc <- pca.xyz(xyz[, ca.inds$xyz], use.svd = TRUE)
plot(pc, col = bwr.colors(nrow(xyz)))
library(shape)
colorlegend(col= bwr.colors(nrow(xyz)), zlim = range(1:1501), main="Time Step(ps)")

pc2 <- pca.xyz(xyz[, ca.inds2$xyz], use.svd = TRUE)
plot(pc2, col = bwr.colors(nrow(xyz)))

# Cluster in PC space
library("adegenet")
clust <- find.clusters(pc$z, n.pca = 100)
clust$Kstat[1:20]
graphics.off()
plot(clust$Kstat[1:30], type="b", col="blue", xlab = "", ylab = "")
title(xlab="Number of Clusters", ylab= "BIC", main = "PKD2 MD Analysis")

hc <- hclust(dist(pc$z[, 1:7]), method = "average")
hc2 <- hclust(dist(pc2$z[, 1:7]), method = "average")

#average can be used as UPGMA

#K-means clustering
#pkd2.pkcl <- kmeans(pc$z,5,7)
#plot(pc, col = pkd2.pkcl$cluster)

#K here can be decided from visualization of possible numbers of clusters from PCs
grps <- cutree(hc, k = 7)
table(grps)
plot(pc, col = grps)

cluster_center = aggregate(pc$z,list(cluster=grps),mean)
pc.z <- cbind(1:nrow(pc$z), pc$z)

cluster4  <- data.frame(pc.z[grps==4,])

index_closest_to_center <- function(n){
  cluster  <- data.frame(pc.z[grps==n,])
  names(cluster_center) <- names(cluster)
  dist =c(0, nrow(cluster))
  angle =c(0, nrow(cluster))
  for (i in 1:nrow(cluster)){
    dist[i] = dist(rbind(cluster[i,-1], cluster_center[n,-1]))
    angle[i] = sum(as.vector(cluster[i,-1])*as.vector(cluster_center[n,-1])) / ( sqrt(sum(as.vector(cluster[i,-1]) * as.vector(cluster[i,-1]))) * sqrt(sum(as.vector(cluster_center[n,-1]) * as.vector(cluster_center[n,-1]))))
  }
  cluster$dist <- dist
  cluster$angle <- angle
  closest_dist <- cluster[which.min(cluster$dist), ]$X1
  closest_angle <- cluster[which.max(cluster$angle), ]$X1
  angle_1 <- cluster[which.max(cluster$angle), ]$angle
  return(c(closest_dist, closest_angle, angle_1))
}

index_closest_to_center(1)
index_closest_to_center(2)
index_closest_to_center(3)
index_closest_to_center(4)
index_closest_to_center(5)
index_closest_to_center(6)
index_closest_to_center(7)


grps2 <- cutree(hc2, k = 4)
plot(pc2, col = grps2)

# check which residues produces variability among PCs

plot.bio3d(pc2$au[,1], ylab="PC1 (A)", xlab="Residue Position", typ="l")
points(pc2$au[,2], typ="l", col="blue")

#visualise those in VMD with tube representation to see which residues contribute to variation more
p1 <- mktrj.pca(pc2, pc=1, b=pc2$au[,1], file="pc1.pdb")
p2 <- mktrj.pca(pc2, pc=2,b=pc2$au[,2], file="pc2.pdb")

#Color scale from blue to red depict low to high atomic displacements.

# Cross-Correlation Analysis
cij <- dccm(xyz[, ca.inds2$xyz])

plot(cij, main= "PKD2 Residues Cross Correlation")
pymol(cij, pdb, type="launch",exefile='/Applications/PyMOL.app/Contents/bin/pymol')


