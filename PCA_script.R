setwd("~/Desktop/Papers-thesis/Data analysis/")

library(bio3d)
#demo("md")

#trajectory file in dcd format instead of xtc
dcdfile <- "PKD2_interface.dcd"
#starting PDB structure
pdbfile <- "firstframe_PKD2.pdb"

#For_own_analysis_remove_system_file
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
plot(rd, typ = "l", ylab = "RMSD", xlab = "Frame No.")
points(lowess(rd), typ = "l", col = "red", lty = 2, lwd = 2)
hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram", xlab="RMSD")
lines(density(rd), col="gray", lwd=3)
summary(rd)

# Root Mean Squared Fluctuations (RMSF)
rf <- rmsf(xyz[, ca.inds$xyz])
plot(rf, ylab = "RMSF", xlab = "Atom Position", typ="l")

rf2 <- rmsf(xyz2[, ca.inds2$xyz])
plot(rf2, ylab = "RMSF", xlab = "Residue Position", typ="l")

# Principal Component Analysis
pc <- pca.xyz(xyz[, ca.inds$xyz])
plot(pc, col = bwr.colors(nrow(xyz)))

pc2 <- pca.xyz(xyz[, ca.inds2$xyz])
plot(pc2, col = bwr.colors(nrow(xyz)))

# Cluster in PC space
hc <- hclust(dist(pc$z[, 1:7]), method = "average")
hc2 <- hclust(dist(pc2$z[, 1:7]), method = "average")

#average can be used as UPGMA

#K here can be decided from visualization of possible numbers of clusters from PCs
grps <- cutree(hc, k = 4)
plot(pc, col = grps)

grps2 <- cutree(hc2, k = 4)
plot(pc2, col = grps2)

# check which residues produces variability among PCs

plot.bio3d(pc$au[,1], ylab="PC1 (A)", xlab="Residue Position", typ="l")
points(pc$au[,2], typ="l", col="blue")

#visualise those in VMD with tube representation to see which residues contribute to variation more
p1 <- mktrj.pca(pc, pc=1, b=pc$au[,1], file="pc1.pdb")
p2 <- mktrj.pca(pc, pc=2,b=pc$au[,2], file="pc2.pdb")

#Color scale from blue to red depict low to high atomic displacements.

# Cross-Correlation Analysis
cij <- dccm(xyz[, ca.inds$xyz])

plot(cij)
pymol(cij, pdb, type="launch",exefile='/Applications/PyMOL.app/Contents/bin/pymol')


