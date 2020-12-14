library(readxl)

PKD1 <- read_excel("~/Desktop/Papers-thesis/Data analysis/PCA_analysis_PKD1.xlsx")
PKD2 <- read_excel("~/Desktop/Papers-thesis/Data analysis/PCA_analysis_PKD2.xlsx")
PKD3 <- read_excel("~/Desktop/Papers-thesis/Data analysis/PCA_analysis_PKD3.xlsx")

PKD1_PCA <- data.matrix(PKD1)
PKD2_PCA <- data.matrix(PKD2)
PKD3_PCA <- data.matrix(PKD3)

# Create custom color palette (yellow to red gradient)
grad <- colorRampPalette(c("yellow","red"))

# Get a vector of colors and deal with negative values in PCA columns
colors <- grad(length(min(PKD1_PCA[,1]):max(PKD1_PCA[,1])))
index <- PKD1_PCA[,1]

plot(PKD1_PCA[,-1], col=colors[index], pch=19, posx = c(0.99, 1), title= "PKD1")

library(shape)
colorlegend(col=colors, zlim=range(PKD1_PCA[,1]), main = "time step", posy=c(0.2,0.9))

colors <- grad(length(min(PKD2_PCA[,1]):max(PKD2_PCA[,1])))
plot(PKD2_PCA[,-1], col=colors[index], pch=19, posx = c(0.99, 1), title= "PKD2")
colorlegend(col=colors, zlim=range(PKD2_PCA[,1]), main = "time step", posy=c(0.2,0.9))


colors <- grad(length(min(PKD3_PCA[,1]):max(PKD3_PCA[,1])))
plot(PKD3_PCA[,-1], col=colors[index], pch=19, posx = c(0.99, 1), title= "PKD3")
colorlegend(col=colors, zlim=range(PKD3_PCA[,1]), main = "time step", posy=c(0.2,0.9))




