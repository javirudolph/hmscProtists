# Get the MEMs for control metacommunities
# These metacommunities consist of 36 patches and are all connected through a dendritic network

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

resetarits_data <- read.table(here::here("data", "Resetarits&al2018FinalData.txt"))
head(resetarits_data)

library(ade4)
library(adespatial)
library(adegraphics)
library(spdep)
library(maptools)
library(vegan)

x <- c(1,1,1,1,1,1,
       2,2,2,1.999,2,2,
       3,3,3,2.998,3,3,
       4,4,4,3.997,4,4,
       5,5,5,4.996,5,5,
       6,6,6,5.995,6,6)
y <- c(5.996,4.997,3.998,2.999,2,1.001,
       5.996,4.997,3.998,2.999,2,1.001,
       5.996,4.997,3.998,2.999,2,1.001,
       5.996,4.997,3.998,2.999,2,1.001,
       5.996,4.997,3.998,2.999,2,1.001,
       5.996,4.997,3.998,2.999,2,1.001)

coords <- coordinates(cbind(x,y))
plot(coords)


# Build spatial neighborhood and using 0.9995 as the threshold
nb1 <- spdep::dnearneigh(coords, 0, 0.9995, longlat=F)
adegraphics::s.label(coords, nb = nb1)

# Define spatial weighing matrices
neighb2.lw <- spdep::nb2listw(nb1, style="B")

# Now as a matrix of 1/0
weight.mat<-spdep::listw2mat(neighb2.lw)

grid.invdist.MEM <- mem(neighb2.lw, store.listw = TRUE)

grid.invdist.MEM
summary(grid.invdist.MEM)
attr(grid.invdist.MEM, which = "value")
barplot(attr(grid.invdist.MEM, which = "value"), main="Eigenvalues of the spatial weighting matrix")

s.value(coords, grid.invdist.MEM[, c(1,3,7,10)])

# Compute moran's I
grid.MEM.Moran <- moran.randtest(grid.invdist.MEM, listw = neighb2.lw, nrepet = 999)
grid.MEM.Moran

head(attr(grid.invdist.MEM, "values") / grid.MEM.Moran$obs)

# Which MEMs are significant?

which(grid.MEM.Moran$pvalue <= 0.05) # MEM with significant spatial correlation
length(which(grid.MEM.Moran$pvalue <= 0.05))

# Now, we want to select those significant ones but only keep the eigenvectors that have positive eigenvalues.

# These are the significant MEMs
sig.MEM <- grid.invdist.MEM[, c(which(grid.MEM.Moran$pvalue <= 0.05))]

# Now, choose the ones with positive eigenvalues
attr(sig.MEM, which = "values") # Seems like they are all positive

# These are the significant MEMs for the control metacommunities
ctrl.mem <- sig.MEM

saveRDS(ctrl.mem, file = "paper/ctrl_mem.RDS")

