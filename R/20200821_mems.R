# Figure out how to get MEMs for this

library(ade4)
library(adespatial)
library(adegraphics)
library(spdep)
library(maptools)
library(vegan)
library(dplyr)

# This is my resource, up to section 5.
# https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html


# Starting with one of the communities from Emlyn\
all_data <- read.table(file = "old_tests/Resetarits&al2018FinalData.txt")
dispersal_data <- read.table(file = "old_tests/dispersal_data.txt")

protists_corrected <- read.table("old_tests/dispersal_data.txt", header=T)
head(protists_corrected)

all_data %>%
  dplyr::filter(type == "E") -> no_removal_communities

# Control metacommunities spatial structure -----

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

# return Euclidean distances for the coordinates based on the spatial neighborhood
grid.thresh.d1 <- nbdists(nb1, coords)

# Spatial weights are defines as a function of distance
fdist <- lapply(grid.thresh.d1, function(x) 1-x/max(dist(coords)))

# This is the spatial weighting matrix
grid.invdist.lw <- nb2listw(nb1, glist = fdist, style="B")
print(listw2mat(grid.invdist.lw)[1:6,1:6], digits=2)

# Worried that we are taking these euclidean distances and coordinates and such, as real distances, but they are based on a sampling scheme, not actual distances?

# MEMs ----

# THIS is the function or step that is creating the MEM based on the distance matrices that come from the spdep::listw objects.
grid.invdist.MEM <- scores.listw(grid.invdist.lw, store.listw = TRUE)
# could also use mem()
grid.invdist.MEM <- mem(grid.invdist.lw, store.listw = TRUE)

grid.invdist.MEM
summary(grid.invdist.MEM)
attr(grid.invdist.MEM, which = "value")
barplot(attr(grid.invdist.MEM, which = "value"), main="Eigenvalues of the spatial weighting matrix")

s.value(coords, grid.invdist.MEM[, c(1,3,7,10)])

# Compute moran's I
grid.MEM.Moran <- moran.randtest(grid.invdist.MEM, listw = grid.invdist.lw, nrepet = 999)
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

grid.invdist.MEM.pos.sig <- sig.MEM

# Those are the significant MEMs for the control metacommunities, the ones with no removals

