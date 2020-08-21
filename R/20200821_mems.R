# Figure out how to get MEMs for this

library(ade4)
library(adespatial)
library(adegraphics)
library(spdep)
library(maptools)
library(vegan)
library(dplyr)

# This is my resource:
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


grid.inv.dist <- lapply(grid.thresh.d1, function(x) 1-x/max(dist(coords)))#item 9.3.2 from Bivand's book:
grid.invdist.lw <- nb2listw(nb1, glist=grid.inv.dist, style="B")
print(listw2mat(grid.invdist.lw)[1:6,1:6], digits=2)

# THIS is the function or step that is creating the MEM based on the distance matrices that come from the spdep::listw objects.
grid.invdist.MEM <- scores.listw(grid.invdist.lw, store.listw = TRUE)

summary(grid.invdist.MEM)
attr(grid.invdist.MEM, which = "value")
barplot(attr(grid.invdist.MEM, which = "value"), main="Barplot of MEM eigenvalues")


# grid.MEM.Moran <- test.scores(grid.invdist.MEM, grid.invdist.lw, 999)
# The function `test.scores()` is deprecated now. To get significance try using:
grid.MEM.Moran <- moran.randtest(grid.invdist.MEM, listw = grid.invdist.lw, nrepet = 999)

# Which MEMs are significant?


which(grid.MEM.Moran$pvalue <= 0.05) # MEM with significant spatial correlation
length(which(grid.MEM.Moran$pvalue <= 0.05))

# Now, we want to select those significant ones but only keep the eigenvectors that have positive eigenvalues.

# These are the significant MEMs
sig.MEM <- grid.invdist.MEM[, c(which(grid.MEM.Moran$pvalue <= 0.05))]

# Now, choose the ones with positive eigenvalues
attr(sig.MEM, which = "values") # Seems like they are all positive

# grid.invdist.MEM.vec <- grid.invdist.MEM$vectors
#
# MEM.Moran.pos <- which(grid.MEM.Moran[,1] > -1/(nrow(grid.invdist.MEM$vectors)-1))
# grid.invdist.MEM.pos <- grid.invdist.MEM.vec[,MEM.Moran.pos]
#
# MEM.Moran.pos.sig <- MEM.Moran.pos[which(grid.MEM.Moran[MEM.Moran.pos,2] <= 0.05)]
#
# grid.invdist.MEM.pos.sig <- grid.invdist.MEM.vec[,MEM.Moran.pos.sig]

grid.invdist.MEM.pos.sig <- sig.MEM

# Save these as an RDS for future HMSC work
saveRDS(grid.invdist.MEM.pos.sig, "sig_ctrl_MEMs.RDS")

head(protists_corrected)


#loop time
results.DF = c("meta", "env_comp", "combo", "spat_comp", "residuals")

for(i in 1:10){


  protists=subset(protists_corrected, landscape == i)
  head(protists)
  mite=protists[14:18] #community data
  mite.env=protists[5] #habitat data, no didinium/mite.d data

  # Standarization of the data
  mite.h <- decostand(mite, "hellinger")


  (mite.varpart <- varpart(mite.h, mite.env, grid.invdist.MEM.pos.sig))
  part <- mite.varpart[[1]]
  ind.fract <- part[[3]]
  adjusted <- ind.fract[,3]
  newrow = c(i, adjusted)

  results.DF = rbind(results.DF,newrow)
  i=i+1
}

results.DF <- data.frame(results.DF)
colnames(results.DF) <- unlist(results.DF[1,])
View(results.DF)

write.table(results.DF, "dispersal_varpart.txt", sep="\t")

