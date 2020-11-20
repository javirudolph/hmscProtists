# Test run of HMSC

library(HMSC)
#library(doParallel)

# Significant MEMs are saved as an RDS file within this directory as:
ctrl.MEMs <- readRDS("old_tests/sig_ctrl_MEMs.RDS")


# Read the data
final.data <- read.table("old_tests/Resetarits&al2018FinalData.txt")

# Subset the data to:
#   - The control of just 1 landscape
# The X matrix = occupancy

prot.occp <- as.matrix(round(final.data[1:36,20:25]))

prot.occp <- ifelse(as.matrix(round(final.data[1:36,20:25])) == 0, 0, 1)


# The env variables
prot.env <- final.data[1:36,c(5,9:10)]
prot.env$sample <- as.factor(prot.env$sample)

#########################################
# HMSC section
#########################################

# Format into hmsc format
protist_hmsc <- as.HMSCdata(Y = prot.occp, X = cbind(prot.env[,2:3], ctrl.MEMs), Random = prot.env$sample)


#protist_hmsc <- as.HMSCdata(Y = prot.occp, X = cbind(prot.env[,2:3], ctrl.MEMs))

# Protist model

protist_model <- hmsc(protist_hmsc, family = "probit", niter = 5000, nburn = 1000, thin = 15)

# Check convergence
protist_model$results$estimation$paramX
dim(protist_model$results$estimation$paramX)

protist_spp_vp <- variPart(protist_model, groupX = c(rep("env", 3), rep("spa", 14)), type = "III", indSite = FALSE, R2adjust = TRUE)

protist_site_vp <- variPart(protist_model, groupX = c(rep("env", 3), rep("spa", 14)), type = "III", indSite = TRUE, R2adjust = TRUE)


# Get the data from VP species:
library(tidyverse)
library(ggtern)

protist_spp_vp %>%
  map(as_tibble) %>%
  bind_cols() %>%
  rownames_to_column() %>%
  set_names(c("species", "c", "b", "a", "e", "f", "d", "g")) %>%
  as.data.frame() -> output
output

rowSums(output[,2:8])

output %>%
  transmute(species = species,
            env = a + f + 0.5 * d + 0.5 * g,
            env = ifelse(env < 0, 0, env),
            spa = b + e + 0.5 * d + 0.5 * g,
            spa = ifelse(spa < 0, 0, spa),
            codist = c,
            codist = ifelse(codist < 0, 0, codist),
            r2 = env + spa + codist) -> plotdata
plotdata

# Plot
plotdata %>%
  ggtern(aes(x = env, z = spa, y = codist, size = r2)) +
  geom_point(aes(color = species), alpha = 0.8) +
  scale_T_continuous(limits=c(0.0,1.0),
                     breaks=seq(0.0,1.0,by=0.1),
                     labels=seq(0.0,1.0,by=0.1)) +
  scale_L_continuous(limits=c(0.0,1),
                     breaks=seq(0,1,by=0.1),
                     labels=seq(0,1,by=0.1)) +
  scale_R_continuous(limits=c(0.0,1.0),
                     breaks=seq(0,1,by=0.1),
                     labels=seq(0,1,by=0.1)) +
  labs(title = "Test Plot",
       x = "E",
       xarrow = "Environment",
       y = "C",
       yarrow = "Co-Distribution",
       z = "S",
       zarrow = "Spatial Autocorrelation") +
  theme_light() +
  theme_showarrows() +
  #scale_colour_brewer(palette = "Set1") +
  #scale_colour_brewer(palette = "Spectral") +
  #scale_color_viridis_d() +
  scale_size_area(limits = c(0,1), breaks = seq(0,1,0.2)) +
  guides(color = guide_legend("none", order = 2),
         size = guide_legend(title = expression(R^2), order = 1)) +
  theme(panel.grid = element_line(color = "darkgrey"),
        axis.text = element_text(size =5),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 12, margin = margin(t = 10, b = -20)),
        tern.axis.arrow = element_line(size = 1))

##########################################################
# Based on the figure, lets try to run the species matrix
##########################################################
library(corrplot)

assoMat <- corRandomEff(protist_model)
siteMean <- apply(assoMat[ , , , 1], 1:2, mean)

siteDrawCol <- matrix(NA, nrow = nrow(siteMean), ncol = ncol(siteMean))
siteDrawCol[which(siteMean > 0.4, arr.ind=TRUE)]<-"red"
siteDrawCol[which(siteMean < -0.4, arr.ind=TRUE)]<-"blue"

# Build matrix of "significance" for corrplot
siteDraw <- siteDrawCol
siteDraw[which(!is.na(siteDraw), arr.ind = TRUE)] <- 0
siteDraw[which(is.na(siteDraw), arr.ind = TRUE)] <- 1
siteDraw <- matrix(as.numeric(siteDraw), nrow = nrow(siteMean), ncol = ncol(siteMean))

Colour <- colorRampPalette(c("blue", "white", "red"))(200)
corrplot(siteMean, method = "color", col = Colour, type = "lower",
         diag = FALSE, p.mat = siteDraw, tl.srt = 45, title = corTitle)

