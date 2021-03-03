# Now that we have the MEMs for control metacommunities, we can fit HMSC

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(purrr)
library(tibble)
library(ggtern)
library(corrplot)

resetarits_data <- read.table(here::here("data", "Resetarits&al2018FinalData.txt"))
head(resetarits_data)


# Filter the data so we have only control metacommunities and the densities

resetarits_data %>%
  dplyr::filter(., type == "E") %>%
  dplyr::select(., landscape, ends_with("_density")) -> ctrl_density

resetarits_data %>%
  dplyr::filter(., type == "E") %>%
  dplyr::select(., landscape, habitat.bi, age) -> ctrl_env

# environmental variable, the 0 is algae

# Load the MEMs for control metacommunity:

ctrl_mem <- readRDS("paper/ctrl_mem.RDS")

library(HMSC)

onectrl <- ctrl_density %>%
  dplyr::filter(., landscape==1) %>%
  dplyr::select(., -landscape) %>%
  as.matrix()

onectrl[onectrl > 0] <- 1

env_vars <- ctrl_env %>%
  dplyr::filter(., landscape==1) %>%
  dplyr::select(., -landscape)

ashmsc <- as.HMSCdata(Y = onectrl, X = cbind(env_vars, ctrl_mem),
                      Random = as.factor(1:36),
                      scaleX = TRUE, interceptX = TRUE)

model <- hmsc(ashmsc, family = "probit", niter = 5000, nburn = 1000, thin = 1)

# Need to write this as a function later, but need to get the model as an mcmc object for the diagnostics

get_mcmc_obj <- function(model){
  parameters <- "meansParamX"
  x <- model

  paramMCMCMat <- x$results$estimation[[parameters]]

  ### Name rows and columns of matrix
  rownames(paramMCMCMat) <- rownames(x$results$estimation[[parameters]])
  colnames(paramMCMCMat) <- colnames(x$results$estimation[[parameters]])

  paramBurnMCMC <- x$results$burning[[parameters]]
  rownames(paramBurnMCMC) <- rownames(x$results$burning[[parameters]])
  colnames(paramBurnMCMC) <- colnames(x$results$burning[[parameters]])
  paramMCMCMat <- rbind(paramBurnMCMC,paramMCMCMat)

  thinSet <- as.numeric(sub("iter", "", rownames(paramMCMCMat)[1:2]))
  modelThin <- thinSet[2] - thinSet[1]

  result <- mcmc(data = paramMCMCMat, thin = modelThin)
}

forconv <- get_mcmc_obj(model)

raftery.diag(forconv)


# Variation partitioning
# Species

vpspp <- variPart(model, groupX = c(rep("env",3),rep("spa",14)),
                  type = "III", R2adjust = TRUE)

# format output
vpspp %>%
  map(as_tibble) %>%
  bind_cols() %>%
  rownames_to_column() %>%
  set_names(c("species", "c", "b", "a", "e", "f", "d", "g")) %>%
  transmute(species = species,
            env = a + f + 0.5 * d + 0.5 * g,
            env = ifelse(env < 0, 0, env),
            spa = b + e + 0.5 * d + 0.5 * g,
            spa = ifelse(spa < 0, 0, spa),
            codist = c,
            codist = ifelse(codist < 0, 0, codist),
            r2 = env + spa + codist) -> sppdata

species_plot <- function(data, plotMain = NULL, colorVar = NULL, colorLegend = "none"){
  data %>%
    ggtern(aes(x = env, z = spa, y = codist, size = r2)) +
    geom_point(aes_string(color = colorVar), alpha = 0.8) +
    scale_T_continuous(limits=c(0,1.0),
                       breaks=seq(0,1,by=0.1),
                       labels=seq(0,1,by=0.1)) +
    scale_L_continuous(limits=c(0.0,1),
                       breaks=seq(0,1,by=0.1),
                       labels=seq(0,1,by=0.1)) +
    scale_R_continuous(limits=c(0.0,1.0),
                       breaks=seq(0,1,by=0.1),
                       labels=seq(0,1,by=0.1)) +
    labs(title = plotMain,
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
    guides(color = guide_legend(colorLegend, order = 2),
           size = guide_legend(title = expression(R^2), order = 1)) +
    theme(panel.grid = element_line(color = "darkgrey"),
          axis.title = element_text(size = 8))
}

species_plot(sppdata)

# Sites
vpsite <- variPart(model, groupX = c(rep("env",3),rep("spa",14)), indSite = TRUE,
                  type = "III", R2adjust = TRUE)

nspp <- 7
c <- rowSums(vpsite$overlap1[,,1])/nspp
b <- rowSums(vpsite$overlap1[,,2])/nspp
a <- rowSums(vpsite$overlap1[,,3])/nspp

e <- rowSums(vpsite$overlap2[,,1])/nspp
f <- rowSums(vpsite$overlap2[,,2])/nspp
d <- rowSums(vpsite$overlap2[,,3])/nspp

g <- rowSums(vpsite$overlap3)/nspp

env <- a + f + 1/2 * d + 1/2 * g
env <- ifelse(env < 0, 0, env)
spa <- b + e + 1/2 * d + 1/2 * g
spa <- ifelse(spa < 0, 0, spa)
random <- c
codist <- ifelse(random < 0, 0, random)
r2 <- env + spa + codist

sitedata <- cbind.data.frame(env, spa, codist, r2)

sites_plot <- function(data, plotMain = NULL, colorVar = NULL, colorLegend = "none"){
  data %>%
    ggtern(aes(x = env, z = spa, y = codist, size = r2)) +
    geom_point(aes_string(color = colorVar), alpha = 0.6) +
    scale_T_continuous(limits=c(0,1.0),
                       breaks=seq(0,1,by=0.1),
                       labels=seq(0,1,by=0.1)) +
    scale_L_continuous(limits=c(0.0,1),
                       breaks=seq(0,1,by=0.1),
                       labels=seq(0,1,by=0.1)) +
    scale_R_continuous(limits=c(0.0,1.0),
                       breaks=seq(0,1,by=0.1),
                       labels=seq(0,1,by=0.1)) +
    labs(title = plotMain,
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
    scale_size_area(limits = c(0, 0.003), breaks = seq(0, 0.003, 0.0005)) +
    guides(color = guide_colorbar(colorLegend),
           size = guide_legend(title = expression(R^2))) +
    theme(panel.grid = element_line(color = "darkgrey"),
          axis.title = element_text(size = 8))
}

sites_plot(sitedata)


# Species correlation

interaction_plot <- function(model, corTitle = NULL){

  assoMat <- corRandomEff(model)
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

}

interaction_plot(model)


