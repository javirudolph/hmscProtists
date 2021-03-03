# Testing Ovaskainen's paper
# Using these vignettes:
# https://cran.r-project.org/web/packages/Hmsc/vignettes/

library(tidyverse)
library(Hmsc)


resetarits_data <- read.table(here::here("data", "Resetarits&al2018FinalData.txt"))
head(resetarits_data)


# Filter the data so we have only control metacommunities and the densities

resetarits_data %>%
  dplyr::filter(., type == "E") %>%
  dplyr::select(., landscape, ends_with("_density")) -> ctrl_density

# The scale for the densities are super different for chilomonas
# Should use log? Yes, much better.
resetarits_data %>%
  dplyr::filter(., type == "E") %>%
  dplyr::select(., landscape, habitat, age, ends_with("_density")) %>%
  pivot_longer(., ends_with("_density"), names_to = "specie", values_to = "density") ->toplot

toplot %>%
  ggplot(., aes(x = specie, y = density, color = habitat)) +
  facet_wrap(~landscape) +
  geom_point() +
  scale_y_log10()


resetarits_data %>%
  dplyr::filter(., type == "E") %>%
  dplyr::select(., landscape, habitat.bi, age) -> ctrl_env

# environmental variable, the 0 is algae

# Load the MEMs for control metacommunity:

ctrl_mem <- readRDS("paper/ctrl_mem.RDS")

onectrl <- ctrl_density %>%
  dplyr::filter(., landscape==1) %>%
  dplyr::select(., -landscape) %>%
  as.matrix()

binary.ctrl <- onectrl
binary.ctrl[binary.ctrl > 0] <- 1

env_vars <- ctrl_env %>%
  dplyr::filter(., landscape==1) %>%
  dplyr::select(., -landscape)


# Using Hmsc
y <- onectrl
model <- Hmsc(Y = y, XData = env_vars, distr = "poisson")

nChains = 2

test.run = FALSE
if (test.run){
  #with this option, the vignette runs fast but results are not reliable
  thin = 1
  samples = 10
  transient = 5
  verbose = 5
} else {
  #with this option, the vignette evaluates slow but it reproduces the results of the
  #.pdf version
  thin = 5
  samples = 1000
  transient = 500*thin
  verbose = 500*thin
}

m = sampleMcmc(model, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)

mpost = convertToCodaObject(m)
summary(mpost$Beta)


preds = computePredictedValues(m)
evaluateModelFit(hM=m, predY=preds)

# HMSC applies the Bayesian framework (and thus e.g. assumes prior distributions and yields credible intervals)

# Posterior trace plots of Beta parameters
plot(mpost$Beta)

# Good results
#   -The chains should look essentially identical
#   -Good mixing of the chains. Go fast up and down with any apparent correlation
#   -Reach stationary distribution. First half of recorded iterations looks statistically identical to the second half.

effectiveSize(mpost$Beta)
# I am not sure of how to interpret these sample sizes

# The potential scale reduction factors are very close to one, which indicates that two chains gave consistent results, as suggested by visual inspection of the trace plots
gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
# I guess except for the last species. Chilomonas.

# In summary, the MCMC diagnostics did not indicate any problems with MCMC convergence. This means
# that the posterior sample is likely to be representative of the true posterior distribution, and thus the inference
# from the model can be trusted.

# Try to stick with 1000 samples, and vary the thinning














