# Testing Ovaskainen's paper
# Using these vignettes:
# https://cran.r-project.org/web/packages/Hmsc/vignettes/

library(tidyverse)
library(Hmsc)
library(corrplot)


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
  pivot_longer(., ends_with("_density"), names_to = "specie", values_to = "density") %>%
  mutate(habitat = ifelse(habitat == "C", "Algae", "Bacteria")) ->toplot

toplot %>%
  ggplot(., aes(x = specie, y = density, color = habitat)) +
  #facet_wrap(~landscape) +
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
y <- log(onectrl+1)
model <- Hmsc(Y = y, XData = env_vars, distr = "poisson")

# Lognormal poisson would allow for more variation around the expectation, unlike the poisson. Perhaps works for the chilomonas thing too

# I guess except for the last species. Chilomonas. Fixed when y <- onectrl
model <- Hmsc(Y = y, XData = env_vars, distr = "lognormal poisson")

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
gelman.diag(mpost$Beta,multivariate=FALSE)$psrf #using log scale for density

# In summary, the MCMC diagnostics did not indicate any problems with MCMC convergence. This means
# that the posterior sample is likely to be representative of the true posterior distribution, and thus the inference
# from the model can be trusted.

# Try to stick with 1000 samples, and vary the thinning



# Vignette N2. Low dimensional multivariate models

model <- Hmsc(Y = y, XData = env_vars, XFormula = ~habitat.bi + age, distr = "lognormal poisson")


nChains = 2
thin = 5
samples = 1000
transient = 500*thin
verbose = 500*thin

m <- sampleMcmc(model, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)

# Check convergence
mpost <- convertToCodaObject(m)
summary(mpost$Beta)

# Apparently effective sample sizes need to be high. Here they are not, so not sure what that means
effectiveSize(mpost$Beta)

plot(mpost$Beta)

gelman.diag(mpost$Beta,multivariate=FALSE)$psrf

preds = computePredictedValues(m)
evaluateModelFit(hM=m, predY=preds)


# Figure 1: Histograms of effective sample sizes and potential scale reduction factors (psrf) for Beta parameters
par(mfrow=c(1,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")


#To assess the model’s explanatory power, we apply the evaluateModelFit function to the posterior predictive
# distribution simulated by the function computePredictedValues.
preds = computePredictedValues(m)
evaluateModelFit(hM = m, predY = preds)

# Evaluate model's predictive power with two-fold cross validation
partition = createPartition(m, nfolds = 2)
preds = computePredictedValues(m, partition = partition)

evaluateModelFit(hM = m, predY = preds)

# Plot for positive or negative responses to environmental parameters included.
# White is no strong statistical support
postBeta = getPostEstimate(m, parName = "Beta")
plotBeta(m, post = postBeta, param = "Support", supportLevel = 0.95)

#  estimating species-to-species residual associations
#  need to include a random effect at the level of the sampling unit
n=36
Y <- y
XData <- env_vars
names(XData) <- c("x1", "x2")
studyDesign = data.frame(sample = as.factor(1:n))
rL = HmscRandomLevel(units = studyDesign$sample)
m = Hmsc(Y = Y, XData = XData, XFormula = ~x1+x2,
         studyDesign = studyDesign, ranLevels = list(sample = rL))
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)

#In a model with random effects, it is important to look at the convergence diagnostics not only for the
#β parameters, but also for the Ω parameters. The matrix Ω is the matrix of species-to-species residual
#covariances.
mpost = convertToCodaObject(m)
par(mfrow=c(2,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")
hist(effectiveSize(mpost$Omega[[1]]), main="ess(omega)")
hist(gelman.diag(mpost$Omega[[1]], multivariate=FALSE)$psrf, main="psrf(omega)")

# Visuallize the beta parameters
postBeta = getPostEstimate(m, parName="Beta")
plotBeta(m, post=postBeta, param="Support", supportLevel = 0.95)

# In addition to the β parameters related to the fixed effects, we can now look at the estimated species-to-species
# associations. We extract them from the model object with the computeAssociations function, which also
# converts the covariances to the more convenient scale of correlation (ranging from -1 to +1). In the script
# below, we choose to plot only those associations for which the posterior probability for being negative or
# positive is at least 0.95. There is no specific function for plotting species-to-species associations in HMSC,
# but such plots can be generated straightforwardly with the corrplot function.

OmegaCor = computeAssociations(m)
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
corrplot(toPlot, method = "color",
         col = colorRampPalette(c("blue","white","red"))(200),
         title = paste("random effect level:", m$rLNames[1]), mar=c(0,0,1,0))
# It basically shows no correlations between species

# Thinking I could also try this with the covariates separately?
