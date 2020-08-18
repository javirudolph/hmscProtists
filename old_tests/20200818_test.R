# We are testing that HMSC works. There was a recent update, R 4.0 and it doesn't work with HMSC, so now trying to run this with previous version of R.

# library(devtools)
#
# devtools::install_github('guiblanchet/HMSC')
library(HMSC)


# Load the data from Emlyn's work
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

# Fit HMSC

# Format into hmsc format
protist_hmsc <- as.HMSCdata(Y = prot.occp, X = cbind(prot.env[,2:3], ctrl.MEMs), Random = prot.env$sample)

protist_model <- hmsc(protist_hmsc, family = "probit",
                      niter = 1000, nburn = 500, thin = 15)


# Check for convergence

