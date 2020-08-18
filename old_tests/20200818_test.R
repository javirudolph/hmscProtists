# We are testing that HMSC works. There was a recent update, R 4.0 and it doesn't work with HMSC, so now trying to run this with previous version of R.

# I keep having problems with Blanchet's HMSC
# Switching and giving a try to the one from Otso's group

library(Hmsc)
library(corrplot)
set.seed(1)


# Starting with one of the communities from Emlyn\
all_data <- read.table(file = "old_tests/Resetarits&al2018FinalData.txt")

# The first 36 rows correspond to one of the control metacommunities, so all 36 patches are present.
Y <- all_data[1:36, 20:26] # density/species matrix
XData <- all_data[1:36, c(5,9,10)] # environmental variables

model <- Hmsc(Y = Y, XData = XData)
