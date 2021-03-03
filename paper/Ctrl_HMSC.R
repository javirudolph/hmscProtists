# Now that we have the MEMs for control metacommunities, we can fit HMSC

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

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
