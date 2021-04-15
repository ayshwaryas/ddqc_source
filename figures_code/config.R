library(ggplot2)
library(ggridges)
library(cowplot)
library(dplyr)
library(tidyr)
source("plotting.R")

DIR.PREFIX_local <- "/Volumes/easystore/primes_storage/"
DIR.PREFIX <- "/Volumes/scqc/"
PATH <- paste0(DIR.PREFIX_local, "figures/")
OUTPUT.DIR <- paste0(DIR.PREFIX, "output_pg/")