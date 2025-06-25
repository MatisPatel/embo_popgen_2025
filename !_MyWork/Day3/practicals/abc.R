library(abc)
library(tidyverse)

setwd("/home/patel/embo_popgen_2025/!_MyWork/Day3/practicals")

dat <- read_csv("tplit_samples.csv")

txtdensity(dat$fst)
