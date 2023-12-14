library(RcppML)
library(Matrix)
#library(Voyager)
#library(CoGAPS)
library(projectR)
library(spatialLIBD)
library(escheR)
library(here)
library(tidyverse)
library(ggforce)

# Try quantile normalization on the Xenium factors and see if that helps with 
# detecting white matter etc. 


# read the model in
k <- 20

# read in the sfe to project factors onto
sfe <- readRDS(here("processed-data", "cindy", "slide-5434", "slide5434-all-samples-spe-with-banksy.RDS"))

sfe_list <- lapply(unique(sfe$region_id), function(x) 
    sfe[, sfe$region_id == x])