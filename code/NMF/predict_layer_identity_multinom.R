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
library(nnet)
library(MASS)

################################################################################
# Use the trained multinomial model to predict layer identity in the Xenium data
################################################################################

## predict using the multinomial model
multinom <- readRDS(here("processed-data", "cindy", "NMF", 
                         sprintf("visium-bayesspace-nmf-k%s-multinom-model.RDS",k)))

slide_number <- "5434"
sfe <- readRDS(here("processed-data", "cindy", sprintf("slide-%s", slide_number), 
                    sprintf("slide%s-all-samples-spe-with-banksy.RDS", slide_number)))

sfe_list <- lapply(unique(sfe$region_id), function(x) 
    sfe[, sfe$region_id == x])



