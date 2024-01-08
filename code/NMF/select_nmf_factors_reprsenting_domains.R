library(RcppML)
library(Matrix)
library(spatialLIBD)
library(escheR)
library(here)
library(tidyverse)
library(ggforce)

################################################################################
# After performing NMF on the reference/source dataset, pick out the top m 
# factors for each spatial domain based on correlation.
################################################################################

model_type <- "bayesspace"
k <- 100
m <- 5
cors <- readRDS(here("processed-data", "cindy", "NMF", model_type,
     sprintf("%s-nmf-correlations-k%s.RDS", model_type, k)))

layer_cor <- cors$layer_cor
n_domains <- ncol(layer_cor)

selected_factors <- c()
for (i in 1:n_domains){
    domain <- layer_cor[,i]
    ord_domain <- domain[order(domain, decreasing=TRUE)]
    top_m_max <- names(ord_domain[1:m])
    selected_factors <- c(selected_factors, top_m_max)
}
