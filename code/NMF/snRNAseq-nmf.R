library(RcppML)
library(Matrix)
library(Voyager)
library(CoGAPS)
library(projectR)
library(spatialLIBD)
library(escheR)
library(here)
library(tidyverse)
library(ggforce)

# Use RcppML::nmf to perform non-negative matrix factorization on the manually
# annotated Visium dataset to get a loadings matrix and a factors matrix with 
# k factors. Plot each of those k factors as well as the layer labels to see
# if any of them are associated with a particular layer.

args <- commandArgs(trailingOnly = TRUE)
spe_path <- args[[1]]
k <- as.numeric(args[[2]])

# get the manually annotated visium data
ehub <- ExperimentHub::ExperimentHub()
vis_anno <- readRDS(spe_path)
# use layer_guess_reordered as the manual annotations

# NMF on the annotated visium data
A <- logcounts(vis_anno) # using logcounts because there are multiple datasets
model <- RcppML::nmf(A, k = k, seed=1237)
patterns <- t(model$h) # these are the factors

rownames(model$w) <- rowData(vis_anno)$gene_name

saveRDS(model, here("processed-data", "cindy", "NMF", sprintf("snRNAseq-nmf-model-k%s", k)))
