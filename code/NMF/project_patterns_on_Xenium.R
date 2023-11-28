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

# read the model in
k <- 20
nmf.mod <- readRDS(here("processed-data", "cindy", "NMF", sprintf("visium-nmf-model-k%s", k)))
loadings <- nmf.mod$w

# read in the sfe to project factors onto
sfe <- readRDS(here("processed-data", "cindy", "slide-5434", "Br8667_Mid_SFE_filt.RDS"))

# subset the loadings matrix to only those genes that are in the Xenium panel
i<-intersect(rownames(sfe),rownames(loadings))
loadings<-loadings[rownames(loadings) %in% i,]
loadings<-loadings[match(rownames(sfe),rownames(loadings)),]

sfe <- scuttle::logNormCounts(sfe)
proj<-project(A=logcounts(sfe), w=loadings, L1=0)

factors <- t(proj)
colnames(factors) <- paste0("NMF", 1:k)
colData(sfe) <- cbind(colData(sfe), factors)
pls <- list()
for (i in 1:k){
    patternName <- paste0("NMF", i)
    pls[[i]] <- make_escheR(sfe) |>
        add_fill(patternName) +
        scale_fill_gradient(low = "white", high = "black")+
        ggtitle(paste(unique(spe$sample_id), patternName))
}
