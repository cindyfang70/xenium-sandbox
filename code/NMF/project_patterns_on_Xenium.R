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
sfe <- readRDS(here("processed-data", "cindy", "slide-5434", "slide5434-sfe-with-banksy.RDS"))

sfe_list <- lapply(unique(sfe$region_id), function(x) 
    sfe[, sfe$region_id == x])
# subset the loadings matrix to only those genes that are in the Xenium panel
i<-intersect(rownames(sfe),rownames(loadings))
loadings<-loadings[rownames(loadings) %in% i,]
loadings<-loadings[match(rownames(sfe),rownames(loadings)),]


pdf(here("plots","NMF", "project-factors-xenium-5434.pdf"),
    height=35, width=35)
for (i in seq_along(sfe_list)){
    sfe <- sfe_list[[i]]
    # log normalize and then project using the loadings
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
            ggtitle(paste(unique(spe$region_id), patternName))
    }
    
    do.call(gridExtra::grid.arrange, c(pls, ncol=3))
}
dev.off()
