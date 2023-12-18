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

# if running locally:
args <- c("snRNA-seq", "20", 
          "processed-data/cindy/NMF/snRNA-seq/snRNAseq-nmf-model-k20.RDS",
          "processed-data/cindy/NMF/snRNA-seq/visium-nmf-snRNA-seq-k20-multinom-model.RDS")

# read in the data
model_type <- args[[1]]
k <- as.numeric(args[[2]])
model <- readRDS(args[[3]])

## predict using the multinomial model
multinom <- readRDS(args[[4]])

slide_number <- "5434"
sfe <- readRDS(here("processed-data", "cindy", sprintf("slide-%s", slide_number), 
                    sprintf("slide%s-all-samples-spe-with-banksy.RDS", slide_number)))

sfe_list <- lapply(unique(sfe$region_id), function(x) 
    sfe[, sfe$region_id == x])

nmf.mod <- readRDS(args[[3]])
ref.factors <- t(nmf.mod$h)
colnames(ref.factors) <- paste0("NMF", 1:k)

for (i in seq_along(sfe_list)){
    sfe <- sfe_list[[i]]
    fname <-  paste0(unique(sfe$region_id), 
                     sprintf("-raw-projected-NMF%sfactors-k%s.RDS", 
                             model_type, k))
    
    proj <- readRDS(here("processed-data", "cindy", "NMF", model_type,
                            fname))
    
    factors <- t(proj)
    colnames(factors) <- paste0("NMF", 1:k)
    
    
    # quantile normalize the factors
    xen.factors.normed <- matrix(,ncol=k, nrow=nrow(factors))
    for (j in 1:k){
        ref <- ref.factors[,j]
        xen <- factors[,j]
        xen.factors.normed[,j] <- normalize.quantiles.use.target(
            x=as.matrix(xen), target=ref)
    }
    xen.factors.normed <- as.data.frame(xen.factors.normed)
    colnames(xen.factors.normed) <- paste0("NMF", 1:k)
    
    probs <- predict(multinom, newdata=xen.factors.normed, type='probs', 
                     na.action=na.exclude)
    
    print(unique(probs))
    preds <- unlist(lapply(1:nrow(probs), function(xx){
        colnames(probs)[which.max(probs[xx,])]
    })) 
    
}


