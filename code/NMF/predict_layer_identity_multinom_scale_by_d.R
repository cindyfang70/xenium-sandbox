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
library(preprocessCore)

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

pls <- list()
for (i in seq_along(sfe_list)){
    sfe <- sfe_list[[i]]
    fname <-  paste0(unique(sfe$region_id), 
                     sprintf("-raw-projected-NMF%sfactors-k%s.RDS", 
                             model_type, k))
    
    proj <- readRDS(here("processed-data", "cindy", "NMF", model_type,
                         fname))
    
    factors <- t(proj)
    colnames(factors) <- paste0("NMF", 1:k)
    
    
    # scale the factors by model$d from NMF 
    xen.factors.normed <- factors/model$d

    xen.factors.normed <- as.data.frame(xen.factors.normed)
    colnames(xen.factors.normed) <- paste0("NMF", 1:k)
    
    probs <- predict(multinom, newdata=xen.factors.normed, type='probs', 
                     na.action=na.exclude)
    
    print(unique(probs))
    preds <- unlist(lapply(1:nrow(probs), function(xx){
        colnames(probs)[which.max(probs[xx,])]
    })) 
    
    preds_name <- sprintf("preds_from_%s_NMF_k%s", model_type, k)
    colData(sfe)[preds_name] <- preds
    sfe_list[[i]] <- sfe
    
    pls[[i]] <- make_escheR(sfe, y_reverse=FALSE)%>%
        add_ground(preds_name)
    
}

# plot each layer individually

pdf(here("plots", "NMF", model_type, sprintf("predicted_layers_scale_by_d_%s_k%s_%s.pdf", 
                                             model_type, k, slide_number)), height=25, width=15)
for (i in seq_along(sfe_list)){
    sfe <- sfe_list[[i]]
    layers <- paste0("L", 1:6)
    layers <- c(layers, "WM")
    
    layer_plts <- list()
    for (l in 1:length(layers)){
        sfe$isLayer <- grepl(layers[[l]], colData(sfe)[[preds_name]])
        layer_plts[[l]] <- make_escheR(sfe, y_reverse=FALSE)%>%
            add_fill("isLayer")+
            scale_fill_manual(values=c("TRUE"="red", "FALSE"="grey"))+
            ggtitle(paste(unique(sfe$region_id), layers[[l]]))
    }
    
    layer_plts[[length(layers) + 1]] <- pls[[i]]
    do.call(gridExtra::grid.arrange, c(layer_plts, ncol=2))
}
dev.off()
