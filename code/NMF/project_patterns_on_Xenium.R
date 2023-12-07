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

# read the model in
k <- 20
nmf.mod <- readRDS(here("processed-data", "cindy", "NMF", sprintf("visium-nmf-model-k%s.RDS", k)))
loadings <- nmf.mod$w
vis.factors <- t(nmf.mod$h)

# read in the sfe to project factors onto
sfe <- readRDS(here("processed-data", "cindy", "slide-5548", "slide5548-all-samples-spe-with-banksy.RDS"))

sfe_list <- lapply(unique(sfe$region_id), function(x) 
    sfe[, sfe$region_id == x])
# subset the loadings matrix to only those genes that are in the Xenium panel
i<-intersect(rownames(sfe),rownames(loadings))
loadings<-loadings[rownames(loadings) %in% i,]
loadings<-loadings[match(rownames(sfe),rownames(loadings)),]


rescaleFactors <- function(query_vec, ref_vec){
    scale_max <- max(ref_vec)
    scale_min <- min(ref_vec)
    
    query_scaled <- scales::rescale(query_vec, to=c(scale_min, scale_max))
    return(query_scaled)
}


for (j in seq_along(sfe_list)){
    sfe <- sfe_list[[j]]
    # log normalize and then project using the loadings
    sfe <- scuttle::logNormCounts(sfe)
    proj<-project(A=logcounts(sfe), w=loadings, L1=0)
    
    fname <- paste0(unique(sfe$region_id), sprintf("-projected-NMF-factors-k%s.RDS", k))
    saveRDS(proj, here("processed-data", "cindy", "NMF", fname))
    
    factors <- t(proj)
    colnames(factors) <- paste0("NMF", 1:k)
    
    # rescale the factors to match the Visium factors' range.
    xen.factors.scaled <- matrix(,ncol=k, nrow=nrow(factors))
    for(i in 1:k){
        xen.factors.scaled[,i] <- rescaleFactors(factors[,i], vis.factors[,i])
    }
    colnames(xen.factors.scaled) <- paste0("NMF", 1:k)
    colData(sfe) <- cbind(colData(sfe), xen.factors.scaled)
    
    sfe_list[[j]] <- sfe
}

pdf(here("plots","NMF", sprintf("project-factors-scaled-xenium-k%s-5548.pdf", k)),
    height=15, width=45)
for (i in 1:k){
    pls <- list()
    for (j in seq_along(sfe_list)){
        sfe <- sfe_list[[j]]
        patternName <- paste0("NMF", i)
        pls[[j]] <- make_escheR(sfe, y_reverse=FALSE) |>
            add_ground("clust_M0_lam0.9_k50_res0.4")|>
            add_fill(patternName) +
            scale_fill_gradient(low = "white", high = "black")+
            ggtitle(paste(unique(sfe$region_id), patternName))
    }
    
    do.call(gridExtra::grid.arrange, c(pls, ncol=3))
}
dev.off()

## predict using the multinomial model
multinom <- readRDS(here("processed-data", "cindy", "NMF", 
                         sprintf("visium-nmf-k%s-multinom-model.RDS",k)))
library(nnet)
library(MASS)
plist <- list()
for (i in seq_along(sfe_list)){
    sfe <- sfe_list[[i]]
    factors <- colData(sfe)[grepl("NMF", colnames(colData(sfe)))]
    
    # predict using the multinomial model
    probs <- predict(multinom, newdata=factors, type='probs', 
                     na.action=na.exclude)
    print(unique(probs))
    preds <- unlist(lapply(1:nrow(probs), function(xx){
             colnames(probs)[which.max(probs[xx,])]
        })) 
    maxprobs <- unlist(lapply(1:nrow(probs), function(xx){
        max(probs[xx,])
    })) 
    
    sfe$max_probs <- maxprobs
    
    sfe$predicted_layers <- preds
    
    sfe_list[[i]] <- sfe
    sfe$counts_MOBP <- counts(sfe)[which(rownames(sfe)=="MOBP"),]
    plist[[i]] <- make_escheR(sfe, y_reverse=FALSE)|>
        add_ground("predicted_layers") |>
        add_fill("counts_MOBP") +
        scale_fill_gradient(low="white", high="black")

}

pdf(here("plots", "NMF", "predicted-layers-nmf-scaled-5548.pdf"),
    height=15, width=35)
for(i in seq_along(plist)){
    p <-  make_escheR(sfe_list[[i]], y_reverse=FALSE)|>
        add_ground("clust_M0_lam0.9_k50_res0.4")
    
    p2 <-  make_escheR(sfe_list[[i]], y_reverse=FALSE)|>
        add_ground("predicted_layers")
    print(plist[[i]] + p + p2)
}
dev.off()



# 
# sfe <- sfe_list[[1]]
# pred <- unlist(lapply(strsplit(sfe$clust_pred, split="r"), "[", 2))
# pred[which(is.na(pred))] <- "7"
# 
# sfe$clust_pred <- as.numeric(pred)
# 
# sfe <- Banksy::connectClusters(sfe)
# 
# p1 <- make_escheR(sfe, y_reverse=FALSE) |>
#     add_fill("clust_pred_num")+
#     scale_fill_discrete()
# 
# p2 <- make_escheR(sfe, y_reverse=FALSE)|>
#     add_fill("clust_M0_lam0.9_k50_res0.4")+
#     scale_fill_discrete()
# p1+p2
