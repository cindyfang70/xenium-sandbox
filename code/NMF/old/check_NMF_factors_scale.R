library(here)
library(scales)
k <- 20
nmf.mod <- readRDS(here("processed-data", "cindy", "NMF", sprintf("visium-nmf-model-k%s.RDS", k)))
loadings <- nmf.mod$w

vis.factors <- t(nmf.mod$h)
colnames(vis.factors) <- paste0("NMF", 1:k)

fname <- paste0("Br8667_Mid_5548", sprintf("-projected-NMF-factors-k%s.RDS", k))
proj <- readRDS(here("processed-data", "cindy", "NMF", fname))

xen.factors <- t(proj)
colnames(xen.factors) <- paste0("NMF", 1:k)


max(xen.factors)
# [1] 6138.568
max(vis.factors)
# [1] 0.0002853049 

# rescale the xenium factors to be the same scale as visium

rescaleFactors <- function(query_vec, ref_vec){
    scale_max <- max(ref_vec)
    scale_min <- min(ref_vec)
    
    query_scaled <- scales::rescale(query_vec, to=c(scale_min, scale_max))
    return(query_scaled)
}

rescaleFactors(xen.factors[,1], vis.factors[,1])
xen.factors.scaled <- matrix(,ncol=k, nrow=nrow(xen.factors))
for(i in 1:k){
    xen.factors.scaled[,i] <- rescaleFactors(xen.factors[,i], vis.factors[,i])
}
colnames(xen.factors.scaled) <- paste0("NMF", 1:k)

multinom <- readRDS(here("processed-data", "cindy", "NMF", 
                         sprintf("visium-nmf-k%s-multinom-model.RDS",k)))

library(nnet)
preds <- predict(multinom, newdata=xen.factors.scaled, type="probs")
