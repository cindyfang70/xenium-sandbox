library(here)
k <- 20
nmf.mod <- readRDS(here("processed-data", "cindy", "NMF", sprintf("visium-nmf-model-k%s.RDS", k)))
loadings <- nmf.mod$w

vis.factors <- t(nmf.mod$h)
colnames(vis.factors) <- paste0("NMF", 1:k)

fname <- paste0("Br8667_Mid_5548", sprintf("-projected-NMF-factors-k%s.RDS", k))
proj <- readRDS(here("processed-data", "cindy", "NMF", fname))

xen.factors <- t(proj)
colnames(xen.factors) <- paste0("NMF", 1:k)
