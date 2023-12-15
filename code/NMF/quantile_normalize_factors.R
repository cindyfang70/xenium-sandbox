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
library(preprocessCore)

# Try quantile normalization on the Xenium factors and see if that helps with 
# detecting white matter etc. 


# read the model in
k <- 20
nmf.mod <- readRDS(here("processed-data", "cindy", "NMF", sprintf("visium-nmf-model-k%s.RDS", k)))
loadings <- nmf.mod$w
vis.factors <- t(nmf.mod$h)

# read in the sfe to project factors onto
sfe <- readRDS(here("processed-data", "cindy", "slide-5434", "slide5434-all-samples-spe-with-banksy.RDS"))

sfe_list <- lapply(unique(sfe$region_id), function(x) 
    sfe[, sfe$region_id == x])

#paste0(unique(sfe$region_id), sprintf("-raw-projected-NMF-factors-k%s.RDS", k))

# make violin plots for the visium factors
colnames(vis.factors) <- paste0("NMF", 1:k)
vis.factors_long <- pivot_longer(as.data.frame(vis.factors), cols=everything())
vis.factors_long$type <- "Visium"

plist <- list()
plist_normed <- list()
for(i in seq_along(sfe_list)){
    sfe <- sfe_list[[i]]
    fct_name <- paste0(unique(sfe$region_id), sprintf("-raw-projected-NMF-factors-k%s.RDS", k))
    proj <- readRDS(here("processed-data", "cindy", "NMF", fct_name))
    factors <- as.data.frame(t(proj))
    colnames(factors) <- paste0("NMF", 1:k)
    factors_long <- pivot_longer(factors, cols=everything())
    factors_long$type <- "Xenium"
    
    factors_long <- rbind(factors_long, vis.factors_long)
    
    # make violin plots of each of the factors before quantile normalization
    plist[[i]] <- ggplot(factors_long, aes(y=value, x=type, colour=type))+
        geom_violin()+
        facet_wrap(~ name)+
        coord_trans(y="log1p")+
        theme(legend.position="none")+
        ggtitle(unique(sfe$region_id))+
        cowplot::theme_cowplot()
    
    # try doing quantile normalization using preprocessCore::normalize.quantiles
    xen.factors.normed <- matrix(,ncol=k, nrow=nrow(factors))
    for (j in 1:k){
        vis <- vis.factors[,j]
        xen <- factors[,j]
        xen.factors.normed[,j] <- normalize.quantiles.use.target(x=as.matrix(xen), target=vis)
    }
    xen.factors.normed <- as.data.frame(xen.factors.normed)
    colnames(xen.factors.normed) <- paste0("NMF", 1:k)
    xen.factors.normed_long <- pivot_longer(xen.factors.normed, cols=everything())
    xen.factors.normed_long$type <- "Xenium"
    
    all_fcts_normed <- rbind(xen.factors.normed_long, vis.factors_long)
    # plot the normed factors dists
    plist_normed[[i]] <- ggplot(all_fcts_normed , aes(x=type, y=value, colour=type))+
        geom_violin()+
        facet_wrap(~ name)+
        coord_trans(y = "log1p")+
        theme(legend.position="none")+
        ggtitle(paste0("Normed ", unique(sfe$region_id)))+
        cowplot::theme_cowplot()
}

pdf(here("plots", "NMF", sprintf("slide5434-NMF-k%s-factorDensityPlots.pdf",k)))
for(i in seq_along(plist)){
    print(plist[[i]])
    print(plist_normed[[i]])
}
dev.off()

