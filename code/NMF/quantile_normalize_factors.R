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

# specify which slide and which type of labels
slide_number <- "5548"
label_type <- "-"


# read the model in
k <- 20
nmf.mod <- readRDS(here("processed-data", "cindy", "NMF", sprintf("visium%snmf-model-k%s.RDS", label_type, k)))
loadings <- nmf.mod$w
vis.factors <- t(nmf.mod$h)

# read in the sfe to project factors onto
sfe <- readRDS(here("processed-data", "cindy", sprintf("slide-%s", slide_number), 
                    sprintf("slide%s-all-samples-spe-with-banksy.RDS", slide_number)))

sfe_list <- lapply(unique(sfe$region_id), function(x) 
    sfe[, sfe$region_id == x])

#paste0(unique(sfe$region_id), sprintf("-raw-projected-NMF-factors-k%s.RDS", k))

# make violin plots for the visium factors
colnames(vis.factors) <- paste0("NMF", 1:k)
vis.factors_long <- pivot_longer(as.data.frame(vis.factors), cols=everything())
vis.factors_long$type <- "Visium"

plist <- list()
plist_normed <- list()
xen.fcts.normed <- list()
for(i in seq_along(sfe_list)){
    sfe <- sfe_list[[i]]
    fct_name <- paste0(unique(sfe$region_id), 
                       sprintf("-raw-projected-NMF%sfactors-k%s.RDS", 
                               label_type, k))
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
    
    xen.fcts.normed <- rlist::list.append(xen.fcts.normed, xen.factors.normed)
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

pdf(here("plots", "NMF", sprintf("slide%s-NMF%sk%s-factorDensityPlots.pdf", 
                                 slide_number, label_type,k)))
for(i in seq_along(plist)){
    print(plist[[i]])
    print(plist_normed[[i]])
}
dev.off()

# Now try to use the multinomial model to predict the layer identity
multinom <- readRDS(here("processed-data", "cindy", "NMF", 
                         sprintf("visium%snmf-k%s-multinom-model.RDS",label_type,k)))
library(nnet)
library(MASS)

layer_pl <- list()
for (i in seq_along(sfe_list)){
    sfe <- sfe_list[[i]]
    fcts <- xen.fcts.normed[[i]]
    
    probs <- predict(multinom, newdata=fcts, type="probs", na.action=na.exclude)
    
    preds <- unlist(lapply(1:nrow(probs), function(xx){
        colnames(probs)[which.max(probs[xx,])]})) 
    
    maxprobs <- unlist(lapply(1:nrow(probs), function(xx){
        max(probs[xx,])})) 
    
    preds_name <- sprintf("predicted_layers%s_NMF_k%s_manual_annot_quant_norm", 
                          label_type, k)
    
    colData(sfe)[preds_name] <- preds
    
    sfe_list[[i]] <- sfe
    
    sfe$counts_MOBP <- counts(sfe)[which(rownames(sfe)=="MOBP"),]
    layer_pl[[i]] <- make_escheR(sfe, y_reverse=FALSE)|>
        add_ground(preds_name) |>
        add_fill("counts_MOBP") +
        scale_fill_gradient(low="white", high="black")+
        ggtitle(unique(sfe$region_id))+
        # theme(legend.title = element_text(size=30), 
        #       legend.text = element_text(size=25),
        #       plot.title = element_text(size=40))+
        guides(color = guide_legend(override.aes = list(stroke = 4)))
    
}

pdf(here("plots", "NMF", 
         sprintf("predicted-layers-visium-nmf%sk%s-quant-norm-%s.pdf",
                 label_type, k, slide_number)),
    height=20, width=25)
print(layer_pl[[1]])
print(layer_pl[[2]])
print(layer_pl[[3]])
dev.off()
