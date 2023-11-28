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
# Try using RcppML::nmf for non-negative matrix factorization to compute a feature
# factor matrix (Amplitude matrix in scCoGAPS terminology) from the manually
# annotated visium dataset. Then, use projectR to project the Xenium datasets
# into the column space of the feature factor matrix. 

# Use RcppML::nmf to perform non-negative matrix factorization on the manually
# annotated Visium dataset to get a loadings matrix and a factors matrix with 
# k factors. Plot each of those k factors as well as the layer labels to see
# if any of them are associated with a particular layer.

# get the manually annotated visium data
ehub <- ExperimentHub::ExperimentHub()
vis_anno <- fetch_data(type = "spe", eh = ehub)
# use layer_guess_reordered as the manual annotations

# NMF on the annotated visium data
k <- 20
A <- logcounts(vis_anno) # using logcounts because there are multiple datasets
model <- RcppML::nmf(A, k = k, seed=1237)
patterns <- t(model$h) # these are the factors

rownames(model$w) <- rowData(vis_anno)$gene_name

saveRDS(model, here("processed-data", "cindy", "NMF", sprintf("visium-nmf-model-k%s", k)))

colnames(patterns) <- paste0("NMF", 1:dim(patterns)[[2]])
colData(vis_anno) <- cbind(colData(vis_anno), patterns)

brains <- unique(vis_anno$sample_id)
pdf(here("plots", "NMF", sprintf("all-visium-samples-NMF-k%s.pdf",k)), 
    height=25, width=20)
for(i in 1:k){
    pls<- list()
    patternName <- paste0("NMF", i)
    for (j in seq_along(brains)){
        spe <- vis_anno[,which(vis_anno$sample_id==brains[[j]])]
        pls[[j]] <- make_escheR(spe) |>
            add_ground("layer_guess_reordered") |>
            add_fill(patternName) +
            scale_fill_gradient(low = "white", high = "black")+
            ggtitle(paste(unique(spe$sample_id), patternName))
    }
    cor.mat <- colData(vis_anno)[,c(paste0("NMF", i), 
                               "layer_guess_reordered", "sample_id")] %>%
        as.data.frame()%>%
        pivot_longer(cols=starts_with("NMF"))
    bp <- ggplot(transform(cor.mat, name=factor(name, levels=colnames(patterns))),
                 aes(x=layer_guess_reordered, y=value, 
                     fill=layer_guess_reordered))+
        geom_boxplot()+
        facet_wrap(~sample_id, scales="free")+
        scale_y_log10()+
        theme(legend.position="none", plot.title=element_text(size=20))+
        theme_minimal()
 
    
    do.call(gridExtra::grid.arrange, c(pls, ncol=3))
    print(bp)
    rm(pls)
    gc()
}
dev.off()


    
# https://rdrr.io/bioc/CoGAPS/man/calcCoGAPSStat-methods.html

# 
# recon <- model$w %*% model$h
# 
# #colData(x)<-cbind(colData(x),patterns)
# 
# rownames(model$w) <- rowData(vis_anno)$gene_name
# 
# proj<-projectR(
#     data=as.matrix(counts(x)),
#     loadings=as.matrix(model$w),
#     full = FALSE,
#     family = "gaussianff"
# )
# 
# dim(proj)
