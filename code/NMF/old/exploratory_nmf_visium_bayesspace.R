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
# Use RcppML::nmf to perform non-negative matrix factorization on the manually
# annotated Visium dataset to get a loadings matrix and a factors matrix with 
# k factors. Plot each of those k factors as well as the layer labels to see
# if any of them are associated with a particular layer.

# get the manually annotated visium data
ehub <- ExperimentHub::ExperimentHub()
vis_anno <- fetch_data(type = "spatialDLPFC_Visium", eh = ehub)
# use BayesSpace_harmony_09 as the manual annotations

# NMF on the annotated visium data
k <- 9
A <- logcounts(vis_anno) # using logcounts because there are multiple datasets
model <- RcppML::nmf(A, k = k, seed=1237)
patterns <- t(model$h) # these are the factors

rownames(model$w) <- rowData(vis_anno)$gene_name

saveRDS(model, here("processed-data", "cindy", "NMF", sprintf("visium-bayesspace-nmf-model-k%s.RDS", k)))

colnames(patterns) <- paste0("NMF", 1:k)
colData(vis_anno) <- cbind(colData(vis_anno), patterns)



brains <- unique(vis_anno$sample_id)
pdf(here("plots", "NMF", sprintf("all-visium-bayesspace-samples-NMF-k%s.pdf",k)), 
    height=35, width=35)
for(i in 1:k){
    pls<- list()
    patternName <- paste0("NMF", i)
    for (j in seq_along(brains)){
        spe <- vis_anno[,which(vis_anno$sample_id==brains[[j]])]
        spe$BayesSpace_harmony_09 <- as.factor(spe$BayesSpace_harmony_09)
        pls[[j]] <- make_escheR(spe) |>
            add_ground("BayesSpace_harmony_09") |>
            add_fill(patternName) +
            scale_fill_continuous(limits=c(min(factors), max(factors)))+
            scale_fill_gradient(low = "white", high = "black")+
            ggtitle(paste(unique(spe$sample_id), patternName))
    }
    cor.mat <- colData(vis_anno)[,c(paste0("NMF", i), 
                                    "BayesSpace_harmony_09", "sample_id")] %>%
        as.data.frame()%>%
        pivot_longer(cols=starts_with("NMF"))
    bp <- ggplot(transform(cor.mat, name=factor(name, levels=colnames(patterns))),
                 aes(x=as.factor(BayesSpace_harmony_09), y=value, 
                     fill=as.factor(BayesSpace_harmony_09)))+
        geom_boxplot()+
        facet_wrap(~sample_id, scales="free_x")+
        scale_y_log10()+
        theme_minimal()
        theme(legend.position="none", plot.title=element_text(size=20))
        
    
    
    do.call(gridExtra::grid.arrange, c(pls, ncol=3))
    print(bp)
    rm(pls)
    gc()
}
dev.off()
