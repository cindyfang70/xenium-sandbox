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

################################################################################
# Annotate the banksy domains using NMF transferred domain labels.
################################################################################
source(here("code", "01_createSCE/", "xenium_helpers.R"))
banksy_sfe <- readRDS(here('processed-data', "cindy", "all-tissues-sfe-with-Banksy-spatial-domains.RDS"))

sfe_5434_nmf <- readRDS(here("processed-data", "cindy", "slide-5434","slide-5434-spe-with-nmf-k200-bayesspace.RDS"))
sfe_5548_nmf <- readRDS(here("processed-data", "cindy", "slide-5548", "slide-5548-spe-with-nmf-k200-bayesspace.RDS"))

sfe_5434_nmf <- SFEtoSPE(sfe_5434_nmf)
sfe_5548_nmf <- SFEtoSPE(sfe_5548_nmf)

# keep only the colData columns that are in both
colData(sfe_5548_nmf) <- colData(sfe_5548_nmf)[colnames(colData(sfe_5434_nmf))]

nmf_sfe <- cbind(sfe_5434_nmf, sfe_5548_nmf)

# make unique cell ids
nmf_sfe$cell_sample_id <- paste(nmf_sfe$cell_id, nmf_sfe$region_id, sep="-")
banksy_sfe$cell_sample_id <- paste(banksy_sfe$cell_id, banksy_sfe$region_id, sep="-")

# make the ordering of cells match in the two separate spe objects
nmf_sfe <- nmf_sfe[,match(banksy_sfe$cell_sample_id, nmf_sfe$cell_sample_id)]
sum(nmf_sfe$cell_sample_id == banksy_sfe$cell_sample_id) # sanity check

nmf_sfe$banksy_domains <- banksy_sfe$clust_M1_lam0.9_k50_res0.45

# proportion of each nmf label in each banksy domain
df <- colData(nmf_sfe)[,c("banksy_domains", "preds_from_bayesspace_NMF_k200")]
prop.table(table(df), 1)

nmf_sfe$nmf_preds <- nmf_sfe$preds_from_bayesspace_NMF_k200

library(RColorBrewer)
myPal <- brewer.pal(12, "Set3")

myPal <- c(myPal, "#83b9c9","#71ab91")
names(myPal) <- levels(colData(nmf_sfe)[["banksy_domains"]])

pdf(here("plots", "cindy", "NMF", "bayesspace", "annotate_banksy_bayesspace_0.45.pdf"),
    height=50, width=50)
for (i in 1:length(unique(nmf_sfe$region_id))){
    sfe <- nmf_sfe[,which(nmf_sfe$region_id == unique(nmf_sfe$region_id)[[i]])]
    
    pls <- list()
    for (j in 1:length(unique(sfe$nmf_preds))){
        pred_label <- unique(sfe$nmf_preds)[[j]]
        
        sfe$is_layer <- sfe$nmf_preds == pred_label
        
        p <- make_escheR(sfe, y_reverse=FALSE) %>%
            add_fill("is_layer", point_size=3) %>%
            add_ground("banksy_domains")+
            scale_fill_manual(values=c("TRUE"="black", "FALSE"="white"))+
            scale_color_manual(values=myPal)+
            ggtitle(paste0(unique(sfe$region_id), pred_label))
        
        pls[[j]] <- p
        
    }
    do.call(gridExtra::grid.arrange, c(pls, ncol=2))
}
dev.off()

props <- prop.table(table(df), 1)
banksy_annot <- unlist(lapply(1:nrow(props), function(xx){
    colnames(props)[which.max(props[xx,])]
})) 

layer_prop <- unlist(lapply(1:nrow(props), function(xx){
    props[xx,which.max(props[xx,])]
})) 

banksy_annot <- as.data.frame(cbind(banksy_clust=rownames(props), banksy_annot, layer_prop))

library(DescTools)

inds <- which(banksy_annot$banksy_annot=="L1")
banksy_clust_with_same_label <- banksy_annot$banksy_clust[inds]

banksy_clust_freqs <- table(df$banksy_domains)
# merge the clusters that get the same annotation
banksy_merged_clust_freqs<- banksy_clust_freqs
banksy_merged_clust_freqs[[inds[[1]]]] <- banksy_clust_freqs[[inds[[1]]]] + banksy_clust_freqs[[inds[[2]]]]
banksy_merged_clust_freqs <- banksy_merged_clust_freqs[-inds[[2]]]

# compare gini coefficient for merged and unmerged
unmerged_probs <- banksy_clust_freqs/sum(banksy_clust_freqs)
gini_unmerged <- 1 - sum(unmerged_probs^2)

merged_probs <- banksy_merged_clust_freqs/sum(banksy_merged_clust_freqs)
gini_merged <- 1 - sum(merged_probs^2)

# 