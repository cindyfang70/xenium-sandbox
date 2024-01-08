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
# Use the NMF loadings from Visium/snRNA-seq data to project the Xenium data 
# and get factors for Xenium.
################################################################################
args <- c("snRNA-seq", "20", 
          "processed-data/cindy/NMF/snRNA-seq/snRNAseq-nmf-model-k20.RDS", "5548")

# read in the model
model_type <- args[[1]]
k <- as.numeric(args[[2]])
 # need to find the path in the script rather than as input
slide_number <- args[[3]]

# get the manually annotated Visium data
ehub <- ExperimentHub::ExperimentHub()
if(model_type=="snRNA-seq"){
    load(here("processed-data", "cindy", "snRNA-seq", "sce_DLPFC.Rdata"))
    layer_name <- "layer_annotation"
}else if(model_type=="manual_annot"){
    sce <- fetch_data(type = "spe", eh = ehub)
    layer_name <- "layer_guess_reordered"
} else if (model_type=="bayesspace"){
    sce <- fetch_data(type = "spatialDLPFC_Visium", eh = ehub)
    layer_name <- "BayesSpace_harmony_09"
}



# read the model in
nmf.path <- here("processed_data", "cindy", "NMF", model_type,
                 sprintf("%s-nmf-model-k%s.RDS", model_type, k))
nmf.mod <- readRDS(nmf.path)
loadings <- nmf.mod$w
vis.factors <- t(nmf.mod$h)

# read in the sfe to project factors onto
sfe <- readRDS(here("processed-data", "cindy", sprintf("slide-%s", slide_number), 
                    sprintf("slide%s-all-samples-spe-with-banksy.RDS", slide_number)))
sfe_list <- lapply(unique(sfe$region_id), function(x) 
    sfe[, sfe$region_id == x])
# subset the loadings matrix to only those genes that are in the Xenium panel
i<-intersect(rownames(sfe),rownames(loadings))
loadings<-loadings[rownames(loadings) %in% i,]
loadings<-loadings[match(rownames(sfe),rownames(loadings)),]


# rescaleFactors <- function(query_vec, ref_vec){
#     scale_max <- max(ref_vec)
#     scale_min <- min(ref_vec)
#     
#     query_scaled <- scales::rescale(query_vec, to=c(scale_min, scale_max))
#     return(query_scaled)
# }


for (j in seq_along(sfe_list)){
    sfe <- sfe_list[[j]]
    # log normalize and then project using the loadings
    sfe <- scuttle::logNormCounts(sfe)
    proj<-project(A=logcounts(sfe), w=loadings, L1=0)
    
    fname <-  paste0(unique(sfe$region_id), 
                     sprintf("-raw-projected-NMF%sfactors-k%s.RDS", 
                             model_type, k))
    saveRDS(proj, here("processed-data", "cindy", "NMF", model_type, fname))
    
    factors <- t(proj)
    colnames(factors) <- paste0("NMF_k",k,"_", 1:k)
    
    # # rescale the factors to match the Visium factors' range.
    # xen.factors.scaled <- matrix(,ncol=k, nrow=nrow(factors))
    # for(i in 1:k){
    #     xen.factors.scaled[,i] <- rescaleFactors(factors[,i], vis.factors[,i])
    # }
    # colnames(xen.factors.scaled) <- paste0("NMF_k",k,"_", 1:k)
    colData(sfe) <- cbind(colData(sfe), factors)

    sfe_list[[j]] <- sfe
}

if(!file.exists(here("plots","NMF", 
                     sprintf("project%sfactors-scaled-xenium-k%s-%s.pdf", 
                             model_type,k, slide_number)))){
pdf(here("plots","NMF", sprintf("project%sfactors-scaled-xenium-k%s-%s.pdf", 
                                model_type,k, slide_number)),
        height=15, width=45)
    for (i in 1:k){
        pls <- list()
        for (j in seq_along(sfe_list)){
            sfe <- sfe_list[[j]]
            patternName <- paste0("NMF_k",k,"_", i)
            factors <- as.matrix(colData(sfe)
                                 [,grepl(sprintf("NMF_k%s",k), colnames(colData(sfe)))])
            pls[[j]] <- make_escheR(sfe, y_reverse=FALSE) |>
                add_ground("clust_M0_lam0.9_k50_res0.4")|>
                add_fill(patternName) +
                scale_fill_gradient(low="white", high="black", 
                                    limits=c(min(factors), max(factors)))+
                ggtitle(paste(unique(sfe$region_id), patternName))
        }
        
        do.call(gridExtra::grid.arrange, c(pls, ncol=3))
    }
    dev.off()
    
}


## predict using the multinomial model
# multinom <- readRDS(here("processed-data", "cindy", "NMF", 
#                          sprintf("visium-nmf-k%s-multinom-model.RDS",k)))
# library(nnet)
# library(MASS)
# plist <- list()
# for (i in seq_along(sfe_list)){
#     sfe <- sfe_list[[i]]
#     factors <- colData(sfe)[grepl(sprintf("NMF_k%s",k), colnames(colData(sfe)))]
#     colnames(factors) <- paste0("NMF", 1:k)
#     
#     # predict using the multinomial model
#     probs <- predict(multinom, newdata=factors, type='probs', 
#                      na.action=na.exclude)
#     print(unique(probs))
#     preds <- unlist(lapply(1:nrow(probs), function(xx){
#              colnames(probs)[which.max(probs[xx,])]
#         })) 
#     maxprobs <- unlist(lapply(1:nrow(probs), function(xx){
#         max(probs[xx,])
#     })) 
#     
#     #sfe$max_probs <- maxprobs
#     
#     preds_name <- sprintf("predicted_layers_NMF_k%s_manual_annot",k)
#     colData(sfe)[preds_name] <- preds
#     
#     sfe_list[[i]] <- sfe
#     sfe$counts_MOBP <- counts(sfe)[which(rownames(sfe)=="MOBP"),]
#     plist[[i]] <- make_escheR(sfe, y_reverse=FALSE)|>
#         add_ground(preds_name) |>
#         add_fill("counts_MOBP") +
#         scale_fill_gradient(low="white", high="black")+
#         ggtitle(unique(sfe$region_id))+
#         theme(legend.title = element_text(size=30), 
#               legend.text = element_text(size=25),
#               plot.title = element_text(size=40))+
#         guides(color = guide_legend(override.aes = list(stroke = 4)))
# 
# }
# 
# sfe_all <- do.call(cbind, sfe_list)
# saveRDS(sfe_all, here("processed-data", "cindy", sprintf("slide-%s", slide_number), 
#                       sprintf("slide%s-all-samples-spe-with-banksy-NMF.RDS", slide_number)))
# if (!file.exists(here("plots", "NMF", sprintf("predicted-layers-nmf-scaled-%s.pdf", slide_number)))){
# pdf(here("plots", "NMF", sprintf("predicted-layers-nmf-scaled-%s.pdf", slide_number)),
#         height=15, width=35)
#     for(i in seq_along(plist)){
#         p <-  make_escheR(sfe_list[[i]], y_reverse=FALSE)|>
#             add_ground("clust_M0_lam0.9_k50_res0.4")+
#             ggtitle(unique(sfe_list[[i]]$region_id))+
#             theme(legend.title = element_text(size=30), 
#                   legend.text = element_text(size=25),
#                   plot.title = element_text(size=40))+
#             guides(color = guide_legend(override.aes = list(stroke = 4)))
#         
#         p2 <-  make_escheR(sfe_list[[i]], y_reverse=FALSE)|>
#             add_ground(preds_name)+
#             theme(legend.title = element_text(size=30), 
#                   legend.text = element_text(size=25))+
#             guides(color = guide_legend(override.aes = list(stroke = 4)))
#         
#         #sfe_list[[i]]$prob_less_than_80 <- sfe_list[[i]]$max_probs < 0.8
#         
#         p4 <- make_escheR(sfe_list[[i]], y_reverse=FALSE)%>%
#             add_fill("max_probs")+
#             scale_fill_viridis_c(limits=c(0, 1))
#             theme(legend.title = element_text(size=30), 
#                   legend.text = element_text(size=25))+
#             ggtitle(unique(sfe_list[[i]]$region_id))
#         
#         p3 <- make_escheR(sfe_list[[i]], y_reverse=FALSE)|>
#             add_ground("clust_M0_lam0.9_k50_res0.4", stroke=3) |>
#             add_fill(preds_name, point_size=1)+
#             theme(legend.title = element_text(size=30), 
#                   legend.text = element_text(size=25))+
#             guides(color = guide_legend(override.aes = list(stroke = 4)),
#                    fill=guide_legend(override.aes=list(size=4)))
#         # 
#         print(ggpubr::ggarrange(plist[[i]], p ,p2, p3, p4, ncol=2))
#     }
#     dev.off()
#     
# }
# 
# 
# all_props <-list()
# plist <- list()
# pdf(here("plots", "NMF", sprintf("predicted-layers-visium-nmf-proportions-in-banksy-%s.pdf", slide_number)))
# for (i in seq_along(sfe_list)){
#     sfe <- sfe_list[[i]]
#     #print(sprintf("-------sfe_%s-------", unique(sfe$region_id)))f
#     
#     sfe_props <- list()
#     for (k in 1:length(unique(sfe$clust_M0_lam0.9_k50_res0.4))){
#         clust <- unique(sfe$clust_M0_lam0.9_k50_res0.4)[[k]]
#         clustk <- sfe[,which(sfe$clust_M0_lam0.9_k50_res0.4 == k)]
#         #print(sprintf("-------clust%s-------", clust))
#         props_k <- round(table(colData(clustk)[preds_name])/
#                   sum(table(colData(clustk)[preds_name])),3)
#         #print(props_k)
#         props_k <- as.data.frame(props_k)
#         colnames(props_k)[[1]] <- "Layer"
#         props_k$clust <- clust
#         props_k$sfe <- unique(sfe$region_id)
#         sfe_props[[k]] <- props_k
#         names(sfe_props)[[k]] <- clust
#     }
#     all_props[[i]] <- do.call(rbind, sfe_props)
#     
#     
#     p <- ggplot(all_props[[i]], aes(x=Layer, y=Freq, fill=Layer))+
#         geom_bar(stat="identity")+
#         facet_wrap(~clust, scales="free_y")+
#         ylab("Proportion")+
#         xlab("Predicted label from NMF")+
#         ggtitle(unique(sfe$region_id))+
#         cowplot::theme_cowplot()+
#         theme(axis.text.x = element_text(angle=90),
#               legend.position="none")
#     plist[[i]] <- p
#     print(p)
# }
# dev.off()

# all_props_df <- do.call(rbind, all_props)
# 
# ggplot(all_props_df, aes(x=Var1, y=Freq, fill=Var1))+
#     geom_bar(stat="identity")+
#     facet_wrap(~sfe + clust, scales="free")+
#     theme(axis.text.x = element_text(angle=90))
