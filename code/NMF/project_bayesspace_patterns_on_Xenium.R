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
nmf.mod <- readRDS(here("processed-data", "cindy", "NMF", sprintf("visium-bayesspace-nmf-model-k%s.RDS", k)))
loadings <- nmf.mod$w
vis.factors <- t(nmf.mod$h)

# read in the sfe to project factors onto
sfe <- readRDS(here("processed-data", "cindy", "slide-5434", "slide5434-all-samples-spe-with-banksy.RDS"))

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
    
    fname <- paste0(unique(sfe$region_id), sprintf("-projected-bayesspace-NMF-factors-k%s.RDS", k))
    saveRDS(proj, here("processed-data", "cindy", "NMF", fname))
    
    factors <- t(proj)
    colnames(factors) <- paste0("NMF", 1:k)
    
    fname <- paste0(unique(sfe$region_id), sprintf("-raw-projected-NMF-bayesspace-factors-k%s.RDS", k))
    saveRDS(proj, here("processed-data", "cindy", "NMF", fname))
    
    # rescale the factors to match the Visium factors' range.
    xen.factors.scaled <- matrix(,ncol=k, nrow=nrow(factors))
    for(i in 1:k){
        xen.factors.scaled[,i] <- rescaleFactors(factors[,i], vis.factors[,i])
    }
    colnames(xen.factors.scaled) <- paste0("NMF", 1:k)
    colData(sfe) <- cbind(colData(sfe), xen.factors.scaled)
    
    sfe_list[[j]] <- sfe
}

pdf(here("plots","NMF", sprintf("project-bayesspace-factors-scaled-xenium-k%s-5434.pdf", k)),
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
                         sprintf("visium-bayesspace-nmf-k%s-multinom-model.RDS",k)))
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
        scale_fill_gradient(low="white", high="black")+
        ggtitle(unique(sfe$region_id))+
        theme(legend.title = element_text(size=30), 
              legend.text = element_text(size=25),
              plot.title = element_text(size=40))+
        guides(color = guide_legend(override.aes = list(stroke = 4)))
    
}

pdf(here("plots", "NMF", "predicted-layers-bayesspace-nmf-scaled-5434.pdf"),
    height=15, width=35)
for(i in seq_along(plist)){
    p <-  make_escheR(sfe_list[[i]], y_reverse=FALSE)|>
        add_ground("clust_M0_lam0.9_k50_res0.4")+
        ggtitle(unique(sfe_list[[i]]$region_id))+
        theme(legend.title = element_text(size=30), 
              legend.text = element_text(size=25),
              plot.title = element_text(size=40))+
        guides(color = guide_legend(override.aes = list(stroke = 4)))
    
    p2 <-  make_escheR(sfe_list[[i]], y_reverse=FALSE)|>
        add_ground("predicted_layers")+
        theme(legend.title = element_text(size=30), 
              legend.text = element_text(size=25))+
        guides(color = guide_legend(override.aes = list(stroke = 4)))
    
    #sfe_list[[i]]$prob_less_than_80 <- sfe_list[[i]]$max_probs < 0.8
    
    # p3 <- make_escheR(sfe_list[[i]], y_reverse=FALSE)%>%
    #     add_fill("prob_less_than_80")+
    #     scale_fill_manual(values=c("TRUE"="red", "FALSE"="grey"))
    
    p3 <- make_escheR(sfe_list[[i]], y_reverse=FALSE)|>
        add_ground("clust_M0_lam0.9_k50_res0.4", stroke=3) |>
        add_fill("predicted_layers", point_size=1)+
        theme(legend.title = element_text(size=30), 
              legend.text = element_text(size=25))+
        guides(color = guide_legend(override.aes = list(stroke = 4)),
               fill=guide_legend(override.aes=list(size=4)))
    
    p4 <- make_escheR(sfe_list[[i]], y_reverse=FALSE)%>%
        add_fill("max_probs")+
        scale_fill_viridis_c(limits=c(0, 1))+
        theme(legend.title = element_text(size=30), 
              legend.text = element_text(size=25))+
        ggtitle(unique(sfe_list[[i]]$region_id))
    
    print(ggpubr::ggarrange(plist[[i]], p ,p2, p3,p4, ncol=2))
}
dev.off()

all_props <-list()
plist <- list()
pdf(here("plots", "NMF", "predicted-layers-visium-nmf-bayesspace-proportions-in-banksy-5434.pdf"))
for (i in seq_along(sfe_list)){
    sfe <- sfe_list[[i]]
    #print(sprintf("-------sfe_%s-------", unique(sfe$region_id)))
    
    sfe_props <- list()
    for (k in 1:length(unique(sfe$clust_M0_lam0.9_k50_res0.4))){
        clust <- unique(sfe$clust_M0_lam0.9_k50_res0.4)[[k]]
        clustk <- sfe[,which(sfe$clust_M0_lam0.9_k50_res0.4 == k)]
        #print(sprintf("-------clust%s-------", clust))
        props_k <- round(table(clustk$predicted_layers)/
                             sum(table(clustk$predicted_layers)),3)
        #print(props_k)
        props_k <- as.data.frame(props_k)
        props_k$clust <- clust
        props_k$sfe <- unique(sfe$region_id)
        sfe_props[[k]] <- props_k
        names(sfe_props)[[k]] <- clust
    }
    all_props[[i]] <- do.call(rbind, sfe_props)
    
    p <- ggplot(all_props[[i]], aes(x=Var1, y=Freq, fill=Var1))+
        geom_bar(stat="identity")+
        facet_wrap(~clust, scales="free_y")+
        ylab("Proportion")+
        xlab("Predicted label from NMF")+
        ggtitle(unique(sfe$region_id))+
        cowplot::theme_cowplot()+
        theme(axis.text.x = element_text(angle=90),
              legend.position="none")
    plist[[i]] <- p
    print(p)
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
