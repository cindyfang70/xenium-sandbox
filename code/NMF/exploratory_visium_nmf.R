library(RcppML)
library(Matrix)
library(Voyager)
library(CoGAPS)
library(projectR)
library(spatialLIBD)
library(escheR)
library(here)

# Try using RcppML::nmf for non-negative matrix factorization to compute a feature
# factor matrix (Amplitude matrix in scCoGAPS terminology) from the manually
# annotated visium dataset. Then, use projectR to project the Xenium datasets
# into the column space of the feature factor matrix. 

# get the manually annotated visium data
ehub <- ExperimentHub::ExperimentHub()
vis_anno <- fetch_data(type = "spe", eh = ehub)
vis_anno <- vis_anno[,which(vis_anno$sample_id == 151673)]
# use layer_guess_reordered as the manual annotations

# read in the Xenium data to project
x <- readRDS("processed-data/cindy/slide-5434/Br6471_Post_SFE_filt.RDS")

# NMF on the annotated visium data, try with just counts first but might
# need to use lognormcounts
k <- 10
model <- RcppML::nmf(counts(vis_anno), k = k, seed=1237)
patterns <- t(model$h)

colData(vis_anno) <- cbind(colData(vis_anno), patterns)

plist <- list()
pdf(here("plots", "NMF", sprintf("visium-raw-counts-nmf-k$s.pdf", k)),
    height=20, width=20)
for (i in 1:dim(patterns)[[2]]){
    patternName <- paste0("V", i)
    p <- make_escheR(vis_anno) |>
        add_ground("layer_guess_reordered") |>
        add_fill(patternName) +
        scale_fill_gradient(low = "white", high = "black")+
        theme(legend.position="none")+
        ggtitle(patternName)
    plist[[i]] <- p
    
}
plist <- rlist::list.append(plist, pl)

do.call(gridExtra::grid.arrange, c(plist,ncol=3))
dev.off()
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
