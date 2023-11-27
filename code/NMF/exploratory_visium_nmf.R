library(RcppML)
library(Matrix)
library(Voyager)
library(CoGAPS)
library(projectR)
library(spatialLIBD)
library(escheR)

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
model <- RcppML::nmf(counts(vis_anno), k = 10)
patterns <- t(model$h)

colData(vis_anno) <- cbind(colData(vis_anno), patterns)

# plot patterns
pl <- make_escheR(vis_anno) |>
    add_ground("layer_guess_reordered")

p1 <- make_escheR(vis_anno) |>
    add_fill("V1") +
    scale_fill_gradient(low = "white", high = "black")

p2 <- make_escheR(vis_anno) |>
    add_fill("V2") +
    scale_fill_gradient(low = "white", high = "black")

p3 <- make_escheR(vis_anno) |>
    add_fill("V3") +
    scale_fill_gradient(low = "white", high = "black")

p8 <- make_escheR(vis_anno) |>
    add_ground("layer_guess_reordered") |>
    add_fill("V8") +
    scale_fill_gradient(low = "white", high = "black")

pl + p1 + p2
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
