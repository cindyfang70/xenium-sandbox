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
################################################################################
# Use RcppML::nmf to perform non-negative matrix factorization on the manually
# annotated Visium dataset to get a loadings matrix and a factors matrix with 
# k factors. Plot each of those k factors as well as the layer labels to see
# if any of them are associated with a particular layer.
################################################################################

args <- commandArgs(trailingOnly = TRUE)
model_type <- args[[1]]

# get the annotated source data
ehub <- ExperimentHub::ExperimentHub()

if(model_type == "manual_annot"){ # use manual annotations
    layer_labs <- "layer_guess_reordered"
    data_type <- "spe"
    vis_anno <- fetch_data(type = data_type, eh = ehub)
}else if(model_type == "bayesspace"){
    layer_labs <- "BayesSpace_harmony_09"
    data_type <- "spatialDLPFC_Visium"
    vis_anno <- fetch_data(type = data_type, eh = ehub)
}else if (model_type=="snRNA-seq"){
    load(here("processed-data", "cindy", "snRNA-seq", "sce_DLPFC.Rdata"))
    vis_anno <- sce # reassign to the same variable so following code doesn't need to be modified
    rm(sce)
    layer_labs <- "layer_annotation"
    vis_anno$sample_id <- vis_anno$Sample
}

# NMF on the annotated visium data
k <- as.numeric(args[[2]])

# trying out GLM-PCA as input to NMF
vis_anno <- scry::nullResiduals(vis_anno, assay="counts", fam="poisson", type="pearson")
vis_anno <- scater::runPCA(vis_anno, ncomponents=50, 
                      ntop = 1000,
                      exprs_values = "poisson_pearson_residuals",
                      scale = TRUE, name = "GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())

A <- reducedDim(vis_anno, "GLM-PCA")

model <- RcppML::nmf(A, k = k, seed=1237, maxit=500)

patterns <- t(model$h) # these are the factors

print(dim(model$w))
print(length(rowData(vis_anno)$gene_name))

rownames(model$w) <- rowData(vis_anno)$gene_name

saveRDS(model, here("processed-data", "cindy", "NMF", "glm-pca", model_type,
                    sprintf("%s-nmf-model-k%s-glm-pca.RDS", model_type, k)))

colnames(patterns) <- paste0("NMF", 1:k)
colData(vis_anno) <- cbind(colData(vis_anno), patterns)

brains <- unique(vis_anno$sample_id)
       add_fill(patternName) +
geom_boxplot()+

# Compute correlation between NMF factors and layer identities 
# colnames(patterns) <- paste0("NMF", 1:k)
patternSums <- colSums(patterns)
patterns <- patterns[,which(patternSums != 0)]

labels <- unique(colData(vis_anno)[[layer_labs]])
n_labels <-length(labels)
factors.ind <- matrix(,nrow=dim(vis_anno)[[2]], ncol=n_labels)
for (i in 1:n_labels){
    label <- labels[[i]]
    factors.ind[,i] <- as.integer(colData(vis_anno)[[layer_labs]]== label)
}

colnames(factors.ind) <- labels
factors.ind <- as.data.frame(factors.ind)

cor.mat <- cbind(patterns, factors.ind)
cor.mat$`NA` <- NULL


M <- cor(cor.mat, use="complete.obs")
M <- M[grepl("NMF", rownames(M)), !grepl("NMF", colnames(M))]

# compute the correlation between NMF factors and the brain IDs
brs <- unique(colData(vis_anno)[["sample_id"]])
n_brs <- length(brs)
brs.ind <- matrix(, nrow=ncol(vis_anno), ncol=n_brs)

for (i in 1:n_brs){
    br <- brs[[i]]
    brs.ind[,i] <- as.integer(colData(vis_anno)[["sample_id"]]==br)
}

colnames(brs.ind) <- brs
brs.ind <- as.data.frame(brs.ind)


brs.cor.mat <- cbind(patterns, brs.ind)

brs.M <- cor(brs.cor.mat, use="complete.obs")
brs.M <- brs.M[grepl("NMF", rownames(brs.M)), !grepl("NMF", colnames(brs.M))]
print(head(brs.M))
print(any(is.na(brs.M)))
print(head(M))
print(any(is.na(M)))

col <- circlize::colorRamp2(seq(-1, 1, length = 3), c("blue", "#EEEEEE", "red"))
fname <- sprintf("%s_NMF_layer_corr_plot_k%s-glm-pca.pdf", model_type, k)
pdf(here("plots", "cindy", "NMF", "glm-pca", model_type, fname), height=30, width=60)
ComplexHeatmap::Heatmap(t(as.matrix(M)), col=col,
                        row_names_gp = grid::gpar(fontsize = 30),
                        column_names_gp = grid::gpar(fontsize = 25),
                        width = ncol(t(as.matrix(M)))*unit(6.5, "mm"), 
                        height = nrow(t(as.matrix(M)))*unit(20, "mm"))
ComplexHeatmap::Heatmap(t(as.matrix(brs.M)), col=col,
                        row_names_gp = grid::gpar(fontsize = 30),
                        column_names_gp = grid::gpar(fontsize = 25),
                        width = ncol(t(as.matrix(brs.M)))*unit(6.5, "mm"), 
                        height = nrow(t(as.matrix(brs.M)))*unit(20, "mm"))
dev.off()

cors <- list(layer_cor=M, sample_cor=brs.M)
saveRDS(cors, here("processed-data", "cindy", "NMF", "glm-pca", model_type,
                   sprintf("%s-nmf-correlations-k%s-glm-pca.RDS", model_type, k)))

