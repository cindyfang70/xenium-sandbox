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

# Path on JHPCE: dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC_annotated

args <- commandArgs(trailingOnly = TRUE)

#args <-c(here("processed-data", "cindy", "snRNA-seq", "sce_DLPFC.Rdata"), 2)
spe_path <- args[[1]]
k <- as.numeric(args[[2]])

# try reading in data based on this: https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/80e285b12b54c66126927363c725f57a1591a308/code/05_explore_sce/07_convert_SCE_HDF5.R#L6C1-L10C90
load(spe_path, verbose = TRUE) 

# get the manually annotated visium data
ehub <- ExperimentHub::ExperimentHub()
#sce <- readRDS(spe_path)
# use layer_guess_reordered as the manual annotations

# NMF on the annotated visium data
A <- logcounts(sce) # using logcounts because there are multiple datasets
model <- RcppML::nmf(A, k = k, seed=1237)
patterns <- t(model$h) # these are the factors

rownames(model$w) <- rowData(sce)$gene_name

saveRDS(model, here("processed-data", "cindy", "NMF", sprintf("snRNAseq-nmf-model-k%s", k)))

colnames(patterns) <- paste0("NMF", 1:k)
colData(sce) <- cbind(colData(sce), patterns)

brains <- unique(sce$Sample)
pdf(here("plots", "cindy", "NMF", sprintf("all-snRNAseq-samples-NMF-k%s.pdf",k)), 
    height=25, width=20)
for(i in 1:k){
    pls<- list()
    patternName <- paste0("NMF", i)
    for (j in seq_along(brains)){
        spe <- sce[,which(sce$Sample==brains[[j]])]
        # pls[[j]] <- make_escheR(spe) |>
        #     add_ground("layer_annotation") |>
        #     add_fill(patternName) +
        #     scale_fill_gradient(low = "white", high = "black")+
        #     ggtitle(paste(unique(spe$Sample), patternName))
    }
    cor.mat <- colData(sce)[,c(paste0("NMF", i), 
                                    "layer_annotation", "Sample")] %>%
        as.data.frame()%>%
        pivot_longer(cols=starts_with("NMF"))
    bp <- ggplot(transform(cor.mat, name=factor(name, levels=colnames(patterns))),
                 aes(x=layer_annotation, y=value, 
                     fill=layer_annotation))+
        geom_boxplot()+
        facet_wrap(~Sample, scales="free")+
        scale_y_log10()+
        theme(legend.position="none", plot.title=element_text(size=20))+
        theme_minimal()
    
    # 
    # do.call(gridExtra::grid.arrange, c(pls, ncol=3))
    print(bp)
    #rm(pls)
    gc()
}
dev.off()

