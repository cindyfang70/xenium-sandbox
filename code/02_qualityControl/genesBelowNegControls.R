suppressPackageStartupMessages({
    library(Voyager)
    library(SpatialFeatureExperiment)
    library(SingleCellExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(stringr)
    library(Matrix)
    library(BiocParallel)
    library(dplyr)
    library(here)
})

#-------------------------------------------------------------------------------
# Check which genes in the Xenium panel have lower expression than the
# negative control probes
#-------------------------------------------------------------------------------

# read in the data
args <- commandArgs(trailingOnly=TRUE)
sfe <- readRDS(args[[1]])


#sfe <- readRDS(here("processed-data/", "cindy", "slide-5434", "Br6522_Post_SFE.RDS"))
genelist <- readxl::read_excel(here("processed-data","XeniumHumanBrainPanelGeneList.xlsx"))

# the negative control probes are denoted in the colData by is_neg
# for each cell, compute the number of negative control probe counts
# for each gene, compute the average counts per cell and compare it to the 
# average counts per cell for the negative control probe counts

negProbe_means <- rowData(sfe)$means[which(grepl("^NegControlProbe", rownames(rowData(sfe))))]
print(negProbe_means)
negProbe_threshold <- max(negProbe_means)
print(negProbe_threshold)

below_threshold_genes_inds <- which(rowData(sfe)$means < negProbe_threshold & !rowData(sfe)$is_neg) 
below_threshold_genes <- rownames(sfe)[below_threshold_genes_inds]
below_threshold_genes_annot <- genelist[which(genelist$Gene %in% below_threshold_genes),]

print(below_threshold_genes)

plist=list()
for (i in 1:nrow(below_threshold_genes_annot)){
    p <- plotSpatialFeature(sfe, below_threshold_genes_annot[[i,1]], colGeometryName = "cellSeg",
                       exprs_values="counts")+
        ggtitle(paste(below_threshold_genes_annot[i,1], below_threshold_genes_annot[i,3], sep="-"))
    plist[[i]] <- p
}

pdfname <- paste0(sfe$region_id[[1]], "-belowThresholdGenes.pdf")
pdf(here("plots", "cindy", "02_qualityControl", pdfname), height=25, width=20)
do.call(gridExtra::grid.arrange, c(plist, ncol=2))
dev.off()
