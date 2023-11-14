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

# recompute rowwise stats
rowData(sfe)$means <- rowMeans(counts(sfe))
rowData(sfe)$vars <- rowVars(counts(sfe))
rowData(sfe)$cv2 <- rowData(sfe)$vars/rowData(sfe)$means^2

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
    gene <- below_threshold_genes_annot[[i,1]]
    sfe$counts <- counts(sfe)[which(rownames(sfe)==gene),]
    p <- make_escheR(sfe, y_reverse=FALSE)|>
        add_fill(var="counts")+
        scale_fill_continuous()+
        ggtitle(paste(gene, below_threshold_genes_annot[i,3], sep="-"))
    plist[[i]] <- p
}

pdfname <- paste0(sfe$region_id[[1]], "-belowThresholdGenes.pdf")
pdf(here("plots", "cindy", "02_qualityControl", pdfname), height=25, width=20)
do.call(gridExtra::grid.arrange, c(plist, ncol=3))
dev.off()
