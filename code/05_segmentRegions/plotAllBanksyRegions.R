suppressPackageStartupMessages({
    library(Voyager)
    library(SpatialFeatureExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(stringr)
    library(dplyr)
    library(here)
    library(gridExtra)
    library(escheR)
    library(RColorBrewer)
})
sfe.paths <- commandArgs(trailingOnly=TRUE)
print(sfe.paths)

sfes <- lapply(sfe.paths, readRDS)

plist <- list()
for (i in 1:length(sfes)){
    sfe <- sfes[[i]]
    # plot using escheR 
    colourCount = nlevels(colData(sfe)[["clust_M1_lam0.9_k50_res1.2"]])
    print(colourCount)
    getPalette = colorRampPalette(brewer.pal(12, "Set3"))
    
    p1 <- make_escheR(sfe) %>%
        add_fill(var="clust_M1_lam0.9_k50_res1.2")+
        scale_fill_manual(values=getPalette(colourCount))+
        ggtitle(unique(sfe$region_id))
    plist <- rlist::list.append(plist, p1)
}

fname <- paste("Banksy", "lambda", 
               0.9, "res", 1.2, sep="-")
pdfname <- paste0(fname, ".pdf")
pdf(here("plots", "cindy", "05_segmentRegions", pdfname), width=15, height=20)
do.call(gridExtra::grid.arrange, c(plist, ncol=2))
dev.off()

