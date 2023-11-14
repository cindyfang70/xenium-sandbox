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
# sfe.paths <- commandArgs(trailingOnly=TRUE)
# print(sfe.paths)
# 
# sfes <- lapply(sfe.paths, readRDS)

sfe.paths1 <- system("ls data/slide-5434", intern=TRUE)
sfe.paths1 <- sfe.paths1[grepl("*filt.RDS", sfe.paths1)]
sfe.paths1 <- paste0("data/slide-5434/", sfe.paths1)

sfe.paths2 <- system("ls data/slide-5548", intern=TRUE)
sfe.paths2 <- sfe.paths2[grepl("*filt.RDS", sfe.paths2)]
sfe.paths2 <- paste0("data/slide-5548/", sfe.paths2)

sfe.paths <- c(sfe.paths1, sfe.paths2)
sfes <- lapply(sfe.paths, readRDS)
plist <- list()
for (i in 1:length(sfes)){
    sfe <- sfes[[i]]
    # plot using escheR 
    colourCount = nlevels(colData(sfe)[["clust_M1_lam0.9_k50_res1.2"]])
    print(colourCount)
    getPalette = colorRampPalette(brewer.pal(12, "Spectral"))
    
    sfe$counts_MOBP <-counts(sfe)[which(rownames(sfe)=="MOBP"),]
    sfe$counts_SST <-counts(sfe)[which(rownames(sfe)=="SST"),]
    p1 <- make_escheR(sfe, y_reverse=FALSE) %>%
        #add_ground(var="clust_M1_lam0.9_k50_res1.2", stroke=0.075)+
        add_fill(var="clust_M1_lam0.9_k50_res1.2", point_size=0.5)+
        scale_fill_discrete()
       # scale_fill_manual(values=getPalette(colourCount))+
        #guides(colour=FALSE)+
        #scale_fill_gradient(low = "white", high = "black")+
        ggtitle(unique(sfe$region_id))
    plist <- rlist::list.append(plist, p1)
}

fname <- paste("Banksy", "lambda", 
               0.9, "res", 1.2, sep="-")
pdfname <- paste0(fname, ".pdf")
#pdf(here("plots", "cindy", "05_segmentRegions", pdfname), width=15, height=20)
pdf(here("plots", "segment", pdfname), width=15, height=20)
do.call(gridExtra::grid.arrange, c(plist, ncol=2))
dev.off()

