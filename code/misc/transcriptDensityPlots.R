sfes <- system("ls data/slide-5548", intern=TRUE)
sfes <- paste0("data/slide-5548/", sfes)

sfes2 <-  system("ls data/slide-5434", intern=TRUE)
sfes2 <- paste0("data/slide-5434/", sfes2)

sfe.paths <- c(sfes, sfes2)
plist <- list()
density.plist <- list()
for (i in 1:length(sfe.paths)){
    path <- sfe.paths[[i]]
    sfe <- readRDS(path)
    sfe$log1p_transcript_counts <- log((sfe$transcript_counts+1))
    p <- plotSpatialFeature(sfe, "log1p_transcript_counts", 
                            colGeometryName = "cellSeg")+
        scale_fill_viridis_c()+
        ggtitle(unique(sfe$region_id))
    plist <- rlist::list.append(plist, p)
    
    sfe$transcript_density <- log(sfe$transcript_counts/sfe$cell_area)
    p2 <- plotSpatialFeature(sfe, "transcript_density",
                             colGeometryName="cellSeg")+
        scale_fill_viridis_c()+
        ggtitle(unique(sfe$region_id))
    
    density.plist <- rlist::list.append(density.plist, p2)
}

pdf("transcriptDensity.pdf")
do.call(gridExtra::grid.arrange, plist)
do.call(gridExtra::grid.arrange, density.plist)
dev.off()


sfe$high_counts <- sfe$log1p_transcript_counts > 6.5
plotSpatialFeature(sfe, "high_counts", colGeometryName="cellSeg")
