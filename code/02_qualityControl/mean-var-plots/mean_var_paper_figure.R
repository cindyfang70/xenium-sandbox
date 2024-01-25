suppressPackageStartupMessages({
    library(SpatialFeatureExperiment)
})

source(here("code/01_createSCE/xenium_helpers.R"))

sfe <- readRDS(here("processed-data/cindy/slide-5434/xenium-0005434-SFE.RDS"))
sfe <- sfe[,sfe$region_id != "Sample"]

plist <- list()
# make the mean-variance plot individually for each tissue region
for (i in 1:length(unique(sfe$region_id))){
    region <- unique(sfe$region_id)[[i]]
    sub.sfe <- sfe[,sfe$region_id==region]
    print(sub.sfe)
    print(head(rowData(sub.sfe)))
    
    # recompute row means and variances within each tissue region
    rowData(sub.sfe)$means <- rowMeans(counts(sub.sfe))
    rowData(sub.sfe)$vars <- rowVars(counts(sub.sfe))
    
    # create the plot titles
    slide_number <- strsplit(unique(sub.sfe$region_id), split="_")[[1]][[3]]
    donor <- strsplit(unique(sub.sfe$region_id), split="_")[[1]][[1]]
    title <- sprintf("Sample %s on Xenium Slide %s", donor, slide_number )
    p <- plotMeanVar(mean_emp = rowData(sub.sfe)$means, var_emp = rowData(sub.sfe)$vars,
                     plotTitle=title)
    plist[[i]] <- p
}

do.call(gridExtra::grid.arrange, c(plist, ncol=3))

