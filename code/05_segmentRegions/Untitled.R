suppressPackageStartupMessages({
    library(Voyager)
    library(SpatialFeatureExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(stringr)
    library(BiocSingular)
    #library(scater)
    library(BiocParallel)
    library(dplyr)
    library(here)
    library(gridExtra)
    library(Banksy)
    library(scuttle)
    library(escheR)
    library(RColorBrewer)
})

#------------------------------------------------------------------------------
# pseudobulking the domains in xenium 
# and seeing how they compare to the domains that you can get from visium. 
# 
# One would hope that similar domains would cluster together in PC space 
# We want to try and quantify any difference we 
# see between the two platforms (assuming an intersecting set of genes)
# Do pseudobulked domains cluster together? do they not?
#------------------------------------------------------------------------------

