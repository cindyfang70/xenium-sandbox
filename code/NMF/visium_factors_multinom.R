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
library(nnet)

# Use a multinomial GLM to fit a layer ~ NMF_factors model in order to see which
# combinations of factors correspond to a layer in the manually annotated visium
# dataset. See example here from Yi Wang: 
# https://www.dropbox.com/scl/fi/4m2y8flc1d06z6prkjnjy/shared_example_glm.r?rlkey=34t4q0ruyjq42uldmnlmuxk0n&dl=0

# read in the model
k <- 20 
model <- readRDS(here("processed-data", "cindy", "NMF", sprintf("visium-nmf-model-k%s", k)))
factors <- t(model$h)
colnames(factors) <- paste0("NMF", 1:dim(factors)[[2]])

# get the manually annotated Visium data
ehub <- ExperimentHub::ExperimentHub()
vis_anno <- fetch_data(type = "spe", eh = ehub)
vis_anno <- vis_anno[,which(vis_anno$sample_id == 151673)]

# create the design matrix
design <- cbind(colData(vis_anno)["layer_guess_reordered"], factors)

# fit the multinomial model
mod <-  multinom(layer_guess_reordered ~ ., data = design)

# predict on the same data
p.fit <- predict(mod, predictors=design, type='probs') 

pred = unlist(lapply(1:nrow(p.fit), function(xx){
    colnames(p.fit)[which.max(p.fit[xx,])]
}))

# compute the prediction accuracy

labs <- design$layer_guess_reordered[!is.na(design$layer_guess_reordered)]
acc = mean(pred==labs,na.rm=T)
acc
                 