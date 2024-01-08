library(RcppML)
library(Matrix)
library(Voyager)
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
model_type <- args[[1]]
k <- args[[2]]

# get the manually annotated Visium data
ehub <- ExperimentHub::ExperimentHub()
if(model_type=="snRNA-seq"){
    load(here("processed-data", "cindy", "snRNA-seq", "sce_DLPFC.Rdata"))
    layer_name <- "layer_annotation"
}else if(model_type=="manual_annot"){
    sce <- fetch_data(type = "spe", eh = ehub)
    layer_name <- "layer_guess_reordered"
} else if (model_type=="bayesspace"){
    sce <- fetch_data(type = "spatialDLPFC_Visium", eh = ehub)
    layer_name <- "BayesSpace_harmony_09"
}

# Read in the NMF model
nmf.path <- here("processed_data", "cindy", "NMF", model_type,
                 sprintf("%s-nmf-model-k%s.RDS", model_type, k))
nmf.mod <- readRDS(nmf.path)
factors <- t(nmf.mod$h)
colnames(factors) <- paste0("NMF", 1:k)

# Read in the list of factors to use
fcts.path <- here("processed-data", "cindy", "NMF", model_type,
                  sprintf("%s-nmf-selected-factors-k%s.RDS", model_type, k))

fcts.use <- readRDS(fcts.path)

# create the design matrix
design <- cbind(colData(sce)[layer_name], factors)
# rename the layer column so it's consistent across all models
colnames(design)[!grepl("^NMF", colnames(design))] <-"layer"

design <- design[,c("layer", fcts.use)]
print(head(design))

# fit the multinomial model
mod <-  multinom(layer ~ ., data = design,
                 na.action=na.exclude, maxit=1000)

# predict on the same data
p.fit <- predict(mod, predictors=design[grepl("NMF", colnames(design))], type='probs') 

pred = unlist(lapply(1:nrow(p.fit), function(xx){
    colnames(p.fit)[which.max(p.fit[xx,])]
}))

# compute the prediction accuracy
labs <- design$layer[!is.na(design$layer)]
acc = mean(pred==labs,na.rm=TRUE)
acc 

saveRDS(mod, here("processed-data", "cindy", "NMF", model_type, sprintf("%s-nmf-k%s-multinom-model.RDS",model_type, k)))
