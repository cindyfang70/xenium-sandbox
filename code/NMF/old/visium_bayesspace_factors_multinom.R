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
model <- readRDS(here("processed-data", "cindy", "NMF", sprintf("visium-bayesspace-nmf-model-k%s.RDS", k)))
factors <- t(model$h)
colnames(factors) <- paste0("NMF", 1:dim(factors)[[2]])

# get the manually annotated Visium data
ehub <- ExperimentHub::ExperimentHub()
vis_anno <- fetch_data(type = "spatialDLPFC_Visium", eh = ehub)

# create the design matrix
design <- cbind(colData(vis_anno)["BayesSpace_harmony_09"], factors)

# fit the multinomial model
mod <-  multinom(BayesSpace_harmony_09 ~ ., data = design,
                 na.action=na.exclude, maxit=1000, model=TRUE)


# predict on the same data
p.fit <- predict(mod, newdata=design, type='probs') 

pred = unlist(lapply(1:nrow(p.fit), function(xx){
    colnames(p.fit)[which.max(p.fit[xx,])]
}))

# compute the prediction accuracy
labs <- design$BayesSpace_harmony_09
acc = mean(pred==labs,na.rm=TRUE)
acc # k=20: 0.7441152, k=15: 0.712545, k=25: 0.7493769
#k=20 all samples: 0.7703311

# compute the layer-specific prediction accuracy
for (i in 1:length(unique(labs))){
    layer <- unique(labs)[[i]]
    print(sprintf("=====layer %s=====",layer))
    preds <- pred[which(labs==layer)]
    acc <- mean(preds == layer, na.rm=TRUE)
    print(acc)
}

saveRDS(mod, here("processed-data", "cindy", "NMF", sprintf("visium-bayesspace-nmf-k%s-multinom-model.RDS",k)))



# # Plot the correlation between the layer labels and the NMF factors
# layers <- matrix(,length(vis_anno$BayesSpace_harmony_09), ncol=length(unique(vis_anno$BayesSpace_harmony_09)))
# for (i in 1:length(unique(vis_anno$BayesSpace_harmony_09))){
#     layers[,i] <- as.integer(vis_anno$BayesSpace_harmony_09 == 
#         unique(vis_anno$BayesSpace_harmony_09)[[i]])
#     
# }
# colnames(layers) <- paste0('layer', 1:length(unique(vis_anno$BayesSpace_harmony_09)))
# library(corrplot)
# mat <- as.matrix(cbind(layers, factors))
# res <- cor(mat)
# res <- res[grepl("layer", rownames(res)), grepl("NMF", colnames(res))]
# corrplot(res)
