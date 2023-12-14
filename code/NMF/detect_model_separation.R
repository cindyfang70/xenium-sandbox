library(brglm2)
library(here)
library(spatialLIBD)

k <- 20
model <- readRDS(here("processed-data", "cindy", "NMF", sprintf("visium-nmf-model-k%s.RDS", k)))
factors <- t(model$h)
colnames(factors) <- paste0("NMF", 1:dim(factors)[[2]])

# get the manually annotated Visium data
ehub <- ExperimentHub::ExperimentHub()
vis_anno <- fetch_data(type = "spe", eh = ehub)

# make the design matrix 
design <- cbind(colData(vis_anno)["layer_guess_reordered"], factors)

mod <- brmultinom(as.factor(design$layer_guess_reordered) ~ factors,
                  ref=1,
                     type = "ML", check_aliasing=FALSE)

saveRDS(mod, here("processed-data", "cindy", "NMF", "visium-brglm2-multinom-k%s.RDS"))
summary(mod)

library(pmlr)
mod <- pmlr(as.factor(design$layer_guess_reordered) ~ factors,
            penalized=TRUE)
