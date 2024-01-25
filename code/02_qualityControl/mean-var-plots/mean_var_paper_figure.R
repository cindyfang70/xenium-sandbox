suppressPackageStartupMessages({
    library(SpatialFeatureExperiment)
})

source(here("code/01_createSCE/xenium_helpers.R"))

sfe <- readRDS(here("processed-data/cindy/slide-5434/xenium-0005434-SFE.RDS"))
sfe <- sfe[,sfe$region_id != "Sample"]

sfe2 <- readRDS(here("processed-data/cindy/slide-5548/xenium-0005548-SFE.RDS"))
sfe2 <- sfe2[,sfe2$region_id != "Sample"]

all.sfes <- list()

for (i in 1:length(unique(sfe$region_id))){
    region <- unique(sfe$region_id)[[i]]
    sub.sfe <- sfe[,sfe$region_id==region]
    all.sfes[[i]] <- sub.sfe
    }

for (i in 1:length(unique(sfe2$region_id))){
    region <- unique(sfe2$region_id)[[i]]
    sub.sfe <- sfe2[,sfe2$region_id==region]
    all.sfes <- rlist::list.append(all.sfes, sub.sfe)
    }

all_mean_vars <- list()
for (i in 1:length(all.sfes)){
    sub.sfe <- all.sfes[[i]]
    print(sub.sfe)
    print(head(rowData(sub.sfe)))
    
    # recompute row means and variances within each tissue region
    rowData(sub.sfe)$means <- rowMeans(counts(sub.sfe))
    rowData(sub.sfe)$vars <- rowVars(counts(sub.sfe))
    
    # create the plot titles
    slide_number <- strsplit(unique(sub.sfe$region_id), split="_")[[1]][[3]]
    donor <- strsplit(unique(sub.sfe$region_id), split="_")[[1]][[1]]
    title <- sprintf("Sample %s on Xenium Slide %s", donor, slide_number )
    
    mean_emp <- rowData(sub.sfe)$means
    var_emp <- rowData(sub.sfe)$vars
    # Plot
    model = lm(var_emp ~ 1*mean_emp + I(mean_emp^2) + 0, tibble(mean_emp, var_emp))
    phi = 1/coef(model)["I(mean_emp^2)"]
    
    
    mean_var_tb <- tibble(mean_emp = mean_emp,
                          var_emp = var_emp,
                          nbinomial = mean_emp + mean_emp^2 * 1/phi,
                          Poisson = mean_emp, 
                          sample=unique(sub.sfe$region_id))
    
    all_mean_vars[[i]] <- mean_var_tb

}

mean_var_tb <- do.call(rbind, all_mean_vars)


mean_var_tb <- mean_var_tb %>% dplyr::select(mean_emp, var_emp, nbinomial, Poisson, sample) %>% 
    tidyr::pivot_longer(cols=-c(mean_emp, sample), names_to = "model", values_to = "var_value")%>%
    rename(Distribution=model)%>%
    mutate(Distribution=case_when(Distribution=="nbinomial" ~ "Negative binomial",
                                  .default=as.character(Distribution)))

print(mean_var_tb)
mean_var_tb$sample <- as.factor(mean_var_tb$sample)

p <-ggplot(mean_var_tb %>% filter(Distribution %in% c("var_emp")), 
           aes(x = mean_emp, y = var_value)) + 
    geom_point(alpha = 0.3) + 
    geom_line(data = mean_var_tb  %>%
                  filter(Distribution %in% c("Negative binomial", "Poisson")),
              aes(x = mean_emp, y = var_value, colour=Distribution), size=1.25) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Log of mean expression",
         y = "Log of variance") +
    #ggtitle("Mean-variance relationship in Xenium data")+
    facet_wrap(~factor(sample, c(unique(mean_var_tb$sample))))+
    theme_bw()+
    theme(strip.background = element_blank(),
          legend.position="bottom",
          strip.text=element_text(size=25),
          title = element_text(size=30),
          axis.title = element_text(size=25),
          legend.text = element_text(size=25),
          legend.key.size = unit(1, "cm"))#+
    #guides(colour = guide_legend(override.aes = list(size=5)))

pdf(here("plots", "cindy", "mean_var_plots", "xenium_mean_var_plots.pdf"),
    height=15, width=25)
print(p)
dev.off()

#do.call(gridExtra::grid.arrange, c(plist, ncol=3))

