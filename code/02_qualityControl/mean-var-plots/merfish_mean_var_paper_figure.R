library(MerfishData)
library(ExperimentHub)
library(ggpubr)

eh <- ExperimentHub()
#query(eh, c("MerfishData", "hypothalamus"))

spe.baysor <- MouseIleumPetukhov2021(segmentation = "baysor")

mean_emp <- rowMeans(counts(spe.baysor))
var_emp <- rowVars(counts(spe.baysor))

model = lm(var_emp ~ 1*mean_emp + I(mean_emp^2) + 0, tibble(mean_emp, var_emp))
phi = 1/coef(model)["I(mean_emp^2)"]


mean_var_tb <- tibble(mean_emp = mean_emp,
                      var_emp = var_emp,
                      nbinomial = mean_emp + mean_emp^2 * 1/phi,
                      Poisson = mean_emp)


mean_var_tb <- mean_var_tb %>% dplyr::select(mean_emp, var_emp, nbinomial, Poisson) %>% 
    tidyr::pivot_longer(cols=-c(mean_emp), names_to = "model", values_to = "var_value")%>%
    rename(Distribution=model)%>%
    mutate(Distribution=case_when(Distribution=="nbinomial" ~ "Negative binomial",
                                  .default=as.character(Distribution)))


p <-ggplot(mean_var_tb %>% filter(Distribution %in% c("var_emp")), 
           aes(x = mean_emp, y = var_value)) + 
    geom_point(alpha = 0.3) + 
    geom_line(data = mean_var_tb  %>%
                  filter(Distribution %in% c("Negative binomial", "Poisson")),
              aes(x = mean_emp, y = var_value, colour=Distribution), size=1.25) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Log of mean expression",
         y = "Log of variance") +
    ggtitle("MERFISH Mouse Ileum")+
    theme_bw()+
    theme(strip.background = element_blank(),
          legend.position="bottom",
          strip.text=element_text(size=25),
          title = element_text(size=30),
          axis.title = element_text(size=25),
          legend.text = element_text(size=25),
          legend.key.size = unit(1, "cm"))#+
#guides(colour = guide_legend(override.aes = list(size=5)))

ggsave(here("plots", "cindy", "mean_var_plots", "merfish_mean_var_plot.pdf"), 
       plot=p,
       height=7.5, width=9)
# print(p)
# dev.off()