library(dplyr)
library(escheR)
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(here)

de.res <- read.csv("processed-data/cindy/ref/TableS8_sig_genes_FDR5perc_enrichment.csv")

n_marker_gene<- 100
gene_df <- de.res |> 
    filter(spatial_domain_resolution == "Sp09") |> 
    group_by(test) |> 
    arrange(fdr, .by_group = TRUE) |> 
    slice_head(n=n_marker_gene)

panel <- readxl::read_excel("processed-data/cindy/ref/XeniumHumanBrainPanelGeneList.xlsx")
gene_df[which(gene_df$gene %in% panel$Gene),] %>% group_by(test) %>% summarize(n=n())

int <- gene_df[which(gene_df$gene %in% panel$Gene),] %>%
    group_by(test) %>%
    mutate(layer = case_when(test == "Sp09D01" | test == "Sp09D02" ~ "L1",
                             test == "Sp09D03" ~ "L2",
                             test == "Sp09D05" ~ "L3",
                             test == "Sp09D08" ~ "L4",
                             test == "Sp09D04" ~ "L5",
                             test == "Sp09D07" ~ "L6",
                             test == "Sp09D06" | test == "Sp09D09" ~ "WM")) 

int <- int %>%
    filter(gene %in% c("MOBP", "AQP4", "CUX2", "RORB", "NTNG1", "HS3ST4", "PVALB"))
slide_number <- "5548"


# read in the sfe to plot the marker genes onto
sfes <- readRDS(here("processed-data", "cindy", sprintf("slide-%s", slide_number), 
                    sprintf("slide%s-all-samples-spe-with-banksy.RDS", slide_number)))
sfe_list <- lapply(unique(sfes$region_id), function(x) 
    sfes[, sfes$region_id == x])

pdf(here('plots', "cindy", "marker_genes", 
        sprintf("slide%s_bayesspace_marker_genes_plots.pdf", slide_number)),
        height=40, width=50)

for (j in 1:length(sfe_list)){
    sfe <- sfe_list[[j]]
    sfe <- scuttle::logNormCounts(sfe)
    for (i in 1:length(unique(int$layer))){
        pls <- list()
        layer <- unique(int$layer)[[i]]
        domain_genes <- unique(int[which(int$layer==layer),]$gene)
        for (k in 1:length(domain_genes)){
            gene <- domain_genes[[k]]
            sfe$marker_logcounts <- logcounts(sfe)[which(rownames(sfe)==gene),]
            p <- make_escheR(sfe, y_reverse=FALSE)%>%
                add_fill("marker_logcounts")+
                ggtitle(paste(gene, layer, sep="_"))+
                theme(plot.title = element_text(size=40), legend.title=element_text(size=30),
                      legend.text = element_text(size=30),
                      legend.key.size = grid::unit(3, "cm"))+
                scale_fill_viridis_c(limits = range(0, 10))+
                guides(guide_legend(size=20))
            pls <- rlist::list.append(pls, p)
        }
        do.call(gridExtra::grid.arrange, c(pls, ncol=3))
    }
}

dev.off()

