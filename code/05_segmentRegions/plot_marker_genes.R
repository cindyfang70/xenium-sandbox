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



sfe <- readRDS("processed-data/cindy/slide-5434/Br6471_Post_SFE_filt.RDS")

pdf(here('plots', "marker_genes", "bayesspace_marker_genes_plots.pdf"),
    height=40, width=50)
for (i in 1:length(unique(int$layer))){
    pls <- list()
    layer <- unique(int$layer)[[i]]
    domain_genes <- int[which(int$layer==layer),]
    for (i in 1:nrow(domain_genes)){
        gene <- domain_genes$gene[[i]]
        sfe$marker_counts <- counts(sfe)[which(rownames(sfe)==gene),]
        p <- make_escheR(sfe, y_reverse=FALSE)%>%
            add_fill("marker_counts")+
            ggtitle(paste(gene, layer, sep="_"))
        pls <- rlist::list.append(pls, p)
    }
    do.call(gridExtra::grid.arrange, c(pls, ncol=3))
}
dev.off()
