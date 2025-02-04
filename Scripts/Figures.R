
library(clusterProfiler)
library(here)
library(Seurat)
library(tidyverse)
library(ggprism)
library(enrichplot)
library(msigdbr)
library(ggh4x)
library(RColorBrewer)
library(circlize)


combined <- readRDS(here("Data","Integrated_Hepatocytes.RDS"))

# Helper functions
# For UMAPs
plot_function <- function(scObj,cells,group_by,output_tiff="umap_tiff"){
  
  axis <- ggh4x::guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(4, "cm")
  )
  DimPlot(object = scObj, group.by = group_by,cells=cells,
          label = F, pt.size = 2, cols = c("lightblue", "red"),order = T)+
    NoLegend()+
    labs(x="UMAP1",y="UMAP2")+
    guides(x = axis, y = axis)+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL)+
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(hjust=0,size=28),
          plot.title=element_blank(),
          axis.line = element_line(arrow = arrow()),
          legend.text = element_text(size=20))+guides(color=guide_legend(override.aes = list(size=8)))
  #ggsave(output_tiff,dpi=300,units='cm',width=30,height=25)
  
}

# For FeaturePlots
featureplot_function <- function(scObj,cells,features,output_tiff = "umap_tiff"){
  
  axis <- ggh4x::guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(4, "cm")
  )
  FeaturePlot(object = scObj, features=features,cells=cells,
              label = F, pt.size = 2, cols = c("lightblue", "red"),order = T)+
    NoLegend()+
    #scale_color_manual(values=colors2)+
    #scale_color_locuszoom()+
    labs(x="UMAP1",y="UMAP2")+
    guides(x = axis, y = axis)+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL)+
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(hjust=0,size=28),
          plot.title=element_blank(),
          axis.line = element_line(arrow = arrow()),
          legend.text = element_text(size=20))+guides(color=guide_legend(override.aes = list(size=8)))
  #ggsave(output_tiff,dpi=300,units='cm',width=30,height=25)
  
}

Idents(combined) <- "fibrosis.stage_simple"
identities <- unique(Idents(combined))
cell_list <- lapply(identities, function(x) {
  WhichCells(combined, idents = x)
})
# Name the list with the identities for easier reference
names(cell_list) <- identities


# Fig 1G
plot_function(combined,cell_list[["0"]],"sen_duke_status")
plot_function(combined,cell_list[["1"]],"sen_duke_status")
plot_function(combined,cell_list[["2"]],"sen_duke_status")
plot_function(combined,cell_list[["3"]],"sen_duke_status")
plot_function(combined,cell_list[["4"]],"sen_duke_status")

# Fig 5G
featureplot_function(combined,cell_list[["0"]],"CDKN1A")
featureplot_function(combined,cell_list[["1"]],"CDKN1A")
featureplot_function(combined,cell_list[["2"]],"CDKN1A")
featureplot_function(combined,cell_list[["3"]],"CDKN1A")
featureplot_function(combined,cell_list[["4"]],"CDKN1A")

# Fig ED 8G
featureplot_function(combined,cell_list[["0"]],"GDF15")
featureplot_function(combined,cell_list[["1"]],"GDF15")
featureplot_function(combined,cell_list[["2"]],"GDF15")
featureplot_function(combined,cell_list[["3"]],"GDF15")
featureplot_function(combined,cell_list[["4"]],"GDF15")

# SHGS DEG
DEG_SHGS <- read.csv(here("Data","DEG_SHGS_cleaned.csv"))

genelist <- DEG_SHGS$avg_log2FC
names(genelist) <- DEG_SHGS$Gene
ranked.gene.list <- sort(genelist, decreasing = TRUE)

# Fig 1 I
gsea_results_reactome <- readRDS(here("Data","gsea_results_reactome_SHGSvsSHGSneg.RDS"))

gseaplot2(gsea_results_reactome,title="REACTOME_CELLULAR_SENESCENCE",
          geneSetID = "REACTOME_CELLULAR_SENESCENCE",
          base_size=24,rel_heights=c(1.5,0.5),
          pvalue_table = F, subplots = 1:2,ES_geom="line")
gseaplot2(gsea_results_reactome,title="REACTOME_SASP",
          geneSetID = "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
          base_size=24,rel_heights=c(1.5,0.5),
          pvalue_table = F, subplots = 1:2,ES_geom="line")

# Fig 1J
gene_set_kegg <- readRDS(here("Data","gsea_results_KEGG_SHGSvsSHGSneg.RDS"))

dotplot(gene_set_kegg, showCategory=10, split=".sign", font.size=12) + 
  scale_size(range = c(3, 10)) +  
  scale_color_gradient(low = "red", high = "blue") +  
  facet_grid(. ~ .sign) + 
  theme(
    strip.text = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5) 
  )         
