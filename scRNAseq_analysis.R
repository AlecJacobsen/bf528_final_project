rm(list = ls())
install.packages('R.filesets')
library(R.filesets)
library(tidyverse)
library(Seurat)

#Loading in data
cells <- loadRDS("/projectnb/bf528/project_4_scrnaseq/GSM2230760_seurat.rda")

#Finding marker genes
cell_markers <- FindAllMarkers(cells, only.pos = T, min.pct = 0.25)
cell_markers_10 <- cell_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#loading in authors marker genes
author_genes <- read.csv('/projectnb/bf528/users/lava_lamp/project_4/daisy_pr4/3_analyst/MarkerGenes.csv')
#alpha = cluster 2, 4, 8
#beta = clusters 1, 6
#delta = cluster 3
#gamma = cluster 0
#Epsilon = none
#ductal = clusters 5, 10
#acinar = clusters 11
#stellate = cluster 9
#vascular = none
#macrophage = none
#Cytotoxic T = none
#mast = none


#visualizing clusters
pdf('/projectnb/bf528/users/lava_lamp/project_4/Alec_final2/clusters_umap.pdf', width = 10, height = 7)
DimPlot(cells,pt.size = 1, label.size = 5, reduction = "umap")
dev.off()

#cluster heat map
pdf('/projectnb/bf528/users/lava_lamp/project_4/Alec_final2/markergene_heatmap.pdf', width = 17, height = 15)
DoHeatmap(cells, features = cell_markers_10$gene) + NoLegend()
dev.off()

#assigning id to clusters
new.cluster.ids <- c('gamma','beta','alpha','delta','alpha','ductal','beta','unidentified','alpha','stellate','ductal','acinar','unidentified')
names(new.cluster.ids) <- levels(cells)
cells <- RenameIdents(cells, new.cluster.ids)
pdf('/projectnb/bf528/users/lava_lamp/project_4/Alec_final2/cells_umap.pdf', width = 10, height = 7)
DimPlot(cells, reduction = "umap", label = TRUE, pt.size = 1, label.size = 5) + NoLegend()
dev.off()

#new heatmap
pdf('/projectnb/bf528/users/lava_lamp/project_4/Alec_final2/cells_heatmap.pdf', width = 17, height = 15)
DoHeatmap(cells, features = cell_markers_10$gene, angle = 90) + NoLegend()
dev.off()


# re-identifying clusters
panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz") %>%  filter(str_detect(species,"Hs")) %>% select(c("official gene symbol", "cell type", "organ"))

#programmatically identifying cell types
cell_ids <- df <- data.frame(matrix(ncol = 3, nrow = 13))
for (i in 0:12) {
  genes <- cell_markers[cell_markers$cluster == i,7]
  matches <- c()
  for (j in 1:length(genes)){
    matches <- c(matches, (panglao[panglao$`official gene symbol` == genes[j],2]%>% pull(`cell type`)) )
  }
  cell_ids[i+1,] <- names(sort(table(matches), decreasing = TRUE)[1:3])
}

newer_cell_ids <- c('alpha','beta','stellate','acinar','alpha','ductal','beta','platelets','alpha','stellate','ductal','acinar','macrophages')
#renaming clusters
names(newer_cell_ids) <- levels(cells)
cells <- RenameIdents(cells, newer_cell_ids)
pdf('/projectnb/bf528/users/lava_lamp/project_4/Alec_final2/better_cells_umap.pdf', width = 10, height = 7)
DimPlot(cells, reduction = "umap", label = TRUE, pt.size = 1, label.size = 5) + NoLegend()
dev.off()

#new heatmap
pdf('/projectnb/bf528/users/lava_lamp/project_4/Alec_final2/better_cells_heatmap.pdf', width = 17, height = 15)
DoHeatmap(cells, features = cell_markers_10$gene, angle = 90) + NoLegend()
dev.off()

#new marker genes
new_cell_markers <- FindAllMarkers(cells, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
new_cell_markers_3 <- new_cell_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
write.csv(new_cell_markers_3,'/projectnb/bf528/users/lava_lamp/project_4/Alec_final2/new_markers.csv', row.names = F)
