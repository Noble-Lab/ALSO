
####################################################################################################################
# merge samples belonging to every stage (embryo, ranging from E9.5 to E13.5)
####################################################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
library(dplyr)
library(irlba)
library(DelayedArray)
set.seed(1234)

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/cross_Nmasked/data"
source("/net/shendure/vol10/projects/cxqiu/nobackup/work/cross_Nmasked/code/help_code/lsi.R")

peak_matrix = readRDS(paste0(work_path, "/0_merge_data/counts_embryo.rds"))
meta = readRDS(paste0(work_path, "/0_merge_data/pd_embryo.rds"))

amulet_res = read.table("/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/analysis_mouse_amulet/AMULET_res.txt", header=T, as.is=T)
rownames(amulet_res) = as.vector(amulet_res$cell_id)
print(sum(!rownames(meta) %in% rownames(amulet_res)))
amulet_res = amulet_res[rownames(meta),]
meta$amulet_doublets = ifelse(amulet_res$q_val < 0.05, "yes", "no")
print(table(meta$amulet_doublets))

### first removing doublets
meta = meta[meta$amulet_doublets == "no",]
peak_matrix = peak_matrix[,colnames(peak_matrix) %in% rownames(meta)]

### binary matrix and remove low express cells and peaks
### filter out peaks which are detected in less than 0.1% cells
### filter out cells which are deteced less than 100 peaks
peak_matrix@x[peak_matrix@x > 0] = 1
peak_matrix = filter_features(peak_matrix, cells=ncol(peak_matrix) * 0.001)
peak_matrix = filter_cells(peak_matrix, 100)

use_meta = meta[colnames(peak_matrix),]

print(dim(use_meta))
print(dim(peak_matrix))

seurat_obj = lsi_workflow(peak_matrix, 
                          dims=2:50, 
                          metadata=use_meta, 
                          log_scale_tf=TRUE, 
                          reduction='pca.l2', 
                          resolution = 0.6, 
                          n.neighbors=50, 
                          n.components=2, 
                          min.dist=0.01, 
                          k.param=20)

seurat_obj = Seurat::FindClusters(seurat_obj, n.start=20, resolution = 1)

saveRDS(seurat_obj, file = paste0(work_path, "/1_embedding_embryo/obj_peaks.rds"))

pd = data.frame(seurat_obj[[]])
emb = Embeddings(seurat_obj, reduction = "umap")
pd$UMAP_1 = emb[,1]
pd$UMAP_2 = emb[,2]
saveRDS(pd, file = paste0(work_path, "/1_embedding_embryo/pd_peaks.rds"))
