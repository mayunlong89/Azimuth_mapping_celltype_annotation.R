library(tidyverse)
library(SeuratObject)
library(anndata)
library(SparseM)

dataset_prefix = "CMC"
data_path = "/path/to/all/files/"

#### Read in data from Ma et al. into a Seurat object
rawRDS_Ma <- readRDS(paset0(data_path,"Ma_Sestan_mat.rds"))
seurat_obj_ref_Ma <- rawRDS_Ma

#### Read in data from BICCN into a Seurat object
rawRDS_BICCN = readRDS(paset0(data_path,"BICCN_mat.RDS"))
metaRDS_BICCN = readRDS(paset0(data_path,"BICCN_meta_share.RDS"))
metaRDS_BICCN = as.data.frame(metaRDS_BICCN)
rownames(metaRDS_BICCN) <- metaRDS_BICCN$sample_id
seurat_obj_ref_BICCN <- CreateSeuratObject(counts = rawRDS_BICCN, project = "dlpfc", meta.data = metaRDS_BICCN,min.cells=1, min.features = 0)
Idents(seurat_obj_ref_BICCN) <- seurat_obj_ref_BICCN$method

### Read in  using read_h5ad########################

### This file naming convention assumes that the dataset has been first prefiltered by 
#### marker gene annotation (see Supplementary Methods), and that this file 
#### has the dataset name as a prefix.
seurat_obj_query <- read_h5ad(paste0(data_path,dataset_prefix,"_Hybrid_filtered.h5ad"))
print(dim(seurat_obj_query$raw$X))
count.data <- t(seurat_obj_query$raw$X) 

#### This is a conversion process to make the raw count matrix processed by Pegasus is
#### accessible in the Seurat object.
count.data <- as(count.data,"matrix.csr")
count.data <- as(count.data,"dgCMatrix")
colnames(count.data) <- row.names(seurat_obj_query$obs)
row.names(count.data) <- row.names(seurat_obj_query$var)
seurat_obj_query <- SeuratObject::CreateSeuratObject(counts = count.data, meta.data = seurat_obj_query$obs)

#######################################################

#### Preprocess and run PCA
seurat_obj_ref_BICCN <- SCTransform(seurat_obj_ref_BICCN)
seurat_obj_ref_Ma <- SCTransform(seurat_obj_ref_Ma)
seurat_obj_query <- SCTransform(seurat_obj_query)

seurat_obj_ref_BICCN <- RunPCA(seurat_obj_ref_BICCN)
seurat_obj_ref_Ma <- RunPCA(seurat_obj_ref_Ma)
seurat_obj_query <- RunPCA(seurat_obj_query)

#####Transfer labels for Ma_Sestan data
transfer_anchors_Ma <- FindTransferAnchors(reference = seurat_obj_ref_Ma, query = seurat_obj_query, dims = 1:30, reduction = "pcaproject")
predictions_Ma <- TransferData(anchorset = transfer_anchors_Ma, refdata = seurat_obj_ref_Ma$subclass, dims = 1:30)

seurat_obj_query <- AddMetaData(seurat_obj_query, metadata = predictions_Ma)

write.csv(seurat_obj_query$predicted.id, paste0(data_path,dataset_prefix,"_Azimuth_predictions_Ma_Sestan.csv"),quote = F, row.names = T)

seurat_obj_query <- RunUMAP(seurat_obj_query, dims = 1:30)
Idents(seurat_obj_query) <- seurat_obj_query$predicted.id

#### Output a UMAP with the Ma et al. predictions alone
png(file=paste0(data_path,dataset_prefix,"_Azimuth_Transferred_UMAP_Ma_Sestan.png"))
DimPlot(seurat_obj_query, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#################################################
#####Transfer labels for BICCN data
transfer_anchors_BICCN <- FindTransferAnchors(reference = seurat_obj_ref_BICCN, query = seurat_obj_query, dims = 1:30, reduction = "pcaproject")
predictions_BICCN <- TransferData(anchorset = transfer_anchors_BICCN, refdata = seurat_obj_ref_BICCN$within_area_subclass, dims = 1:30)

seurat_obj_query <- AddMetaData(seurat_obj_query, metadata = predictions_BICCN)

write.csv(seurat_obj_query$predicted.id, paste0(data_path,dataset_prefix,"_Azimuth_predictions_BICCN.csv"),quote = F, row.names = T)

seurat_obj_query <- RunUMAP(seurat_obj_query, dims = 1:30)
Idents(seurat_obj_query) <- seurat_obj_query$predicted.id

#### Output a UMAP with the BICCN alone
png(file=paste0(data_path,dataset_prefix,"_Azimuth_Transferred_UMAP_BICCN.png"))
DimPlot(seurat_obj_query, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#################################################

#### The predictions for each stage are reconciled and combined in the next step.

dev.off()
