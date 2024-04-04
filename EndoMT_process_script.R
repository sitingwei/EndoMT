library(tidyverse)
library(Seurat)
library(data.table)
library(grid)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)
library(clustree)
library(ggh4x)
library(ggrepel)
library(SingleR)
library(celldex)

setwd("/Users/swei/Documents/project_PHD/CEA_snRNA-seq")

# 1. read ---------
floders <- list.dirs("/Users/swei/Documents/project_PHD/CEA_snRNA-seq/matrix_and_report") %>%
  str_subset(., pattern = "matrix_10X$") 
samp <- c("patient-1","patient-2","patient-3","patient-4","patient-5","patient-6")
AE_list <- lapply(1:length(floders),function(x){
  tmp <- CreateSeuratObject(counts = Read10X(floders[x]),project = samp[x],min.cells = 3, min.features = 200)
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-") 
  return(tmp)
})
names(AE_list) <- samp

# 2. integration with CCA ------------
## 2.1 filter ---------
AE_list_f <- lapply(AE_list, function(tmp){
  subset(tmp, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 15)
})

## * data for supfigure -------
AE_merge <- merge(AE_list[[1]],y=c(AE_list[[2]],AE_list[[3]],AE_list[[4]],AE_list[[5]],AE_list[[6]]))
AE_merge_1 <- merge(AE_list_f[[1]],y=c(AE_list_f[[2]],AE_list_f[[3]],AE_list_f[[4]],AE_list_f[[5]],AE_list_f[[6]]))
save(AE_merge,AE_merge_1,file = "process_data/sf1.D1.QC_merge_6samples_data_bef_and_af_filter.Rda")

## 2.2 normalization ----------
AE_list_f <- lapply(AE_list_f, function(x){
  tmp <- NormalizeData(x)
  tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
  return(tmp)
})

## 2.3 integrate the data ---------
features <- SelectIntegrationFeatures(object.list = AE_list_f)
AE.anchors <- FindIntegrationAnchors(object.list = AE_list_f, anchor.features = features)
AE_integrate <- IntegrateData(anchorset = AE.anchors)

save(AE_integrate,file="process_data/D1.AE_integrate_with_CAA_no_annotation_6_samples.Rda")

DefaultAssay(AE_integrate) <- "integrated"

## 2.4 reduction and cluter --------
AE_integrate <- ScaleData(AE_integrate, features = rownames(AE_integrate))

AE_integrate <- RunPCA(AE_integrate,npcs = 30,verbose = F)
AE_integrate <- RunUMAP(AE_integrate, reduction = "pca", dims = 1:30)
AE_integrate <- RunTSNE(AE_integrate, reduction = "pca", dims = 1:30)
AE_integrate <- FindNeighbors(AE_integrate, reduction = "pca", dims = 1:30)

### resolution selecting
sce <- AE_integrate
sce <- FindClusters(object = sce,resolution = c(seq(0.2,0.8,0.1)))
clustree(sce@meta.data,prefix = "integrated_snn_res.")

AE_integrate <- FindClusters(AE_integrate, resolution = 0.5)

saveRDS(AE_integrate,file = "process_data/Integrated_cleaned_no_annotation_seurat_of_6_samples_snRNA-seq.RDS")

## * data for clutering ---------
save(AE_integrate,file = "process_data/sf1.D2.integrated_data_with_clustering_and_no_annotation.Rds")

# 3. annotation ----------
celltype <- list(
  SMC <- c("ACTA2"),
  T_cells <- c("CD3E","CD4","CD8A","CD3G","CD3D"),
  Myeloid_cells <- c("CD68","CD14","CD163"),
  EC <- c("CD34","VWF"),
  B_cells <- c("CD79A","CD79B"),
  Mast_cells <- c("KIT","CD63"),
  #immune <- c("PTPRC"),
  CD4T_NK <- c("CD7","IL7R")
) 
names(celltype) <- c("SMC","T_cell","Myeloid_cell","EC","B_cell","Mast_cell","CD4T/NK")
DefaultAssay(AE_integrate) <- "RNA"
## * data for clutering annotation -----------
save(AE_integrate,celltype,file = "process_data/sf1.D3.integrated_data_with_clustering_and_annotation.Rds")

# 4. select subcluster of SMC&EC -------------
DefaultAssay(AE_integrate) <- "integrated"
sub_cluster <- subset(AE_integrate,idents = c(0,1,4,5,7,8,10,12,3))

sub_cluster@meta.data <- sub_cluster@meta.data %>% 
  mutate(cluster = case_when(seurat_clusters %in% c(0,1,4,5,7,8,10,12) ~ "SMC",
                             seurat_clusters == 3 ~ "EC"))

sub_cluster_1 <- ScaleData(sub_cluster,features = rownames(sub_cluster))
sub_cluster_1 <- RunPCA(sub_cluster_1,npcs=30)
sub_cluster_1 <- FindNeighbors(sub_cluster_1,dims=1:30)

### resolution choice
sce_1 <- sub_cluster_1
sce_1@meta.data <- dplyr::select(sce_1@meta.data,orig.ident:percent.mt,cluster)
sce_1 <- FindClusters(object = sce_1, resolution = c(seq(0, 1, 0.1)))
pdf(file = "plot_out/D0326.clustree_for_resolution_choice_with_SMC&EC.pdf", height = 6, width = 5)
clustree(sce_1@meta.data, prefix = "integrated_snn_res.")
dev.off()

sub_cluster_1 <- FindClusters(object = sub_cluster_1, resolution = 0.8, verbose = FALSE)
sub_cluster_1 <- RunTSNE(sub_cluster_1, dims = 1:30)
sub_cluster_1 <- RunUMAP(sub_cluster_1, dims = 1:30)

sub_cluster_1@meta.data <- sub_cluster_1@meta.data %>%
  mutate(sub_cluster = case_when(integrated_snn_res.0.8%in%c(1,2,3,4,6,7,8) ~ paste0("SMC",integrated_snn_res.0.8),
                                 integrated_snn_res.0.8 == 0 ~ "SMC5",
                                 integrated_snn_res.0.8 ==10 ~ "SMC9",
                                 integrated_snn_res.0.8 ==11 ~ "SMC10",
                                 integrated_snn_res.0.8 ==12 ~ "SMC11",
                                 integrated_snn_res.0.8 == 5 ~ "EC1",
                                 integrated_snn_res.0.8 == 9 ~ "EC2"))

## * data for subcluster clustering -------------
umap.res <- data.frame(cbind(sub_cluster_1@reductions$umap@cell.embeddings,sub_cluster_1@meta.data))
save(sub_cluster_1,umap.res,file = "process_data/mf1.D1.selected_data_of_EC_and_SMC_no_annotation.Rda")

# 5. monocle of SMC&EC with RNA$counts array ---------------------
library(monocle)
expr_matrix <- as(as.matrix(sub_cluster_1@assays$RNA@counts),'sparseMatrix')
p_data <- sub_cluster_1@meta.data #usually including cell type and sample origin
f_data <- data.frame(gene_short_name = rownames(sub_cluster_1@assays$RNA),
                     row.names = rownames(sub_cluster_1@assays$RNA)) #including gene annotation and function

# creat CDS object
pd <- new('AnnotatedDataFrame',data=p_data)
fd <- new('AnnotatedDataFrame',data=f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily = negbinomial.size(),
                      lowerDetectionLimit = 0.5)

# estimate size factor and dipersion
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# filter low quality cells 
cds <- detectGenes(cds,min_expr = 0.1) #fData(cds)

expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10)) #16103

# Trajectory-defined gene selection
## a.Selecting differential genes between subclusters
diff <- differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr = "~integrated_snn_res.0.8", cores = 1)
      # fullModeFormulaStr is the variable for doing the differential analysis, can be the column in p_data
deg <- subset(diff, qval < 0.01) 
deg <- deg[order(deg$qval, decreasing = F)[1:2000], ]

ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene) #stored in cds@featureData@data[["use_for_ordering"]]

#reduction
cds <- reduceDimension(cds,max_components = 2,method="DDRTree")

#constructing pseudotime trajectories
ordercds <- orderCells(cds)
#orcercds <- orderCells(ordercds, root_state = 6) #changing the root point

## * data of pseudotime -----------
save(ordercds,file = "process_data/mf1.D2.pseudotime_data_of_SMC&EC.Rda")

# 6. calculate module score ------------
geneset <- fread("/Users/swei/Documents/project_PHD/CEA_snRNA-seq/other_data/Copy of 73 (76) midstage genes_10323_EndMTpaper.csv") %>%
  .$GeneID %>% list()
DefaultAssay(sub_cluster_1) <- "RNA"
mod_subcluster <- AddModuleScore(sub_cluster_1,features = geneset,name="midstage_genes")

mod_score <- mod_subcluster@meta.data %>% 
  dplyr::select(midstage_genes1)
mod_score[mod_score<(-0.1)] <- (-0.1)
mod_score[mod_score>0.3] <- 0.3
meta_dat <- mod_score %>%
  rownames_to_column(var="cell") %>%
  inner_join(pData(ordercds) %>% rownames_to_column(var="cell"),.,by="cell") %>%
  column_to_rownames(var="cell")
pData(ordercds) <- meta_dat 

## * data of pseudotime of module genes -----------
save(ordercds,mod_subcluster,file = "process_data/mf1.D3.pseudotime_data_of_module_score_of_midgenes.Rda")
