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
library(RColorBrewer)
library(tidydr)
library(MASS)
        
setwd("/Users/swei/Documents/project_PHD/CEA_snRNA-seq")

# sf.QC plot -----------------
## sf1.a. violin plot bef & after filter -------------
load("process_data/sf1.D1.QC_merge_6samples_data_bef_and_af_filter.Rda")
samp_col <- c("#B5D66E","#FBB463","#80B1D3","#BEBDDA","#8DD3C6","#F28072")

p1 <- VlnPlot(AE_merge,features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
              group.by = "orig.ident", pt.size = 0,cols = samp_col,combine = T)
p2 <- VlnPlot(AE_merge_1,features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
              group.by = "orig.ident", pt.size = 0,cols = samp_col,combine = T)
pdf("plot_out/sf1.a.violin_QC_plot_before_filter.pdf",width = 8,height = 5.5)
p1/p2
dev.off()

## sf1.b. cell number bef & after filter --------------
load("process_data/sf1.D1.QC_merge_6samples_data_bef_and_af_filter.Rda")
cell_count <- as.data.frame(table(AE_merge$orig.ident)) %>%
  merge(.,as.data.frame(table(AE_merge_1$orig.ident)),by = "Var1")
colnames(cell_count) <- c("samples","before_filter","after_filter")
cell_count <- pivot_longer(cell_count,-samples) %>% mutate()
cell_count$name <- factor(cell_count$name,levels = c("before_filter","after_filter"))

p1 <- cell_count %>% as.data.frame() %>%
  ggplot(aes(x=samples,y=value,fill=name)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.8) +
  scale_fill_manual(
    limits = c("before_filter","after_filter"),
    values = c("#799BB8", "#EF9EA3")
  )+
  theme_classic() +
  geom_text(aes(label=value,color=name),position=position_dodge(0.9), vjust=0)+
  scale_color_manual(
    limits = c("before_filter","after_filter"),
    values = c("#799BB8", "#EF9EA3")
  )+
  labs(y = "Cell Number")+
  theme(
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(5, "pt"),
    axis.text.x = element_text(color = "black", size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.y = element_text(color = "black", size = 18),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    #legend.key.size = unit(0.7, "cm"),
    #legend.text = element_text(size = 12),
    axis.line = element_line(size = 0.8, arrow = arrow(length = unit(0.2, "cm")))
  ) 
#guides(color = guide_legend(override.aes = list(size = 2)))
p1
pdf(file = "plot_out/sf1.b.barplot_of_cell_number_bef_and_af_filter.pdf", width = 7.3, height = 4.3)
p1
dev.off()

## sf1.c. umap plot of samples and clustering ---------------
load("process_data/sf1.D2.integrated_data_with_clustering_and_no_annotation.Rds")
p1 <- DimPlot(AE_integrate,reduction = "umap",group.by = "orig.ident",pt.size = 0.1)+
  theme_classic() +
  scale_color_manual(
    limits = c("patient-1","patient-2","patient-3","patient-4","patient-5","patient-6"),
    values = samp_col
  ) +
  theme(
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(5, "pt"),
    axis.text = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 18),
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.text = element_text(size = 12),
    axis.line = element_line(size = 0.8, arrow = arrow(length = unit(0.2, "cm")))
  ) +
  labs(x = "UMAP_1", y = "UMAP_2",title = NULL) +
  guides(color = guide_legend(override.aes = list(size = 2)))
p1
pdf(file = "plot_out/sf1.c.UMAP_plot_of_samples.pdf", width = 5, height = 3.6)
p1
dev.off()

umap_col <- c("#9E9DCB","#A3CFBE","#DBB0D4","#57B3C3","#DEA368","#CADD93","#A79798",
              "#929F74","#E68552","#EFE3C3","#679AC3","#BDD9B6","#7B537D","#C1747B")
p2 <- DimPlot(AE_integrate,reduction = "umap",pt.size = 0.2,label = T,label.size = 5)+
  theme_classic() +
  scale_color_manual(
    values = umap_col) +
  theme(
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(5, "pt"),
    axis.text = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 18),
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.text = element_text(size = 12),
    axis.line = element_line(size = 0.8, arrow = arrow(length = unit(0.2, "cm")))
  ) +
  labs(x = "UMAP_1", y = "UMAP_2",title = NULL) +
  guides(color = guide_legend(override.aes = list(size = 2)))
p2
pdf(file = "plot_out/sf1.d.UMAP_plot_of_clusters.pdf", width = 5, height = 3.8)
p2
dev.off()

## sf1.d. bubble plot of marker genes ----------------
load("process_data/sf1.D3.integrated_data_with_clustering_and_annotation.Rds")
p1 <- DotPlot(AE_integrate, features=celltype,cols = c("#E54924", "#498EA4"), cluster.idents = T)+
  theme(
    axis.text.x  = element_text(color="black",size=11,angle = 45,vjust = 1, hjust=1),
    panel.border = element_rect(color="black",size = 1), 
    panel.spacing = unit(1, "mm"),
    axis.line = element_blank(),
    strip.text.x = element_text(size=10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm") 
  )+
  labs(x="", y="")
pdf(file = "plot_out/sf1.e.Dotplot_with_marker_genes.pdf", width = 7, height = 3.7)
p1
dev.off()

# mf --------------------
## mf1.a. umap plot of subcluster EC&SMC -----------
load("process_data/mf1.D1.selected_data_of_EC_and_SMC_no_annotation.Rda")
umap_col <- c("#F6C8A8","#DAA87C","#E89DA0","#88CEE6","#E73D4D","#B696B6","#80C1C4",
              "#B2D3A4","#28996b","#E6CECf","#6d7876","#5D9FE5","#92A5D1")
celltype <- umap.res %>%
  group_by(sub_cluster) %>%
  summarise(UMAP_1=median(UMAP_1),
            UMAP_2=median(UMAP_2))
p1 <- ggplot(umap.res,aes(x=UMAP_1,y=UMAP_2,color=sub_cluster))+
  geom_point(size = 0.4) +
  scale_color_manual(
    limits=c(paste0("SMC",1:11),paste0("EC",1:2)),
    values=umap_col)+
  theme(panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"),
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=10), 
        legend.key.size=unit(0.5,'cm'))+
  theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.10, "inches"), type = "closed")
  )+
  geom_label_repel(data = celltype,aes(label=sub_cluster),size=3,
                   fontface="bold",point.padding=unit(0.1, "lines"))+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
p1
pdf(file = "plot_out/mf1.a.UMAP_plot_of_subcluster_2.pdf", width = 4.6, height = 4.5)
p1
dev.off()

## mf1.b. density plot of SMC&EC ------------
p1 <- ggplot(umap.res, aes(x=UMAP_1,y=UMAP_2)) + 
  geom_point(color="lightgrey", shape=16,size=1) +
  stat_density2d(data=subset(umap.res, cluster == "SMC"), aes(x=UMAP_1, y=UMAP_2, alpha=..level..), 
                 fill="#C25759", size=2, bins=5, geom='polygon',adjust = 0.9) +
  theme(panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"),
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=10), 
        legend.key.size=unit(0.5,'cm'))+
  theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.10, "inches"), type = "closed")
  )+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")+
  scale_alpha_continuous(range=c(0.15,1))
p1

p2 <- ggplot(umap.res, aes(x=UMAP_1,y=UMAP_2)) + 
  geom_point(color="lightgrey", shape=16,size=1) +
  stat_density2d(data=subset(umap.res, cluster == "EC"), aes(x=UMAP_1, y=UMAP_2, alpha=..level..), 
                 fill="#599CB4", size=2, bins=5, geom='polygon',adjust = 0.9) +
  theme(panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"),
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=10), 
        legend.key.size=unit(0.5,'cm'))+
  theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.10, "inches"), type = "closed")
  )+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")+
  scale_alpha_continuous(range=c(0.15,1))
p2
pdf(file = "plot_out/mf1.b.UMAP_plot_of_density_of_SMC&EC.pdf", width = 4.6, height = 8.8)
p1/p2
dev.off()

## mf1.c. feature plot of ACTA2, VWF ----------
DefaultAssay(sub_cluster_1) <- "RNA"
p3 <- FeaturePlot(sub_cluster_1,features = "ACTA2",reduction = "umap",
                  cols = colorRampPalette(c("#f8eced","#d07c7e","#C25759" ))(15)) +
  theme(panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"),
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=10), 
        legend.key.size=unit(0.5,'cm'))+
  theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.10, "inches"), type = "closed")
  )+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold"))
p3  
p4 <- FeaturePlot(sub_cluster_1,features = "VWF",reduction = "umap",
                  cols = colorRampPalette(c("#edf4f7","#7eb2c5","#599CB4" ))(15)) +
  theme(panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"),
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=10), 
        legend.key.size=unit(0.5,'cm'))+
  theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.10, "inches"), type = "closed")
  )+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold"))
p4
pdf(file = "plot_out/mf1.b.UMAP_plot_and_Feature_plot_of_density_of_SMC&EC.pdf", width = 9.2, height = 8)
(p1+p3)/(p2+p4)
dev.off()

##mf1.d. pseudotime -----------
load("process_data/mf1.D2.pseudotime_data_of_SMC&EC.Rda")
p1 = plot_cell_trajectory(ordercds,color_by = "Pseudotime",cell_size=0.5, show_backbone = TRUE)+
  scale_color_continuous(trans = scales::reverse_trans())
p2 = plot_cell_trajectory(ordercds,color_by = "cluster",cell_size=0.5,show_backbone = TRUE) + 
  scale_color_manual(
    limits=c("EC","SMC"),
    values = c("#599CB4","#C25759"))
p3 = plot_cell_trajectory(ordercds,color_by = "sub_cluster",cell_size=0.5, show_backbone = TRUE)+ 
  labs(color="subcluster")+
  scale_color_manual(
    limits=c(paste0("SMC",1:11),paste0("EC",1:2)),
    values = c("#F6C8A8","#DAA87C","#E89DA0","#88CEE6","#E73D4D","#B696B6","#80C1C4",
              "#B2D3A4","#28996b","#E6CECf","#6d7876","#5D9FE5","#92A5D1"))
p4 = plot_cell_trajectory(ordercds,color_by = "State",cell_size=0.5, show_backbone = TRUE) + scale_color_npg()
p5 = plot_complex_cell_trajectory(ordercds,x=1,y=2,color_by = "integrated_snn_res.0.8")+
  scale_color_d3("category20") +
  theme(legend.title = element_blank())
pdf("plot_out/mf1.c.pseudotime_mapping_with_DEGs.pdf",width = 7.2,height = 8)
(p1+p2)/(p4+p3)
dev.off()

## mf1.e. module score of midgenes -------------
load("process_data/mf1.D3.pseudotime_data_of_module_score_of_midgenes.Rda")
p1 <- plot_cell_trajectory(ordercds,color_by = "midstage_genes1",cell_size=0.5, show_backbone = TRUE)+
  scale_color_viridis_c(option="viridis")
pdf(file = "plot_out/mf1.d.pseudotime_trajectory_of_midstage_score.pdf",width = 3.6,height = 4)
p1
dev.off()

all_gene_meta <- mod_subcluster@meta.data %>%
  mutate(final_cluster = case_when(str_detect(sub_cluster,"SMC")~ "SMC",
                                   str_detect(sub_cluster,"EC")~ "EC"))
all_gene_meta$sub_cluster <- as.character(all_gene_meta$sub_cluster) %>%
  factor(.,levels=c("EC1","EC2","SMC7","SMC6","SMC5","SMC4","SMC8","SMC3","SMC2","SMC9","SMC1","SMC11","SMC10"))
p1 <- all_gene_meta %>%
  ggplot(aes(x=sub_cluster, y=midstage_genes1)) +
  geom_violin(aes(fill=sub_cluster),trim = F) +
  geom_boxplot(width=0.2,outlier.shape = NA) +
  facet_grid2(.~final_cluster,
              strip = strip_themed(background_x = elem_list_rect(fill=c("#599CB4","#C25759"))),
              scales="free_x",space="free")+
  scale_fill_manual(limits=c(paste0("SMC",1:11),paste0("EC",1:2)),
                    values = c("#F6C8A8","#DAA87C","#E89DA0","#88CEE6","#E73D4D","#B696B6","#80C1C4",
                               "#B2D3A4","#28996b","#E6CECf","#6d7876","#5D9FE5","#92A5D1"))+
  labs(x=NULL,y="midstage genes score",fill="subcluster")+
  theme(axis.text.x=element_text(angle = 45,vjust=1,hjust=1,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        legend.key.size = unit(0.4, "cm"))
pdf(file = "plot_out/mf1.f.violin_plot_of_module_score_of_midstage_gene_with_all_genes.pdf",width = 8,height = 3)
p1
dev.off()

