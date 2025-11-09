rm(list=ls())
gc()
library(Seurat);library(tidyverse);library(Nebulosa);
library(clustree);library(GEOquery);library(scCustomize);library(scRNAtoolVis)
library(UCell);library(clusterProfiler);library(harmony)
##Import------
data<-getGEO('GSE131907',file ='/Users/huzixin/Brain/download/GSE131907/GSE131907_series_matrix.txt.gz',getGPL = F )
data<-pData(data)
colnames(data)
data<-data[,c(1,2,46:48)]
colnames(data)[3:5]<-c("Patient id","tissue","Stage")
data$tissue<-paste(data$`Patient id`,data$tissue,sep='.')

feature_summary <- read_excel("~/Brain/download/GSE131907/GSE131907_Lung_Cancer_Feature_Summary.xlsx", skip = 1)
as.character(read_excel("~/Brain/download/GSE131907/GSE131907_Lung_Cancer_Feature_Summary.xlsx", skip = 1)[59,1])
#* nLung, normal lung; tLung, early stage tumor lung; 
# tL/B, advanced stage tumor lung;
# mLN, lymph node metastases; 
# nLN, normal lymph node; 
# PE, pleural effusion; 
# mBrain, brain metastases; 
# ADC, Adenocarcinoma; 
# Cur, Current smoker; Ex, Ex-smoker; Never, Never smoker; 
# PD, poorly differentiated; MD, moderately differentiated;, WD, well differentiated
feature_summary <-feature_summary[-59,]
feature_summary <-feature_summary[,-1]
feature_summary$tissue<-paste(feature_summary$`Patient id`,feature_summary$`Tissue origins`,sep='.')
colnames(data)<-str_replace_all(colnames(data),' ','_')
table(data$Tissue_origins)
colnames(feature_summary)
data<-data %>% inner_join(feature_summary[,-1])
dir.create('./01seuRegular/')
data.table::fwrite(data,
                   './01seuRegular/clinic.txt',sep='\t',
                   quote = F,row.names = F)
##Brain-----

data2<-data %>% filter(Tissue_origins%in%c('mBrain'))#,'tL/B'
table(data2$Pathology)
data2$EGFR_status<-ifelse(data2$EGFR%in%c('WT','na'),data2$EGFR,'Mut')

table(data2$EGFR_status)
count<-data.table::fread('/Users/huzixin/Brain/download/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz') %>% as.data.frame() %>% column_to_rownames('Index')
#count_rds<-readRDS('/Users/huzixin/Brain/download/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.rds')
library(readxl)
count[3]%>% tail()

cell_anno<-data.table::fread('~/Brain/download/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz')
cell_anno<-cell_anno %>% as.data.frame() %>% filter(Sample%in%data2$title)
rownames(cell_anno)<-cell_anno$Index
colnames(cell_anno)[3]<-'title'
cell_anno<-cell_anno %>% left_join(data2)
rownames(cell_anno)<-cell_anno$Index
count2<-count[,cell_anno$Index]
table(cell_anno$Cell_type)
count[1:5,1:5]

seu <- CreateSeuratObject(counts = count2,
                          meta.data = cell_anno,
                          project = 'GSE131907_brain')

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu),reduction.name = "pca")

table(seu$title)
seu  <- RunHarmony(seu ,
                   reduction = "pca",
                   group.by.vars = "geo_accession",
                   reduction.save = "harmony")
ElbowPlot(seu, reduction = "pca", ndims = 30)
xx <- cumsum(seu[["harmony"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.90) 
ndim = 30

seu  <- RunUMAP(seu , 
                reduction = "harmony", 
                dims = 1:30,
                reduction.name = "umap")
seu <- RunTSNE(seu, reduction = "harmony", 
               dims = 1:30,
               reduction.name = "tsne", 
               perplexity = 150)


DimPlot(seu, reduction = "umap", group.by = "title",split.by = "EGFR_status")
DimPlot(seu, reduction = "umap", group.by = "Cell_type",label=T)
DimPlot(seu, reduction = "umap", group.by = "Cell_subtype",label=T)
DimPlot(seu, reduction = "umap", group.by = "Cell_type.refined",label=T)

resolution = seq(0.2,2,0.1)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20, k.param = 20)
seu <- FindClusters(seu, resolution = seq(0.2,2,0.1))
pdf('./01seuRegular/seu_brain_umap_resolution.pdf',width=18,height=18)
DimPlot(seu, reduction = "umap", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu, reduction = "tsne", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
dev.off()
seu$seurat_clusters<-as.character(seu$RNA_snn_res.0.2)
table(seu$seurat_clusters)
seu$seurat_clusters[seu$seurat_clusters=='4']<-'3'
seu$seurat_clusters[as.character(seu$RNA_snn_res.1.7)=='12']<-'4'
seu$seurat_clusters[as.character(seu$RNA_snn_res.1.7)=='29']<-'13'

Idents(seu)<-seu$seurat_clusters

(DimPlot(seu,group.by = 'seurat_clusters',label=T)+NoLegend())+
  (DimPlot(seu, reduction = "umap", group.by = "Cell_type",label=T)+NoLegend())&ggsci::scale_color_d3('category20')


marker_genes <- c("KRT8","EPCAM")
marker_genes <- c("ACTA2", "CSPG4", "RGS5")
marker_genes <- c("PECAM1", "VWF", "CLDN5")

marker_genes <- c("GAL3ST1","GJC2","MBP")
marker_genes <- c("GFAP","S100B","ALDH1L1")


marker_genes <- c("MRC1", "CD163", "TREM2")
marker_genes <- c("AIF1","TMEM119","CSF1R")
marker_genes <- c("PLAC8",'CAMP','S100A9')
marker_genes <- c("CX3CR1","CCR7","PTPRC")
marker_genes <- c("IGKC","ITGAL","ITGA4")
marker_genes <- c("CST3","CLEC9A","ITGAM")
marker_genes <- c("FCGR1A","CLEC10A","CD1C")# DC 

# c( 'CAMP','S100A9') #neutrophils


marker_genes <- c("CD3D", "CD28", "CD2")
marker_genes <- c( "NKG7","GNLY")
marker_genes <- c("IGKC","MS4A1", "CD79A")# NK &B
marker_genes<-c( "CCDC53", "RARRES2","PPBP")

seu$celltype<-seu$Cell_type
seu$celltype[seu$seurat_clusters=='13']<-'Pericytes'
seu$celltype[seu$seurat_clusters=='12']<-'Epivascular'
table(seu$celltype)

pdf('./01seuRegular/TumorEpi.pdf',height=10,width=10) 
pdf('./01seuRegular/oligoastro.pdf',height=10,width=10) 
pdf('./01seuRegular/myloid.pdf',height=10,width=10) 
pdf('./01seuRegular/Lym.pdf',height=10,width=10) 
pdf('./01seuRegular/Ependymal.pdf',height=10,width=10) 


p1<-(DimPlot(seu,group.by = 'seurat_clusters',label=T)+NoLegend())+
  (DimPlot(seu, reduction = "umap", group.by = "Cell_type",label=T)+NoLegend())&ggsci::scale_color_d3('category20')

p2<-FeaturePlot(seu,  reduction = "umap",features = marker_genes, order = TRUE,ncol=3)
p3<-VlnPlot(seu, features = marker_genes, slot = "counts", log = TRUE,ncol=3)
p<-(p1/p3/p2)
p
dev.off()
Idents(seu)<-seu$celltype
table(seu$celltype)
Idents(seu)<-factor(seu$celltype,levels = c('Epithelial cells','Myeloid cells',
                                            'Oligodendrocytes',
                                            'Epivascular','Pericytes','T lymphocytes','NK cells','B lymphocytes','MAST cells','Fibroblasts'))
AverageHeatmap(object = seu,
               markerGene = c('KRT8', 'EPCAM',
                              'AIF1', 'MRC1',
                              #'CD1C','CLEC10A',
                              'GFAP', 
                              'GAL3ST1', 
                              'PECAM1','CLDN5',
                              'ACTA2', 'RGS5',
                              'CD3D','CD2',
                              'NKG7','GNLY',
                              'CD79A','MS4A1',
                              'TPSAB1','TPSB2',
                              'COL1A1','COL1A2'),
               clusterAnnoName = F,
               #myanCol = color.set[levels(Idents(seu))],
               annoCol = F,width = 8,height=8)
qs::qsave(seu,'./01seuRegular/seu.anno.qs')

##Dimplot Figure1E,F----
seu<-qs::qread('./01seuRegular/seu.anno.qs')
d <- c(
  "Epithelial cells"="Tumor",#"#6A3D9A"
  "T lymphocytes" = "T cells", #"#C99B7D"
  "NK cells" = "NK cells",# "#FF7F00"
  "B lymphocytes" = "B cells",#"#33A02C"
  "MAST cells"="Mast cells", #"#e1a8b8"
  "Pericytes" = "Pericytes", #"#1F78B4"
  "Epivascular"="Vascular endothelial cells", #"#68AB9F"
  "Myeloid cells"="TAMs & DC", #"#A6CEE3"
  "Oligodendrocytes"="Oligodendrocytes", #"#B2DF8A"
  "Fibroblasts"="Fibroblasts"
)
seu$celltype_plot<-d[seu$celltype]
table(seu$celltype_plot %>% is.na)
color.lung<-c("#6A3D9A","#C99B7D","#FF7F00","#33A02C",
              "#e1a8b8","#1F78B4","#68AB9F","#A6CEE3","#B2DF8A","#FF6F61")
table(seu$celltype)
names(color.lung)<-d
Idents(seu)<-seu$celltype_plot
seu$celltype_plot<-factor(seu$celltype_plot,levels=c('Tumor','TAMs & DC',
                                                     'Oligodendrocytes',
                                                     'Vascular endothelial cells','Pericytes',
                                                     'T cells','NK cells','B cells','Mast cells',
                                                     'Fibroblasts'))

pdf('./01seuRegular/seu.all.dimplot.brain.pdf',height=8,width=8)

clusterCornerAxes(object = seu,reduction = 'umap',
                  clusterCol = 'celltype_plot',
                  noSplit = T,cellLabel = F,
                  keySize = 8,aspect.ratio = 1,
                  relLength = 0.5)& 
  scale_color_manual(values = color.lung) 
dev.off()
Idents(seu)<-factor(seu$celltype_plot,
                    levels=c('Tumor','TAMs & DC',
                             'Oligodendrocytes',
                             'Vascular endothelial cells','Pericytes',
                             'T cells','NK cells','B cells','Mast cells',
                             'Fibroblasts'))
table(seu$celltype_plot)
pdf('./01seuRegular/averagePlot.brain.pdf',height=8,width=8)
AverageHeatmap(object = seu,
               markerGene = c('KRT8', 'EPCAM',
                              'AIF1', 'CLEC10A',
                              'GFAP',
                              'GAL3ST1',# 'GJC2',
                              'PECAM1', 'CLDN5',
                              'ACTA2', 'RGS5',
                              'CD3D', 'CD2',
                              'GNLY','NKG7',
                              'IGKC','MZB1', 
                              'TPSAB1','TPSB2',
                              'COL1A1','COL1A2'),
               clusterAnnoName = F,
               myanCol = color.lung[levels(Idents(seu))],
               annoCol = TRUE,width = 8,height=8)
dev.off()
data.table::fwrite(table(seu$geo_accession,seu$Histology) %>% as.data.frame(),
                   './01seuRegular/Brain.title.txt',sep='\t',
                   quote = F,row.names = F)
##add primary----
dir.create('./05addPrimary')
rm(list=ls());gc()
data<-data.table::fread('./01seuRegular/clinic.txt')
table(data$Tissue_origins)
data2<-data %>% filter(Tissue_origins%in%c('mBrain','tL/B','tLung'))#,'tL/B'
count<-data.table::fread('/Users/huzixin/Brain/download/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz') %>% as.data.frame() %>% column_to_rownames('Index')
cell_anno<-data.table::fread('~/Brain/download/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz')
cell_anno<-cell_anno %>% as.data.frame() %>% filter(Sample%in%data2$title)
rownames(cell_anno)<-cell_anno$Index
colnames(cell_anno)[3]<-'title'
cell_anno<-cell_anno %>% left_join(data2)
rownames(cell_anno)<-cell_anno$Index
count2<-count[,cell_anno$Index]
table(cell_anno$Cell_type)
count[1:5,1:5]

seu <- CreateSeuratObject(counts = count2,
                          meta.data = cell_anno,
                          project = 'GSE131907_brainaddprimary')

rm(count2);rm(count);gc()
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu),reduction.name = "pca")

table(seu$title)
seu  <- RunHarmony(seu ,
                   reduction = "pca",
                   group.by.vars = "geo_accession",
                   reduction.save = "harmony")
ElbowPlot(seu, reduction = "pca", ndims = 30)
xx <- cumsum(seu[["harmony"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.90) 
ndim = 30

seu  <- RunUMAP(seu , 
                reduction = "harmony", 
                dims = 1:30,
                reduction.name = "umap")
seu <- RunTSNE(seu, reduction = "harmony", 
               dims = 1:30,
               reduction.name = "tsne", 
               perplexity = 150)
DimPlot(seu,group.by = 'Cell_type',label=T,split.by = 'Tissue_origins')
seu$Tissue_origins %>% table()
resolution = seq(0.2,2,0.1)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20, k.param = 20)
seu <- FindClusters(seu, resolution = seq(0.2,2,0.1))
pdf('./05addPrimary/seu_umap_resolution.pdf',width=18,height=18)
DimPlot(seu, reduction = "umap", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu, reduction = "tsne", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
dev.off()
DimPlot(seu,group.by='Cell_type')
seu$seurat_clusters<-as.character(seu$RNA_snn_res.0.2)
table(seu$seurat_clusters)
seu$seurat_clusters[as.character(seu$RNA_snn_res.1.4)=='27']<-'22'
Idents(seu)<-seu$seurat_clusters

marker_genes <- c("KRT8","EPCAM")
marker_genes <- c("ACTA2", "CSPG4", "RGS5")
marker_genes <- c("PECAM1", "VWF", "CLDN5")

marker_genes <- c("GAL3ST1","GJC2","MBP")
marker_genes <- c("GFAP","S100B","ALDH1L1")


marker_genes <- c("MRC1", "CD163", "TREM2")
marker_genes <- c("AIF1","TMEM119","CSF1R")
marker_genes <- c("PLAC8",'CAMP','S100A9')
marker_genes <- c("CX3CR1","CCR7","PTPRC")
marker_genes <- c("IGKC","ITGAL","ITGA4")
marker_genes <- c("CST3","CLEC9A","ITGAM")
marker_genes <- c("FCGR1A","CLEC10A","CD1C")# DC 

# c( 'CAMP','S100A9') #neutrophils


marker_genes <- c("CD3D", "CD28", "CD2")
marker_genes <- c( "NKG7","GNLY")
marker_genes <- c("IGKC","MS4A1", "CD79A")# NK &B
marker_genes<-c( "CCDC53", "RARRES2","PPBP")

seu$celltype<-seu$Cell_type
seu$celltype[seu$seurat_clusters=='13']<-'Pericytes'
seu$celltype[seu$seurat_clusters=='12']<-'Epivascular'
table(seu$celltype)

pdf('./05addPrimary/TumorEpi.pdf',height=10,width=10) 
pdf('./05addPrimary/oligoastro.pdf',height=10,width=10) 
pdf('./05addPrimary/myloid.pdf',height=10,width=10) 
pdf('./05addPrimary/Lym.pdf',height=10,width=10) 
pdf('./05addPrimary/Ependymal.pdf',height=10,width=10) 


p1<-(DimPlot(seu,group.by = 'seurat_clusters',label=T)+NoLegend())+
  (DimPlot(seu, reduction = "umap", group.by = "Cell_type",label=T)+NoLegend())

p2<-FeaturePlot(seu,  reduction = "umap",features = marker_genes, order = TRUE,ncol=3)
p3<-VlnPlot(seu, features = marker_genes, slot = "counts", log = TRUE,ncol=3)
p<-(p1/p3/p2)
p
dev.off()

seu$celltype<-seu$Cell_type
seu$celltype[seu$seurat_clusters=='22']<-'Pericytes'
seu$celltype[seu$seurat_clusters=='10']<-'EpiVascular'
head(markers_25$gene,30)
seu$celltype[seu$seurat_clusters=='23']<-'Red_cells'
seu$celltype[seu$seurat_clusters%in%c('18','19','0','21','17')]<-'T_cells'
seu$celltype[seu$seurat_clusters%in%c('25')]<-'NK_cells'
seu$celltype[seu$seurat_clusters%in%c('26','14')]<-'Tumor'
seu$celltype[seu$seurat_clusters%in%c('4','9')]<-'B_cells'
seu$celltype[seu$seurat_clusters%in%c('7')]<-'Mast_cells'
seu$celltype[seu$seurat_clusters%in%c('8')]<-'Fibroblasts'
seu$celltype[seu$seurat_clusters%in%c('11')]<-'Oligo_Astro'
qs::qsave(seu,'./05addPrimary/seu.anno.qs')
##compare markers Brm vs primary------
seu<-qs::qread('./05addPrimary/seu.anno.qs')
table(seu$celltype)
table(seu$Tissue_origins)
seu$resource<-ifelse(seu$Tissue_origins=='mBrain','Brmt','Primary')
table(seu$resource)
seu$celltype_resource<-paste(seu$resource,seu$celltype,sep='.')
Idents(seu) <-seu$celltype_resource
table(seu$celltype_resource)
celltype <- setdiff(unique(seu$celltype),c('Red_cells','Oligo_Astro'))
group <- unique(seu$resource)

DefaultAssay(seu)<-'RNA'

de.list <- lapply(celltype, function(ct) {
  message(glue::glue("processing {ct} ..."))
  ct1 <- paste(group[1],ct,sep='.')
  ct2 <- paste(group[2],ct,sep='.')
  
  de <- FindMarkers(seu, ident.1 = ct2, 
                    ident.2 = ct1, 
                    test.use = "wilcox",
                    logfc.threshold = 0)
  de$change <- ifelse(de$avg_log2FC > 0,
                      paste("up in", group[2]),
                      paste("up in", group[1]))
  de$group <- ct
  return(de)
})
names(de.list) <- celltype
filtered_list <- de.list[!sapply(de.list, is.null)]
filtered_list<-lapply(seq_along(filtered_list),function(xx){
  filtered_list[[xx]]<-filtered_list[[xx]] %>% 
    filter(p_val<0.05) %>% 
    rownames_to_column('gene')
  return(filtered_list[[xx]])
})
de.gene <- purrr::reduce(filtered_list,rbind) %>% 
  arrange(group,desc(avg_log2FC)) %>% 
  filter(abs(avg_log2FC)>=1)
View(de.gene[de.gene$group=='TAM_DC',])
data.table::fwrite(de.gene,
                   './05addPrimary/Diff_brainvsprimary_allcelltype.txt',sep='\t',
                   quote = F,row.names = F)