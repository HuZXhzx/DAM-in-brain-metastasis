library(Seurat);library(tidyverse);library(Nebulosa);library(clustree);library(GEOquery)
rm(list=ls())
##Routine -----
file<-list.files('/Users/huzixin/Brain/download/GSE143423/')
counts<-data.table::fread('/Users/huzixin/Brain/download/GSE143423/GSE143423_lbm_scRNAseq_gene_expression_counts.csv.gz') %>% column_to_rownames('V1')
meta.data<-data.table::fread('/Users/huzixin/Brain/download/GSE143423/GSE143423_lbm_scRNAseq_metadata.csv.gz',header = T) %>% column_to_rownames('V1')
samples<-colnames(counts)
geo<-GEOquery::getGEO('GSE143423',getGPL = F,file='../download/GSE143423/GSE143423_series_matrix.txt.gz')
geo<-pData(geo)
colnames(geo)
geo<-geo[,c(1,2,39:41)]
colnames(geo)<-c('title','geo_accession','age','disease','sex')
table(geo$disease)
table(meta.data$samples)
colnames(meta.data)[2]<-'title'
meta.data<-meta.data %>% left_join(geo)
rownames(meta.data)<-meta.data$cell.id
intersect(samples,meta.data$cell.id) %>% length()

seu <- CreateSeuratObject(counts = counts,
                          meta.data = meta.data,
                          project = 'GSM4259357')
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
VlnPlot(seu, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"), 
        ncol = 3)
table(seu$percent.mt<=20)
table(seu$nFeature_RNA>=100&seu$nFeature_RNA<=6000&seu$percent.mt<=20,seu$geo_accession) 
quantile(seu$nFeature_RNA)
seu <- scutilsR::MarkDoublets(seu, split.by = "orig.ident")
seu <-subset(seu,nFeature_RNA>=200 & nFeature_RNA<=5000 & percent.mt<=20)



seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu),reduction.name = "pca")

table(seu$geo_accession)
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
seu  <- RunUMAP(seu , 
                reduction = "pca", 
                dims = 1:30,
                reduction.name = "umap_pca")

DimPlot(seu, reduction = "umap_pca", group.by = "geo_accession")+
  DimPlot(seu, reduction = "umap_pca", split.by = "geo_accession")

DimPlot(seu, reduction = "umap", group.by = "geo_accession")+
  DimPlot(seu, reduction = "umap", split.by = "geo_accession")
DimPlot(seu, reduction = "tsne", group.by = "geo_accession")+
  DimPlot(seu, reduction = "tsne", split.by = "geo_accession")


resolution = seq(0.2,2,0.1)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, k.param = 20)
seu <- FindClusters(seu, resolution = seq(0.2,2,0.1))

colnames(seu@meta.data)
pdf('./dataset6/Lung/01merge/seu_umap_resolution.pdf',width=18,height=18)
DimPlot(seu, reduction = "umap", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu, reduction = "tsne", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
dev.off()


seu$seurat_clusters <- as.character(seu$RNA_snn_res.1.2)


seu$seurat_clusters[as.character(seu$RNA_snn_res.1.1)=='17']<-'16'
seu$seurat_clusters[as.character(seu$RNA_snn_res.1.7)=='28']<-'17'
table(seu$seurat_clusters)
Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
Idents(seu) <-factor(Idents(seu),levels=0:22 )
(DimPlot(seu, reduction = "tsne", group.by = "seurat_clusters", label = T)+NoLegend())+
  (DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = T)+NoLegend()) 



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

marker_genes <- c("FCER1A","CLEC9A","CD1C")

marker_genes <- c("CD3D", "CD28", "CD2")
marker_genes <- c( "NKG7","GNLY")
marker_genes <- c("IGKC","MS4A1", "CD79A")# NK &B
marker_genes<-c( "CCDC53", "RARRES2","PPBP") 
marker_genes<-c('COL1A1','COL1A2')
pdf('./dataset6/Lung//01merge/TumorEpi.pdf',height=10,width=10) 
pdf('./dataset6/Lung/01merge/oligoastro.pdf',height=10,width=10) 
pdf('./dataset6/Lung/01merge/myloid.pdf',height=10,width=10) 
# pdf('./dataset6/Lung/01merge/DC.pdf',height=10,width=10) 

pdf('./dataset6/Lung/01merge/Lym.pdf',height=10,width=10) 
pdf('./dataset6/Lung/01merge/Ependymal.pdf',height=10,width=10) 
pdf('./dataset6/Lung/01merge/fibro.pdf',height=10,width=10) 


p1<-(DimPlot(seu, reduction = "umap", label = T)+NoLegend()) +
  (DimPlot(seu, reduction = "umap", group.by='geo_accession',label = T)+NoLegend())
p2<-FeaturePlot(seu,  reduction = "umap",features = marker_genes, order = TRUE,ncol=3)
p3<-VlnPlot(seu, features = marker_genes, slot = "counts", log = TRUE,ncol=3)
p<-(p1/p3/p2)
p
dev.off()

Idents(seu)
markers<-FindAllMarkers(seu, only.pos = TRUE)
markers <-markers%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(avg_log2FC))%>%
  ungroup()
markers_filter <-markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:10) %>% 
  ungroup()
markers$gene[as.character(markers$cluster)=='19'] %>% head(40)
data.table::fwrite(markers,
                   './dataset6/Lung/01merge/all_seurat_cluster_markers.txt',sep='\t',
                   quote = F,row.names = F)
#0,1,2,3,5,6,8,9,10,13,18,22 Tumor
#7,4,14,21,17 TAM 
#在TAM中聚类后发现14为肺巨噬细胞 17为肿瘤细胞
#11 pericytes 12 epivascular 
#15 ﬁbroblasts
#B 16
#19 oligo_astro?
#20 CD8+T NK 
d <- c(
  "0" = "Tumor",
  "1"="Tumor", 
  "2"="Tumor", 
  "3" = "Tumor",
  "4" = "TAM", 
  "5" = "Tumor",
  "6" = "Tumor", 
  "7" = "TAM", 
  "8"="Tumor", 
  "9"="Tumor", 
  "10"="Tumor" ,
  "11"="Pericytes",
  "12"="EpiVascular",
  "13"="Tumor",
  "14"="TAM",
  "15"="Fibroblasts",
  "16"="B_cells",
  "17"="Tumor",
  "18"="Tumor",
  "19"="Oligo_Astro",
  "20"="T_NK_cells",
  "21"="TAM",
  "22"="Tumor"
)
seu$celltype <- d[seu$seurat_clusters]
table(seu$celltype)
d <- c(
  "Tumor"="Tumor",#"#6A3D9A"
  "T_NK_cells" = "T cells & NK cells",# "#FF7F00"
  "B_cells" = "B cells",#"#33A02C"
  "Pericytes" = "Pericytes", #"#1F78B4"
  "EpiVascular"="Vascular endothelial cells", #"#68AB9F"
  "TAM"="TAMs & DC", #"#A6CEE3"
  "Oligo_Astro"="Oligodendrocytes", #"#62A3CB"
  "Fibroblasts"="Fibroblasts"
)
##Dimplot Figure1E,F-------
seu$celltype_pl <- d[seu$celltype]
table(seu$celltype_pl %>% is.na)
table(seu$celltype)
Idents(seu)<-seu$celltype_pl
color.lung<-c("#6A3D9A","#FF7F00","#33A02C",
              "#1F78B4","#68AB9F","#A6CEE3","#B2DF8A","#FF6F61")
names(color.lung)<-d
library(scCustomize);library(scRNAtoolVis)
seu$celltype_pl<-factor(seu$celltype_pl,levels=c("Tumor",
                                                 "TAMs & DC",
                                                 "Oligodendrocytes",
                                                 "Vascular endothelial cells" ,
                                                 "Pericytes",
                                                 "T cells & NK cells",
                                                 "B cells",
                                                 "Fibroblasts"))
pdf('./Lung/03addTNBC/seu.all.dimplot.pdf',height=8,width=8)
clusterCornerAxes(object = seu,reduction = 'umap',
                  clusterCol = 'celltype_pl',
                  noSplit = T,cellLabel = F,
                  keySize = 8,aspect.ratio = 1,
                  relLength = 0.5)& 
  scale_color_manual(values = color.lung) 
dev.off()
Idents(seu)<-factor(Idents(seu),levels=c("Tumor",
                                         "TAMs & DC",
                                         "Oligodendrocytes",
                                         "Vascular endothelial cells" ,
                                         "Pericytes",
                                         "T cells & NK cells",
                                         "B cells",
                                         "Fibroblasts"))
pdf('./Lung/03addTNBC/averagePlot.pdf',height=8,width=8)

AverageHeatmap(object = seu,
               markerGene = c('KRT8', 'EPCAM',
                              'AIF1', 'CLEC10A',
                              'GFAP', 
                              'GAL3ST1', 
                              'PECAM1','CLDN5',
                              'ACTA2', 'RGS5',
                              'CD3D',
                              'NKG7',
                              'JCHAIN','MZB1',
                              'COL1A1','COL1A2'),
               clusterAnnoName = F,
               myanCol = color.lung[levels(Idents(seu))],
               annoCol = T,width = 8,height=8)
dev.off()
qs::qsave(seu,'./dataset6/Lung/03addTNBC/seu.all.lung.breast.qs')

