library(GEOquery)
library(tidyverse)
library(arrayQualityMetrics)
library(affy)
library(limma)
library(Seurat)
library(harmony)
library(Nebulosa)
##Import and Regular-------
file<-list.files('../download/GSE186344/GSE186344_RAW')
file<-file[str_detect(file,'barcodes|mtx|features')]
sample <- unique(str_extract(file,'GSM\\d+'))
pattern <- ".*?(barcodes\\.tsv\\.gz|features\\.tsv\\.gz|matrix\\.mtx\\.gz)"
BM<-unique(str_extract(file[str_detect(file,'Lung|Breast')],'GSM\\d+'))
dir.create('../download/GSE186344/BM')
for (i in 1:18) {
  if(sample[i]%in%LUBM){
    file0 <- file[str_detect(file,sample[i])]
    new_dir<-paste0('../download/GSE186344/BM/',sample[i])
    dir.create(new_dir)
    file.copy(from=paste0('../download/GSE186344/GSE186344_RAW/',file0),
              to=new_dir)
    file0<-list.files(new_dir,full.names = T)
    file.rename(file0,to=paste0(new_dir,"/",sub(pattern, "\\1", file0)))
  }
}
file0<-paste0('../download/GSE186344/BM/',BM,'/')
seu.list <- pbapply::pblapply(file0, function(sn) {
  counts <- Read10X(sn)
  sn <- gsub("_", "-", sn) 
  colnames(counts) <- paste(str_extract(sn,'GSM\\d+'), colnames(counts), sep = "_")
  seu <- CreateSeuratObject(counts = counts,project = unique(str_extract(sn,'GSM\\d+'))[i],)
  return(seu)
})
seu <- base::Reduce(f = merge, x = seu.list)
data<-getGEO('GSE186344',getGPL = F)
data<-pData(data$`GSE186344-GPL24676_series_matrix.txt.gz`)

metadata<-seu@meta.data %>%
  mutate(cellname=rownames(.),
         geo_accession=orig.ident) %>%
  inner_join(data) %>%
  column_to_rownames('cellname')
seu@meta.data<-metadata
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

VlnPlot(seu, features = c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0) +
  geom_hline(yintercept = c(1000, 5000), linetype = "dashed", color = "blue")

VlnPlot(seu, features = c("percent.mt"), group.by = "orig.ident", pt.size = 0, log = T) +
  geom_hline(yintercept = c(10), linetype = "dashed", color = "blue")
table(seu$percent.mt<=20)
seu <-subset(seu,nFeature_RNA>=200 & nFeature_RNA<=5000 & percent.mt<=20)
quantile(seu$nFeature_RNA)

seu <- scutilsR::MarkDoublets(seu, split.by = "orig.ident")
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu),reduction.name = "pca")

table(seu$type,seu$geo_accession)
seu  <- RunHarmony(seu ,
                   reduction = "pca",
                   group.by.vars = "geo_accession",
                   reduction.save = "harmony")
ElbowPlot(seu, reduction = "harmony", ndims = 30)
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



resolution = seq(0.2,2,0.1)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, k.param = 20)
seu <- FindClusters(seu, resolution = seq(0.2,2,0.1))

colnames(seu@meta.data)
pdf('./dataset2/Lung/01merge/seu_umap_resolution.pdf',width=18,height=18)
DimPlot(seu, reduction = "umap", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu, reduction = "tsne", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
dev.off()


seu$seurat_clusters <- as.character(seu$RNA_snn_res.0.8)


seu$seurat_clusters[as.character(seu$RNA_snn_res.1.1)=='17']<-'16'
seu$seurat_clusters[as.character(seu$RNA_snn_res.1.7)=='28']<-'17'
table(seu$seurat_clusters)
Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
Idents(seu) <-factor(Idents(seu),levels=0:21 )


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
markers$gene[as.character(markers$cluster)=='21'] %>% head(20)

file<-list.files('../download/GSE186344/GSE186344_RAW/')
file<-paste0('../download/GSE186344/GSE186344_RAW/',
             file[str_detect(file,'Cell_Types')&str_extract(file,'GSM\\d+')%in%unique(seu$geo_accession)])


metadata<-pbapply::pblapply(file, function(sn){
  data <- data.table::fread(sn)
})
file<-list.files('../download/GSE186344/GSE186344_RAW/')
file<-file[str_detect(file,'Cell_Types')&str_extract(file,'GSM\\d+')%in%unique(seu$geo_accession)]
sample<-str_extract(file,'GSM\\d+')
anno<-data.frame(barcodes=c(paste(sample[1],as.data.frame(metadata[[1]][-1,])[,1],sep=''),
                            paste(sample[2],as.data.frame(metadata[[2]][-1,])[,1],sep='')),
                 celltype=c(as.data.frame(metadata[[1]][-1,])[,2],
                            as.data.frame(metadata[[2]][-1,])[,2]))
class(as.data.frame(metadata[[1]][-1,])[,1])
colnames(anno)[1:2]<-c('barcodes','celltype')      

anno<-anno %>% 
  distinct(barcodes,.keep_all = T)

anno$barcodes<-paste0(anno$barcodes,'-1')
sapply(strsplit(rownames(seu@meta.data),'-'),'[',2) %>% table()

meta_data<-seu@meta.data %>% 
  rownames_to_column('barcodes') %>% 
  left_join(anno) %>% 
  column_to_rownames('barcodes')
table(is.na(meta_data$celltype))

seu@meta.data<-meta_data

data.table::fwrite(markers,
                   './dataset2/Lung/01merge/seurat_cluster_markers.txt',sep='\t',
                   quote = F,row.names = F)

seu2<-subset(seu,seurat_clusters%in%c('7','21'),invert=T )

d <- c(
  "0" = "T_cells",
  "1"="T_NK_cells", 
  "2"="T_NK_cells", 
  "3" = "Pericytes",
  "4" = "T_cells", 
  "5" = "Tumor",
  "6" = "T_cells", 
  #"7" = "T_cells", 
  "8"="EpiVascular", 
  "9"="T_NK_cells", 
  "10"="Plasma_cells" ,
  "11"="TAM",
  "12"="B_cells",
  "13"="Pericytes",
  "14"="Tumor",
  "15"="Epi",
  "16"="Tumor",
  "17"="Plasma_cells",
  "18"="Pericytes",
  "19"="Mast_cells",
  "20"="Oligo_Astro"
)
seu$celltype_cl <- d[seu$seurat_clusters]
qs::qsave(seu,'./dataset2/Lung/03addTNBC/seu.all.lung.breast.qs an')

##Dimplot Figure1E,F------
rm(list=ls());gc()
seu<-qs::qread('./Lung/03addTNBC/seu.all.lung.breast.qs')
DimPlot(seu,group.by = 'celltype_cl',label=T)
rownames(seu@meta.data)<-seu$sample
table(seu$celltype_cl)
table(seu$celltype_pl)
d <- c(
  
  "Tumor"="Tumor",#"#6A3D9A"
  "T_cells" = "T cells", #"#C99B7D"
  "T_NK_cells" = "T cells & NK cells",# "#FF7F00"
  "T_NK"= "T cells & NK cells",
  "B_cells" = "B cells",#"#33A02C"
  'B'= "B cells",
  "Plasma_cells"="B cells", #"#2e8327"
  "Epi"="Alveolar epithelium",#"#4798b3"
  
  "Mast_cells"="Mast cells", #"#e1a8b8"
  "Pericytes" = "Pericytes", #"#1F78B4"
  "pericytes" = "Pericytes",
  "EpiVascular"="Vascular endothelial cells", #"#68AB9F"
  "TAM"="TAMs", #"#A6CEE3"
  "DC"="DC",
  "Oligo_Astro"="Oligodendrocytes", #"#62A3CB" #Astrocytes &
  "MSC_like"="Fibroblasts"
)

seu$celltype_pl <- d[seu$celltype_cl]
table(seu$celltype_pl %>% is.na)
table(seu$celltype_cl)
Idents(seu)<-seu$celltype_pl
color.lung<-c("#6A3D9A","#C99B7D","#FF7F00","#33A02C",
              "#4798b3", "#e1a8b8","#1F78B4","#68AB9F","#A6CEE3",'#C48244',"#B2DF8A","#FF6F61")
names(color.lung)<-unique(d)
library(scCustomize);library(scRNAtoolVis)
seu$celltype_pl<-factor(seu$celltype_pl,levels=c("Tumor",
                                                 "TAMs",
                                                 "DC",
                                                 "Oligodendrocytes",
                                                 "Vascular endothelial cells" ,
                                                 "Pericytes",
                                                 "T cells",
                                                 "T cells & NK cells",
                                                 "B cells",
                                                 #"Plasma cells",
                                                 "Mast cells",
                                                 "Fibroblasts",
                                                 "Alveolar epithelium"))
pdf('./Lung/03addTNBC/seu.all.dimplot.pdf',height=8,width=8)
clusterCornerAxes(object = seu,reduction = 'umap',
                  clusterCol = 'celltype_pl',
                  noSplit = T,cellLabel = F,
                  keySize = 8,aspect.ratio = 1,
                  relLength = 0.5)& 
  scale_color_manual(values = color.lung) 
dev.off()
Idents(seu)<-factor(Idents(seu),levels=c("Tumor",
                                         "TAMs",
                                         "DC",
                                         "Oligodendrocytes",
                                         "Vascular endothelial cells" ,
                                         "Pericytes",
                                         "T cells",
                                         "T cells & NK cells",
                                         "B cells",
                                         #"Plasma cells",
                                         "Mast cells",
                                         "Alveolar epithelium",
                                         "Fibroblasts"))
pdf('./Lung/03addTNBC/averagePlot.pdf',height=8,width=8)

AverageHeatmap(object = seu,
               markerGene = c('KRT8', 'EPCAM',
                              'AIF1', 'TMEM119',
                              'CLEC10A','CD1C',
                              'GFAP', 
                              'GAL3ST1', 
                              'PECAM1','CLDN5',
                              'ACTA2', 'RGS5',
                              'CD3D','CD2',
                              'GNLY','NKG7',
                              'IGKC','MS4A1', 
                             # 'IGHG1','JCHAIN',
                              'TPSAB1','TPSB2',
                              'SFTPA1', 'SFTPA2',
                             'COL1A1','COL1A2'),
               clusterAnnoName = F,
               myanCol = color.lung[levels(Idents(seu))],
               annoCol = T,width = 8,height=8)

dev.off()

data.table::fwrite(table(seu$orig.ident,seu$group) %>% as.data.frame(),
                   './Lung/03addTNBC/title.txt',sep='\t',
                   quote = F,row.names = F)
