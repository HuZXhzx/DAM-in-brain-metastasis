library(GEOquery)
library(tidyverse)
library(arrayQualityMetrics)
library(affy)
library(limma)
library(Seurat)
library(harmony)
library(Nebulosa)
##Import and Regular-----
file<-list.files('../download/GSE234832/GSE234832_RAW/')
file[1]
data<-getGEO('GSE234832',file = '/Users/huzixin/Brain/download/GSE234832/GSE234832_series_matrix.txt',getGPL = F)
data<-pData(data)
data$geo_accession
table(str_extract(file,'GSM\\d+'))
intersect(unique(str_extract(file,'GSM\\d+')),data$geo_accession)
strsplit(file0,'_([A-Z]{6}\\d+.)')
for (i in 1:5) {
  sample <- unique(str_extract(file,'GSM\\d+'))
  file0 <- file[str_detect(file,sample[i])]
  new_dir<-paste0('../download/GSE234832/GSE234832_RAW/',sample[i],'/')
  dir.create(new_dir)
  file.copy(from=paste0('../download/GSE234832/GSE234832_RAW/',file0),
            to=new_dir)
  file0<-list.files(new_dir,full.names = T)
  
  file.rename(file0,
              to=paste0(new_dir,'/',sapply(strsplit(file[str_detect(file,sample[i])],'_([A-Z]{6}\\d+.)'),'[',2)))
}
file0<-paste0('../download/GSE234832/GSE234832_RAW/',unique(str_extract(file,'GSM\\d+')),'/')
sclist<-{}
for(i in 1:5){
  scdata <- Read10X(data.dir =file0[i])
  scobj <- CreateSeuratObject(counts = scdata, 
                              project = unique(str_extract(file,'GSM\\d+'))[i], 
                              min.cells = 3, 
                              min.features = 200)
  metadata = scobj@meta.data
  scobj@meta.data$group = unique(str_extract(file,'GSM\\d+'))[i]
  scobj<-RenameCells(scobj,
                     add.cell.id=unique(str_extract(file,'GSM\\d+'))[i])
  sclist[[i]]<-scobj
  
}

seu <- merge(x=sclist[[1]], 
               y=sclist[2:length(sclist)])
metadata<-seu@meta.data %>%
  mutate(cellname=rownames(.),
         geo_accession=orig.ident) %>%
  inner_join(data[,c('title','geo_accession')]) %>%
  mutate(group=str_extract(title,'([A-Z]{6})')) %>% 
  column_to_rownames('cellname')

seu@meta.data<-metadata

seu <-seu[[2]]
View(seu@meta.data)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <-subset(seu,nFeature_RNA>=200 & nFeature_RNA<=5000 & percent.mt<=20)
seu <- scutilsR::MarkDoublets(seu, split.by = "orig.ident")
rownames(seu@meta.data)
table(seu$type)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu),reduction.name = "pca")
seu$geo_accession %>% table()
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

DimPlot(seu, reduction = "umap", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu, reduction = "tsne", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)


seu$seurat_clusters <- as.character(seu$RNA_snn_res.0.4)


table(seu$seurat_clusters)
Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
Idents(seu) <-factor(Idents(seu),levels=0:9 )

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
markers$gene[as.character(markers$cluster)=='3'] %>% head(20)

d <- c(
  "0" = "Tumor",
  "1"="Tumor", 
  "2"="TAM", 
  "3" = "Oligo_Astro",
  "4" = "T_NK_cells", 
  "5" = "Tumor",
  "6" = "TAM", 
  "7" = "Tumor", 
  "8"="Oligo_Astro", 
  "9"="Pericytes", 
  "10"="EpiVascular",
  "11"="B_cells",
  "12"="Oligo_Astro")
seu$celltype<- d[seu$seurat_clusters]
table(seu$celltype)

qs::qsave(seu,'./dataset2/Lung/03addTNBC/seu.all.lung.breast.qs')

##Dimplot Figure1E,F----
table(seu$celltype)
DimPlot(seu,group.by = 'celltype',label=T)
d <- c(
  
  "Tumor"="Tumor",#"#6A3D9A"
  "T_NK_cells" = "T cells & NK cells", #"#C99B7D"
  "T"="T cells & NK cells",
  "NK"="T cells & NK cells",
  "B_cells" = "B cells",#"#33A02C"
  "B"= "B cells",
  "Pericytes" = "Pericytes",
  "pericytes"="Pericytes",#"#1F78B4"
  "EpiVascular"="Vascular endothelial cells", #"#68AB9F"
  "Epivascular"="Vascular endothelial cells", 
  "TAM"="TAMs & DC", #"#A6CEE3"
  "DC"="TAMs & DC",
  "Oligo_Astro"="Oligodendrocytes" ,
  "Astro"="Oligodendrocytes" ,
  "Oligo"="Oligodendrocytes" 
)
color.set<-c("#6A3D9A","#FF7F00","#33A02C","#1F78B4","#68AB9F","#A6CEE3","#B2DF8A")
seu$celltype_pl <- d[seu$celltype]
table(seu$celltype_pl %>% is.na())
table(seu$celltype_pl)
names(color.set)<-unique(d)
Idents(seu)<-seu$celltype_pl
library(scCustomize);library(scRNAtoolVis)
seu$celltype_pl<-factor(seu$celltype_pl,levels=c("Tumor",
                                                 "TAMs & DC",
                                                 
                                                 "Oligodendrocytes",
                                                 "Vascular endothelial cells" ,
                                                 "Pericytes",
                                                 "T cells & NK cells",
                                                 "B cells"))
pdf('./Lung/03addBreast/seu.all.dimplot.pdf',height=8,width=8)
clusterCornerAxes(object = seu,reduction = 'umap',
                  clusterCol = 'celltype_pl',
                  noSplit = T,cellLabel = F,
                  keySize = 8,aspect.ratio = 1,
                  relLength = 0.5)& 
  scale_color_manual(values = color.set) 
dev.off()
Idents(seu)<-factor(Idents(seu),levels=c("Tumor",
                                         "TAMs & DC",
                                         "DC",
                                         "Oligodendrocytes",
                                         "Vascular endothelial cells" ,
                                         "Pericytes",
                                         "T cells & NK cells",
                                         "B cells"))
pdf('./Lung/03addBreast/averagePlot.pdf',height=8,width=8)
AverageHeatmap(object = seu,
               markerGene = c('KRT8', 'EPCAM',
                              'AIF1',
                              'CLEC10A',
                              'GFAP', 
                              'GAL3ST1', 
                              'PECAM1','CLDN5',
                              'ACTA2', 'RGS5',
                              'CD3D','NKG7',
                              'IGKC','MS4A1'),
               clusterAnnoName = F,
               myanCol = color.set[levels(Idents(seu))],
               annoCol = T,width = 8,height=8)

dev.off()

