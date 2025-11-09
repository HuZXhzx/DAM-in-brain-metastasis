rm(list=ls())
library(GEOquery)
library(tidyverse)
library(arrayQualityMetrics)
library(affy)
library(limma)
library(Seurat)
library(harmony)
library(Nebulosa)
##Import and Routine Primary and brain------
data<-getGEO('GSE123902',file ='/Users/huzixin/Brain/download/GSE123902/GSE123902_series_matrix.txt.gz',getGPL = F )
data<-pData(data)
colnames(data)
data[1,]
data<-data[,c(1,2,8,11:13)]
table(data$source_name_ch1)
colnames(data)<-c('title','geo_accession','resource','Stage','histology','chemotherapy')
data<-data %>% as.data.frame() %>% filter(str_detect(resource,'BRAIN')|resource=='Patient primary lung adenocarcinoma')
table(data$resource)
table(data$Stage)
test<-data.table::fread('../download/GSE123902/GSE123902_RAW/GSM3516663_MSK_LX661_PRIMARY_TUMOUR_dense.csv.gz')
rm(test)
files<-list.files('../download/GSE123902/GSE123902_RAW',full.names = T)
files<-data.frame(files=files,geo_accession=str_extract(files,'GSM\\d+'))
files<-files %>% inner_join(data)

for (i in 1:11) {
  
  file0 <- files$files[i]
  new_dir<-paste0('../download/GSE123902/new/',files$geo_accession[i])
  dir.create(new_dir)
  file.copy(from=paste0('../download/GSE123902/GSE123902_RAW/',file0),
            to=new_dir)
  file0<-list.files(new_dir,full.names = T)
  file.rename(file0,to=paste0(new_dir,"/matrix.mtx.gz",))
  
}

seu.list <- pbapply::pblapply(files$files, function(sn) {
  
  counts <-data.table::fread(sn) %>% column_to_rownames('V1')
  rownames(counts)<-paste0(str_extract(sn,'GSM\\d+'),'-',rownames(counts))
  counts<-t(counts)
  seu <- CreateSeuratObject(counts = counts,
                            project = str_extract(sn,'GSM\\d+'),)
  return(seu)
})
seu <- base::Reduce(f = merge, x = seu.list)
seu@meta.data$geo_accession<-seu$orig.ident
seu$sample<-rownames(seu@meta.data)
seu@meta.data<-seu@meta.data %>% inner_join(data)
rownames(seu@meta.data)<-seu$sample
rownames(seu@meta.data)
table(seu$resource)
View(seu@meta.data)
seu$resource[str_detect(seu$resource,'BRAIN')]<-'brain'
seu$resource[!(str_detect(seu$resource,'brain'))]<-'primary'

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT")

VlnPlot(seu, features = c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0) +
  geom_hline(yintercept = c(1000, 5000), linetype = "dashed", color = "blue")

VlnPlot(seu, features = c("percent.mt"), group.by = "orig.ident", pt.size = 0, log = T) +
  geom_hline(yintercept = c(10), linetype = "dashed", color = "blue")
table(seu$percent.mt<=20)
quantile(seu$nFeature_RNA)

seu <- scutilsR::MarkDoublets(seu, split.by = "orig.ident")
seu <-subset(seu,nFeature_RNA>=200 & nFeature_RNA<=5000 & percent.mt<=20)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu),reduction.name = "pca")

table(seu$resource,seu$geo_accession)
seu  <- RunHarmony(seu ,
                   reduction = "pca",
                   group.by.vars = "geo_accession",
                   reduction.save = "harmony")
ElbowPlot(seu, reduction = "harmony", ndims = 30)
xx <- cumsum(seu[["pca"]]@stdev^2)
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
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20, k.param = 20)
seu <- FindClusters(seu, resolution = seq(0.2,2,0.1))

colnames(seu@meta.data)
pdf('./Lung_to_brain/GSE123902/01merge/seu_umap_resolution.pdf',width=18,height=18)
DimPlot(seu, reduction = "umap", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu, reduction = "tsne", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
dev.off()


seu$seurat_clusters <- as.character(seu$RNA_snn_res.0.5)
seu$seurat_clusters[as.character(seu$RNA_snn_res.1.1)=='17']<-'16'
seu$seurat_clusters[as.character(seu$RNA_snn_res.1.7)=='28']<-'17'
table(seu$seurat_clusters)
Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
Idents(seu) <-factor(Idents(seu),levels=0:17 )
(DimPlot(seu, reduction = "tsne", group.by = "seurat_clusters", label = T)+NoLegend())+
  (DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = T)) & ggsci::scale_color_d3("category20")
table(seu$DF.classifications,seu$seurat_clusters)

seu$seurat_clusters[as.character(seu$RNA_snn_res.1.7)=='12']<-'18'

Idents(seu)<-seu$seurat_clusters
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
data.table::fwrite(markers,
                   './Lung_to_brain/GSE123902/01merge/seurat_cluster_markers.txt',sep='\t',
                   quote = F,row.names = F)
pdf('./Lung_to_brain/GSE123902/01merge/Seurat_clusters_resource.pdf',height=4,width=8) 
(DimPlot(seu, reduction = "umap", split.by = 'resource',group.by = "seurat_clusters", label = T)) & ggsci::scale_color_d3("category20")
dev.off()
d <- c(
  "0" = "T_cells",
  "1"="T_cells", #"#CAB2D6"
  "2"="NK_cells", #"#62A3CB"
  "3" = "TAM",#"#C99B7D"
  "4" = "Tumor", #"#8DC594"
  "5" = "B_cells",#"#33A02C"
  "6" = "TAM", #"#B2DF8A"
  "7" = "Tumor", #"#1F78B4"
  "8"="Pericytes", #"#68AB9F"
  "9"="T_cells", #"#A6CEE3"
  "10"="DC" ,#"#FDBF6F"
  "11"="Mast_cells",
  "12"="Plasma_cells",
  "13"="EpiVascular",
  "14"="Tumor",
  "15"="Tumor",
  "16"="TAM",
  "17"="Brain_Epi"
)
seu$celltype <- d[seu$seurat_clusters]
table(seu$celltype)

##Dimplot ----

d <- c(
  "Tumor"="Tumor", #"#CAB2D6"
  "T_cells" = "T cells & NK cells",#"#C99B7D"
  "NK_cells" = "T cells & NK cells", #"#8DC594"
  "B_cells" = "B cells",#"#33A02C"
  "Plasma_cells"="B cells", #"#2e8327"
  # "Brain_Epi"="Brain endothelial cells", #"#4798b3"
  "Mast_cells"="Mast cells", #"#e1a8b8"
  "Pericytes" = "Pericytes", #"#1F78B4"
  "EpiVascular"="Vascular endothelial cells", #"#68AB9F"
  "TAM"="TAMs & DC", #"#A6CEE3"
  "DC"="TAMs & DC"#"#C48244"
  
)
seu$celltype_pl <- d[seu$celltype]
table(seu$celltype_pl %>% is.na())

color.lung<-c("#6A3D9A","#FF7F00","#33A02C",
              "#e1a8b8","#1F78B4","#68AB9F","#A6CEE3")
names(color.lung)<-unique(d)
Idents(seu)<-seu$celltype_pl
seu$celltype_pl<-factor(seu$celltype_pl,levels=c('Tumor',
                                                 'TAMs & DC',
                                                 'Vascular endothelial cells',
                                                 'Pericytes',
                                                 #'Plasma cells',
                                                 'T cells & NK cells',
                                                 'B cells',
                                                 'Mast cells'))
pdf('./01merge/seu.all.dimplot.brain.pdf',height=8,width=8)
clusterCornerAxes(object = seu,reduction = 'umap',
                  clusterCol = 'celltype_pl',
                  noSplit = T,cellLabel = F,
                  keySize = 8,aspect.ratio = 1,
                  relLength = 0.5)& 
  scale_color_manual(values = color.lung) 
dev.off()

table(seu$celltype_pl %>% is.na())
Idents(seu)<-factor(seu$celltype_pl,
                    levels=c('Tumor',
                             'TAMs & DC',
                             'Vascular endothelial cells',
                             'Pericytes',
                             #'Plasma cells',
                             'T cells & NK cells',
                             'B cells',
                             'Mast cells'))
pdf('./01merge/averagePlot.brain.pdf',height=8,width=8)
AverageHeatmap(object = seu,
               markerGene = c('KRT8', 'EPCAM',
                              'AIF1',#'TMEM119',
                              'CLEC10A',#'CD1C',
                              # 'GFAP','ALDOC'
                              # 'GAL3ST1', 'GJC2',
                              'PECAM1', 'CLDN5',
                              'ACTA2', 'RGS5',
                              
                              #'IGHG1','JCHAIN',
                              'CD3D', 'NKG7',
                              'IGKC','MZB1', 
                              'TPSAB1','TPSB2'),
               clusterAnnoName = F,
               myanCol = color.lung[levels(Idents(seu))],
               annoCol = TRUE,width = 8,height=8)
dev.off()

data.table::fwrite(table(seu$geo_accession,seu$histology) %>% as.data.frame(),
                   './01merge/Brain.title.txt',sep='\t',
                   quote = F,row.names = F)
qs::qsave(seu,'./01merge/seu.brain.all.qs')


##findmarkers brain meta and primary=====
rm(list=ls());gc()
seu<-qs::qread('./01merge/seu_human_annotate.qs')
table(seu$celltype)
seu$celltype_resource<-paste(seu$resource,seu$celltype,sep='.')
Idents(seu) <-seu$celltype_resource

celltype <- unique(seu$celltype)
group <- unique(seu$resource)

seu2 <- subset(seu, resource == group[1])
Idents(seu2) <- seu2$celltype
baseline.levels <- AverageExpression(seu2, assays = "RNA")$RNA
DefaultAssay(seu)<-'RNA'

de.list <- lapply(celltype, function(ct) {
  message(glue::glue("processing {ct} ..."))
  ct1 <- paste(group[1],ct,sep='.')
  ct2 <- paste(group[2],ct,sep='.')
  n1<-length(colnames(subset(seu,
                             celltype_resource==ct1)))
  n2<-length(colnames(subset(seu,
                             celltype_resource==ct2)))
  if(n1>20&n2>20){
    de <- FindMarkers(seu, ident.1 = ct2, 
                      ident.2 = ct1, 
                      test.use = "wilcox",
                      logfc.threshold = 0)
    de$change <- ifelse(de$avg_log2FC > 0,
                        paste("up in", group[2]),
                        paste("up in", group[1]))
    de$diff_rate <- de$avg_log2FC / baseline.levels[rownames(de), ct]
    de$group <- ct
    return(de)
  }
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
View(de.gene[de.gene$group=='Tumor',])
data.table::fwrite(de.gene,
                   './01merge/Diff_brainvsprimary_allcelltype.txt',sep='\t',
                   quote = F,row.names = F)

##Brain Figure1E,F----
seu<-subset(seu,resource=='brain')
table(seu$celltype)
seu<-subset(seu,celltype=='Brain_Epi',invert=T)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu),reduction.name = "pca")

table(seu$resource,seu$geo_accession)
seu  <- RunHarmony(seu ,
                   reduction = "pca",
                   group.by.vars = "geo_accession",
                   reduction.save = "harmony")
ElbowPlot(seu, reduction = "harmony", ndims = 30)
xx <- cumsum(seu[["pca"]]@stdev^2)
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

d <- c(
  "Tumor"="Tumor", #"#CAB2D6"
  "T_cells" = "T cells & NK cells",#"#C99B7D"
  "NK_cells" = "T cells & NK cells", #"#8DC594"
  "B_cells" = "B cells",#"#33A02C"
  "Plasma_cells"="B cells", #"#2e8327"
  # "Brain_Epi"="Brain endothelial cells", #"#4798b3"
  "Mast_cells"="Mast cells", #"#e1a8b8"
  "Pericytes" = "Pericytes", #"#1F78B4"
  "EpiVascular"="Vascular endothelial cells", #"#68AB9F"
  "TAM"="TAMs & DC", #"#A6CEE3"
  "DC"="TAMs & DC"#"#C48244"
  
)
seu$celltype_pl <- d[seu$celltype]
table(seu$celltype_pl %>% is.na())

color.lung<-c("#6A3D9A","#FF7F00","#33A02C",
              "#e1a8b8","#1F78B4","#68AB9F","#A6CEE3")
names(color.lung)<-unique(d)
Idents(seu)<-seu$celltype_pl
seu$celltype_pl<-factor(seu$celltype_pl,levels=c('Tumor',
                                                 'TAMs & DC',
                                                 'Vascular endothelial cells',
                                                 'Pericytes',
                                                 #'Plasma cells',
                                                 'T cells & NK cells',
                                                 'B cells',
                                                 'Mast cells'))
pdf('./01merge/seu.all.dimplot.brain.pdf',height=8,width=8)
clusterCornerAxes(object = seu,reduction = 'umap',
                  clusterCol = 'celltype_pl',
                  noSplit = T,cellLabel = F,
                  keySize = 8,aspect.ratio = 1,
                  relLength = 0.5)& 
  scale_color_manual(values = color.lung) 
dev.off()

table(seu$celltype_pl %>% is.na())
Idents(seu)<-factor(seu$celltype_pl,
                    levels=c('Tumor',
                             'TAMs & DC',
                             'Vascular endothelial cells',
                             'Pericytes',
                             #'Plasma cells',
                             'T cells & NK cells',
                             'B cells',
                             'Mast cells'))
pdf('./01merge/averagePlot.brain.pdf',height=8,width=8)
AverageHeatmap(object = seu,
               markerGene = c('KRT8', 'EPCAM',
                              'AIF1',#'TMEM119',
                              'CLEC10A',#'CD1C',
                              # 'GFAP','ALDOC'
                              # 'GAL3ST1', 'GJC2',
                              'PECAM1', 'CLDN5',
                              'ACTA2', 'RGS5',
                              
                              #'IGHG1','JCHAIN',
                              'CD3D', 'NKG7',
                              'IGKC','MZB1', 
                              'TPSAB1','TPSB2'),
               clusterAnnoName = F,
               myanCol = color.lung[levels(Idents(seu))],
               annoCol = TRUE,width = 8,height=8)
dev.off()

data.table::fwrite(table(seu$geo_accession,seu$histology) %>% as.data.frame(),
                   './01merge/Brain.title.txt',sep='\t',
                   quote = F,row.names = F)
qs::qsave(seu,'./01merge/seu.brain.all.qs')


