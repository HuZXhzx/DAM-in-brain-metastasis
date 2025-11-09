library(GEOquery)
library(tidyverse)
library(arrayQualityMetrics)
library(affy)
library(limma)
library(Seurat)
library(harmony)
##Import------
data<-getGEO('GSE202371',file ='/Users/huzixin/Brain/download/GSE202371/GSE202371_series_matrix.txt',getGPL = F )
data<-pData(data)
colnames(data)
data[1,]
data<-data[,c(1,2,40,41)]
colnames(data)<-c('title','geo_accession','diagnosis','location')
table(data$diagnosis)
data<-data %>% filter(diagnosis=='Lung adenocarcinoma')


files<-list.files('/Users/huzixin/Brain/download/GSE202371/GSE202371_RAW',full.names = T)
files<-files[str_extract(files,'GSM\\d+')%in%data$geo_accession]
files_cellanno<-files[str_detect(files,'cellname.list')]
files_count<-files[str_detect(files,'counts')]

test<-data.table::fread(files[2])
rm(test)
test2<-data.table::fread(files[1])

files_cellanno<-data.frame(files=files_cellanno,
                           geo_accession=str_extract(files_cellanno,'GSM\\d+')) %>% 
  column_to_rownames('geo_accession')
files_cellanno<-files_cellanno[data$geo_accession,]
files_count<-data.frame(files=files_count,
                           geo_accession=str_extract(files_count,'GSM\\d+'))%>% 
  column_to_rownames('geo_accession')
files_count<-files_count[data$geo_accession,]



seu.list <- pbapply::pblapply(1:10, function(sn) {
  
  counts <-data.table::fread(files_count[sn]) %>% column_to_rownames('gene')
  counts<-t(counts)
  cellindex<-rownames(counts)
  metadata<-data.table::fread(files_cellanno[sn]) 
  
  metadata$CellName<-paste0(str_extract(files_cellanno[sn],'GSM\\d+'),'-',metadata$CellName)
  metadata<-metadata %>% as.data.frame()
  rownames(metadata)<-metadata$CellIndex
  metadata<-metadata[cellindex,]
  
  rownames(counts)<-metadata$CellName
  rownames(metadata)<-metadata$CellName
  seu <- CreateSeuratObject(counts = t(counts),
                            meta.data = metadata,
                            project = str_extract(files_count[sn],'GSM\\d+'),)
  return(seu)
})




seu <- base::Reduce(f = merge, x = seu.list)
#cellanno  <- purrr::reduce(cellanno ,rbind) 
seu@meta.data$geo_accession<-seu$orig.ident
seu$sample<-rownames(seu@meta.data)
seu@meta.data<-seu@meta.data %>% inner_join(data)
rownames(seu@meta.data)<-seu$sample
rownames(seu@meta.data)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT")

VlnPlot(seu, features = c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0) +
  geom_hline(yintercept = c(1000, 5000), linetype = "dashed", color = "blue")

VlnPlot(seu, features = c("percent.mt"), group.by = "orig.ident", pt.size = 0, log = T) +
  geom_hline(yintercept = c(100), linetype = "dashed", color = "blue")
table(seu$percent.mt<=80)

quantile(seu$nFeature_RNA)

seu <- scutilsR::MarkDoublets(seu, split.by = "orig.ident")

seu <-subset(seu,nFeature_RNA>=200 & nFeature_RNA<=5000 & percent.mt<=20)

dir.create('./01merge')
qs::qsave(seu,'./01merge/seu_human.qs')
##seu regular-----
rm(list = ls())
seu<-qs::qread('./01merge/seu_human.qs')
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu),reduction.name = "pca")

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


DimPlot(seu, reduction = "umap", group.by = "title")
DimPlot(seu, reduction = "umap", split.by = "title")

resolution = seq(0.2,2,0.1)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20, k.param = 20)
seu <- FindClusters(seu, resolution = seq(0.2,2,0.1))

colnames(seu@meta.data)
pdf('./01merge/seu_umap_resolution.pdf',width=18,height=18)
DimPlot(seu, reduction = "umap", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu, reduction = "tsne", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
dev.off()

DimPlot(seu,group.by='RNA_snn_res.0.8',label=T)
seu$seurat_clusters <- as.character(seu$RNA_snn_res.0.8)

table(seu$seurat_clusters)
Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
Idents(seu) <-factor(Idents(seu),levels=0:23)
(DimPlot(seu, reduction = "tsne", group.by = "seurat_clusters", label = T)+NoLegend())+
  (DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = T)) 
table(seu$DF.classifications,seu$seurat_clusters)
##Findmarkers annotate
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
marker_genes <- c('TPSAB1','TPSB2')
marker_genes <- c('SPARC','GNG11')#brain epi
marker_genes <- c('SFTPA1', 'SFTPA2') #lung epi
marker_genes<-c( "CCDC53", "RARRES2","PPBP")

pdf('./01merge/TumorEpi.pdf',height=10,width=10) 
pdf('./01merge/oligoastro.pdf',height=10,width=10) 
pdf('./01merge/myloid.pdf',height=10,width=10) 
pdf('./01merge/Lym.pdf',height=10,width=10) 


p1<-(DimPlot(seu, reduction = "umap", label = T)+NoLegend()) +
  (DimPlot(seu, reduction = "tsne", label = T)+NoLegend()) 

p2<-FeaturePlot(seu,  reduction = "umap",features = marker_genes, order = TRUE,ncol=3)
p3<-VlnPlot(seu, features = marker_genes, slot = "counts", log = TRUE,ncol=3)
p<-(p1/p3/p2)
p
dev.off()
DimPlot(seu,group.by='RNA_snn_res.1.3',label=T)
seu$seurat_clusters[as.character(seu$RNA_snn_res.1.3)==21&seu$seurat_clusters=='11'] <-'24'
d<-c('0'='T_NK_cells',
    '1'='Tumor',
     '2'='Oligo_Astro',#
     '3'='TAM',
    '4'='T_NK_cells',
     '5'='TAM',
     '6'='TAM',
     '7'='Tumor',
    '8'='T_cells',#
     '9'='TAM',
     '10'='TAM',
     '11'='Tumor', #recycling
     '12'='Tumor',
     '13'='Oligo_Astro',#
     '14'='Pericytes',
     '15'='DC',#
      '16'='T_NK_cells',
     '17'='B',#
     '18'='Tumor',
     '19'='B',
     '20'='EpiVascular',
     '21'='Pericytes',
     '22'='Mast_cells',
     '23'='Unknown', #去掉
     '24'='Tumor')#?
Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
Idents(seu) <-factor(Idents(seu),levels=0:24)
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
View(markers[as.character(markers$cluster)=='11',] %>% head(100))
markers$gene[as.character(markers$cluster)=='23'] %>% head(20)
data.table::fwrite(markers,
                   './01merge/seurat_cluster_markers.txt',sep='\t',
                   quote = F,row.names = F)



data <- FetchData(seu, 
                  vars = c('UMAP_1','UMAP_2', 'seurat_clusters'))
data$pt.col <- ifelse(data[['seurat_clusters']] =='23', "red", "#DFDFDF")
data$pt.size <- ifelse(data[['seurat_clusters']]=='2' , 0.2, 0.1)
table(data$pt.col)

ggplot(data, aes(UMAP_1, UMAP_2)) +
  geom_point(size=data$pt.size, color=data$pt.col) +
  theme_bw(base_size = 12) +
  #ggtitle(ct) +
  xlab('umap1') + ylab('umap2') +
  annotate("text",x=Inf,y=Inf,hjust=1.1,vjust=1.5,
           label='brain',color="red",size=8) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black")
  )+
  (DimPlot(seu,reduction='umap',group.by = 'seurat_clusters',label=T)+NoLegend())

seu@meta.data$celltype<-d[as.character(seu$seurat_clusters)]
table(seu$celltype)
pdf('./01merge/Dimplot_seu_all.lung.pdf',height=5,width=10)
(DimPlot(seu,group.by = 'seurat_clusters',label=T)+NoLegend())+ 
(DimPlot(seu,group.by = 'celltype',label=T)&ggsci::scale_color_d3('category20'))
dev.off()

qs::qsave(seu,'./01merge/seu.anno.qs')





##Dimplot Figure1E,F------
rm(list=ls())
seu<-qs::qread('./01merge/seu.anno.qs')
table(seu$celltype)
DimPlot(seu,label = T,group.by='celltype')
seu<-subset(seu,celltype=='Unknown',invert=T)
d <- c(
  "Tumor"="Tumor",#"#6A3D9A"
  "T_cells" = "T cells", #"#C99B7D"
  "T_NK_cells" = "T cells & NK_cells",# "#FF7F00"
  "B" = "B cells",#"#33A02C"
  "Mast_cells"="Mast cells", #"#e1a8b8"
  "Pericytes" = "Pericytes", #"#1F78B4"
  "EpiVascular"="Vascular endothelial cells", #"#68AB9F"
  "TAM"="TAMs & DC", #"#A6CEE3"
  "DC"="TAMs & DC",
  "Oligo_Astro"= "Oligodendrocytes" #"#B2DF8A" #Astro
)
seu$celltype_pl<-d[seu$celltype]
table(seu$celltype_pl %>% is.na())
color.lung<-c("#6A3D9A","#C99B7D","#FF7F00",
              "#33A02C",
              "#e1a8b8","#1F78B4","#68AB9F","#A6CEE3",
            #  "#C48244",
              "#B2DF8A")
names(color.lung)<-unique(d)
seu$celltype_pl<-factor(seu$celltype_pl,levels=c('Tumor','TAMs & DC','Oligodendrocytes',
                                                 'Vascular endothelial cells','Pericytes',
                                                 'T cells','T cells & NK_cells','B cells','Mast cells'))
pdf('./01merge/seu.all.dimplot.brain.pdf',height=8,width=8)
clusterCornerAxes(object = seu,reduction = 'umap',
                  clusterCol = 'celltype_pl',
                  noSplit = T,cellLabel = F,
                  keySize = 8,aspect.ratio = 1,
                  relLength = 0.5)& 
  scale_color_manual(values = color.lung) 
dev.off()

Idents(seu)<-factor(seu$celltype_pl,
                    levels=c('Tumor','TAMs & DC','Oligodendrocytes',
                             'Vascular endothelial cells','Pericytes',
                             'T cells','T cells & NK_cells','B cells','Mast cells'))
table(seu$celltype_pl)
pdf('./01merge/averagePlot.brain.pdf',height=8,width=8)
AverageHeatmap(object = seu,
               markerGene = c('KRT8', 'EPCAM',
                              'AIF1',#'TMEM119',
                              'CLEC10A',#'CD1C',
                              'GFAP',
                              'GAL3ST1',# 'GJC2',
                              'PECAM1', 'CLDN5',
                              'ACTA2', 'RGS5',
                             
                              'CD3D', 'CD2','GNLY',
                              'NKG7',
                              'IGKC','MZB1', 
                              'TPSAB1','TPSB2'),
               clusterAnnoName = F,
               myanCol = color.lung[levels(Idents(seu))],
               annoCol = TRUE,width = 8,height=8)
dev.off()
data.table::fwrite(table(seu$geo_accession,seu$diagnosis) %>% as.data.frame(),
                   './01merge/Brain.title.txt',sep='\t',
                   quote = F,row.names = F)
