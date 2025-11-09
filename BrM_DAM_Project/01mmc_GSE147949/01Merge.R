library(Seurat);library(tidyverse);library(Nebulosa);library(clustree);library(GEOquery)
file<-list.files('../download/GSE147949/GSE147942_RAW/')
sample <- unique(str_extract(file,'GSM\\d+'))
pattern <- ".*?(barcodes\\.tsv\\.gz|features\\.tsv\\.gz|matrix\\.mtx\\.gz)"
dir.create('../download/GSE147949/new/')
for (i in 1:6) {

  file0 <- file[str_detect(file,pattern)&str_detect(file,sample[i])]
  new_dir<-paste0('../download/GSE147949/new/',sample[i])
  dir.create(new_dir)
  file.copy(from=paste0('../download/GSE147949/GSE147942_RAW/',file0),
            to=new_dir)
  file0<-list.files(new_dir,full.names = T)
  file.rename(file0,to=paste0(new_dir,"/",sub(pattern, "\\1", file0)))

}
file0<-paste0('../download/GSE147949/new/',sample,'/')
seu.list <- pbapply::pblapply(file0, function(sn) {
  
  counts <- Read10X(sn)
  sn <- gsub("_", "-", sn) # 注意"_"在`CreateSeuratObject()`里有特殊的意义
  colnames(counts) <- paste(str_extract(sn,'GSM\\d+'), colnames(counts), sep = "_")
  seu <- CreateSeuratObject(counts = counts,meta.data = ,project = unique(str_extract(sn,'GSM\\d+'))[i],)
  return(seu)
})
seu <- base::Reduce(f = merge, x = seu.list)

rownames(seu@meta.data)

## metadata

file0 <- paste0('../download/GSE147949/GSE147942_RAW/',file[str_detect(file,'metadata')])
metadata<-pbapply::pblapply(file0, function(sn){
  if(str_detect(sn,'Con')){
    data <- data.table::fread(sn)
    data<-data %>% mutate(`231BR_clusters`=NA)
  }else{
    data <- data.table::fread(sn)
  }
})

x<- base::Reduce(f = rbind, x = metadata) 
table(str_extract(x$barcodes,'mouse_C\\d+|mouse_Met\\d+'))
x$barcodes <- str_replace_all(x$barcodes, c('mouse_C1' = 'GSM4450693', 
                                            'mouse_C2' = 'GSM4450694', 
                                            'mouse_C3' = 'GSM4450695', 
                                            'mouse_Met1' = 'GSM4450696', 
                                            'mouse_Met2' = 'GSM4450697', 
                                            'mouse_Met3' = 'GSM4450698'))

table(str_extract(x$barcodes,'GSM\\d+'))
meta_data<-seu@meta.data %>% 
  rownames_to_column('barcodes') %>% 
  mutate(barcodes=sapply(strsplit(barcodes,'-'),"[",1),
         geo_accession=orig.ident) %>%
  left_join(x[,c(1,6:12)]) %>% 
  mutate(orig.ident=str_replace_all(orig.ident, c('GSM4450693'='Con1', 
                                                  'GSM4450694'='Con2', 
                                                  'GSM4450695'='Con3', 
                                                  'GSM4450696'='Met1', 
                                                  'GSM4450697'='Met2', 
                                                  'GSM4450698'='Met3')),
         barcodes=paste0(barcodes,'-1')) %>% 
  column_to_rownames('barcodes')
seu@meta.data<-meta_data
table(str_detect(rownames(seu@assays$RNA@counts),'^grch38-MT'))

length(rownames(seu@assays$RNA@counts))
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^grch38-MT")

VlnPlot(seu, features = c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0) +
  geom_hline(yintercept = c(1000, 5000), linetype = "dashed", color = "blue")

VlnPlot(seu, features = c("percent.mt"), group.by = "orig.ident", pt.size = 0, log = T) +
  geom_hline(yintercept = c(10), linetype = "dashed", color = "blue")
seu.backup <- seu
dir.create('dataset3/')
dir.create('dataset3/tmp/')
qs::qsave(seu.backup, "./dataset3/tmp/01.seurat.init.qs")
table(seu$percent.mt<=10)
table(seu$nFeature_RNA<=9000)

seu <- subset(seu.backup, nFeature_RNA >= 200 & seu$nFeature_RNA<5000)
seu <- subset(seu, percent.mt <= 10)
seu.backup

qs::qsave(seu, "./dataset3/tmp/02-2.seurat.filtered.qs")

seu<-qs::qread('./dataset3/tmp/02-2.seurat.filtered.qs')
seu.list<-SplitObject(seu, split.by = "geo_accession")
for (i in 1:length(seu.list)) {
  seu.list[[i]] <- NormalizeData(seu.list[[i]], verbose = FALSE)
  seu.list[[i]] <- FindVariableFeatures(seu.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}
names(seu.list)



seu.anchors <- FindIntegrationAnchors(object.list = seu.list, dims = 1:20,
                                      anchor.features = 2000)

seu.integrated <- IntegrateData(anchorset = seu.anchors, dims = 1:20)
DefaultAssay(seu.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
seu.integrated <- RunPCA(seu.integrated, npcs = 30, verbose = FALSE)

ElbowPlot(seu.integrated, reduction = "pca", ndims = 30)
xx <- cumsum(seu.integrated[["pca"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.90) # 20 PCs解释了95%的方差，假设5%的方差来自于噪声
ndim = 20

## try a larger perplexity!
seu.integrated <- RunTSNE(seu.integrated, reduction = "pca", dims = 1:20, perplexity = 150)
seu.integrated <- RunUMAP(seu.integrated, reduction = "pca", dims = 1:20)
seu.integrated$type<-str_extract(seu.integrated$orig.ident,'[A-Z|a-z]{3}')
table(seu.integrated$type)

esolution = seq(0.2,1.2,0.1)
seu.integrated <- FindNeighbors(seu.integrated, reduction = "pca", dims = 1:20, k.param = 20)
seu.integrated <- FindClusters(seu.integrated, resolution = seq(0.2,1.2,0.1))

pdf('./dataset3/FigCellAnno/seu_integrated_umap_resolution.pdf',width=12,height=12)
DimPlot(seu.integrated, reduction = "umap", group.by = paste0("integrated_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu.integrated, reduction = "tsne", group.by = paste0("integrated_snn_res.", resolution), ncol = 3, label = T)
dev.off()

seu.integrated@meta.data$seurat_clusters <- seu.integrated@meta.data$integrated_snn_res.0.7
Idents(seu.integrated) <- "seurat_clusters"
DefaultAssay(seu.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
seu.integrated <- NormalizeData(seu.integrated, verbose = FALSE)

(DimPlot(seu.integrated, reduction = "tsne", group.by = "seurat_clusters", label = T)+NoLegend())+
  (DimPlot(seu.integrated, reduction = "umap", group.by = "seurat_clusters", label = T)) & ggsci::scale_color_d3("category20")

(DimPlot(seu.integrated, reduction = "umap", group.by = "seurat_clusters", label = T)+NoLegend())+
(DimPlot(seu.integrated, reduction = "umap", group.by = "combined_cell_type", label = T)) & ggsci::scale_color_d3("category20")

qs::qsave(seu.integrated, "./dataset3/tmp/02-3.seurat_process.seurat.qs")

##Divide MMC

seu.integrated.mmc<-subset(seu.integrated,seurat_clusters%in%c('8','10','16','12') & (combined_cell_type=='231BR'|is.na(combined_cell_type)),invert=T)
DimPlot(seu.integrated.mmc,label=T)

colnames(seu.integrated.mmc)
rownames(seu@assays$RNA@counts)
seu.mmc<-CreateSeuratObject(counts=seu@assays$RNA@counts[rownames(seu@assays$RNA@counts)[str_detect(rownames(seu@assays$RNA@counts),'mm10')],
                                                         colnames(seu.integrated.mmc@assays$RNA@counts)],
                            meta.data=seu.integrated.mmc@meta.data[colnames(seu.integrated.mmc@assays$RNA@counts),c(1:7,13,25)])
table(rownames(seu.mmc@meta.data)==colnames(seu.mmc@assays$RNA@counts)) 

table(is.na(sapply(strsplit(rownames(seu.mmc@assays$RNA@data),'---'),'[',2)))
anyDuplicated(sapply(strsplit(rownames(seu.mmc@assays$RNA@data),'---'),'[',2))
table(sapply(strsplit(rownames(seu.mmc@assays$RNA@counts),'---'),'[',2)==sub("^mm10---", "", rownames(seu.mmc@assays$RNA@data)))
rownames(seu.mmc@assays$RNA@counts)<-sapply(strsplit(rownames(seu.mmc@assays$RNA@counts),'---'),'[',2)
rownames(seu.mmc@assays$RNA@data)<-sapply(strsplit(rownames(seu.mmc@assays$RNA@data),'---'),'[',2)
seu.mmc[["RNA"]]@meta.features <- data.frame(row.names = rownames(seu.mmc[["RNA"]]))
seu.mmc$geo_accession

seu.mmc@assays$RNA@data
View(seu.mmc@meta.data)
DefaultAssay(seu.mmc)
seu.list<-SplitObject(seu.mmc, split.by = "geo_accession")
for (i in 1:length(seu.list)) {
  seu.list[[i]] <- NormalizeData(seu.list[[i]], verbose = FALSE)
  seu.list[[i]] <- FindVariableFeatures(seu.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}
names(seu.list)
seu.anchors <- FindIntegrationAnchors(object.list = seu.list, dims = 1:20,
                                      anchor.features = 2000)

seu.mmc.integrated <- IntegrateData(anchorset = seu.anchors, dims = 1:20)
DefaultAssay(seu.mmc.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
seu.mmc.integrated <- ScaleData(seu.mmc.integrated, verbose = FALSE)
seu.mmc.integrated <- RunPCA(seu.mmc.integrated, npcs = 30, verbose = FALSE)
ElbowPlot(seu.mmc.integrated, reduction = "pca", ndims = 30)
xx <- cumsum(seu.mmc.integrated[["pca"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.95) # 20 PCs解释了95%的方差，假设5%的方差来自于噪声
ndim = 20

## try a larger perplexity!
seu.mmc.integrated <- RunTSNE(seu.mmc.integrated, reduction = "pca", dims = 1:20, perplexity = 150)
seu.mmc.integrated <- RunUMAP(seu.mmc.integrated, reduction = "pca", dims = 1:20)
table(seu.mmc.integrated$type)
seu.mmc.integrated$type<-str_extract(seu.mmc.integrated$orig.ident,'[A-Z|a-z]{3}')

DimPlot(seu.mmc.integrated, reduction = "umap", group.by = "geo_accession")

resolution = seq(0.2,1.2,0.1)
seu.mmc.integrated <- FindNeighbors(seu.mmc.integrated, reduction = "pca", dims = 1:20, k.param = 20)
seu.mmc.integrated <- FindClusters(seu.mmc.integrated, resolution = seq(0.2,1.2,0.1))
pdf('./dataset3/FigCellAnno/seu_mmc_umap_resolution.pdf',width=12,height=12)
DimPlot(seu.mmc.integrated, reduction = "umap", group.by = paste0("integrated_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu.mmc.integrated, reduction = "tsne", group.by = paste0("integrated_snn_res.", resolution), ncol = 3, label = T)
dev.off()
(DimPlot(seu.mmc.integrated, reduction = "umap", group.by = "integrated_snn_res.0.9", label = T)+NoLegend())+
  (DimPlot(seu.mmc.integrated, reduction = "umap", group.by = "combined_cell_type", label = T) & ggsci::scale_color_d3("category20"))
(DimPlot(seu.mmc.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.9", label = T)+NoLegend())+
  (DimPlot(seu.mmc.integrated, reduction = "tsne", group.by = "combined_cell_type", label = T) & ggsci::scale_color_d3("category20"))
(DimPlot(seu.mmc.integrated, reduction = "umap", group.by = "integrated_snn_res.0.8", label = T)+NoLegend())+
  (DimPlot(seu.mmc.integrated, reduction = "umap", group.by = "combined_cell_type", label = T) & ggsci::scale_color_d3("category20"))
(DimPlot(seu.mmc.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.8", label = T)+NoLegend())+
  (DimPlot(seu.mmc.integrated, reduction = "tsne", group.by = "combined_cell_type", label = T) & ggsci::scale_color_d3("category20"))

seu.mmc.integrated@meta.data$seurat_clusters <- seu.mmc.integrated@meta.data$integrated_snn_res.0.8

#seu.mmc.integrated$integrated_snn_res.0.9==16
Idents(seu.mmc.integrated) <- "seurat_clusters"
DefaultAssay(seu.mmc.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
seu.mmc.integrated <- NormalizeData(seu.mmc.integrated, verbose = FALSE)

qs::qsave(seu.mmc.integrated, "./dataset3/tmp/02-4.seurat_mmc_integrated.qs")
##Markers====
seu<-qs::qread("./dataset3/tmp/02-4.seurat_mmc_integrated.qs")

Idents(seu)<-'seurat_clusters'
BRmarkers <- FindAllMarkers(seu, only.pos = TRUE)
BRmarkers <-BRmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(abs(avg_log2FC)))%>% 
  ungroup()
write.csv(BRmarkers,'./dataset3/231BRmusmarkers.csv',quote = F,row.names = F)
Idents(seu)<-'celltype_cluster'
BRmarkers_celltype <- FindAllMarkers(seu, only.pos = TRUE)
BRmarkers_celltype <-BRmarkers_celltype %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(abs(avg_log2FC)))%>%
  ungroup()
write.csv(BRmarkers_celltype,'./dataset3/231BRmusmarkers_celltype.csv',quote = F,row.names = F)


cluster.annot<-data.frame(cluster=c(0:25),
                          annotation=c('microglia',
                                       'microglia','astrocytes','microglia',
                                       'microglia','astrocytes','microglia','astrocytes','astrocytes','astrocytes',
                                       'oligo','astrocytes','astrocytes','microglia',
                                       'T_NK','myeloid_cells_DC ','microglia','astrocytes','matureDC','pericytes',
                                       'T_NK','Ependymal_cells','B','23','EpiVascular','Ependymal_cells'
                                       
                          ))


seu.mmc.integrated$celltype_cluster <- plyr::mapvalues(x = as.character(seu.mmc.integrated$seurat_clusters),
                                                       from = as.character(cluster.annot$cluster),
                                                       to = cluster.annot$annotation)



qs::qsave(seu.mmc.integrated, "./dataset3/tmp/02-4.seurat_mmc_anno.seurat.qs")



## seu.mmc dimplot------
rm(list=ls())
seu<-qs::qread('./PriorResults/tmp/02-4.seurat_mmc_anno.seurat.qs')
seu<-qs::qread('./V1/01seu.qs')

Idents(seu)<-'celltype_cluster'
DimPlot(seu,label=T)
SCpubr::do_DimPlot(seu)

##Dimplot------

library(scales)
show_col(color.use$meta.cluster)
length(unique(seu$celltype_cluster))
set.seed(20)
color.set<-color.use$meta.cluster[sample(length(unique(seu$celltype_cluster)))]
show_col(color.set)
color.set<-c('oligo'= "#B2DF8A", 
             'microglia'="#FB9A99",   
             'astrocytes'= "#62A3CB",
             'myeloid_cells_DC'="#A6CEE3",
             'Ependymal_cells'="#EF595A",
             'EpiVascular'="#68AB9F", 
             'pericytes'="#1F78B4",
             'T_NK'="#C99B7D", 
             'B'= "#33A02C", 
             'matureDC'="#C48244", 
             'T cells & NK cells'="#FF7F00")


d <- c(
  "microglia" = "Microglia",
  "astrocytes"="Astrocytes",
  "T_NK" = "T cells & NK cells",
  "B" = "B cells",
  "oligo" = "Oligodendrocytes",
  "pericytes" = "Pericytes",
  "EpiVascular"="Vascular endothelial cells",
  "myeloid_cells_DC "="BMDM",
  "matureDC"="Mature DC",
  "Ependymal_cells"="Ependymal cells" 
)
names(color.set)<-d[names(color.set)]

names(color.set)<-unique(seu$celltype_cluster)
Idents(seu)<-seu$celltype_cluster

table(seu$celltype_cluster)


seu$celltype_plot <- d[seu$celltype_cluster]
Idents(seu)<-seu$celltype_plot
names(color.set)<-d[names(color.set)]
table(seu$celltype_plot %>% is.na)
library(scRNAtoolVis)
seu$celltype_plot<-factor(seu$celltype_plot,levels=c('Microglia','BMDM','Mature DC',
                                                     'Astrocytes','Oligodendrocytes',
                                                     'Vascular endothelial cells','Pericytes',
                                                     'T cells & NK cells','B cells','Ependymal cells'))
pdf('./Figuresubmit/seu.mmc.BrTME.pdf',height=8,width=8)
clusterCornerAxes(object = seu,reduction = 'umap',
                  clusterCol = 'celltype_plot',
                  noSplit = T,cellLabel = F,
                  keySize = 8,aspect.ratio = 1,
                  relLength = 0.5)& 
  scale_color_manual(values = color.set) 

dev.off()

##cellabundant----
df<-data.frame(table(seu$celltype_plot,seu$type)) 
df<- df %>% group_by(Var2) %>% 
  mutate(percent=Freq/sum(Freq)*100) 

df$label = paste0(sprintf("%.2f", df$percent), "%")                

pdf('./submit/dataset3/01celltype.abundance.pdf',height=6,width=6)
ggplot(df,aes(Var2,percent,fill=Var1))+
  geom_bar(stat="identity",position = position_stack())+
  scale_fill_manual(values = color.set)+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="group",y="Celltype Percentage")+
  geom_text(aes(label=label),vjust=3,color="black")+
  
  theme_classic()+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))

dev.off()



pdf('./Figuresubmit/celltype.abundance.pie.pdf')
names <- table(subset(seu,type=='Met')$celltype_plot) %>% names()
ratio <- table(subset(seu,type=='Met')$celltype_plot) %>% as.numeric()
pielabel <- paste0(names," (", round(ratio/sum(ratio)*100,2), "%)")
pie(ratio, labels=pielabel,
    radius = 1.0,clockwise=T,
    main = "Met",col = color.set[names])
names <- table(subset(seu,type=='Con')$celltype_plot) %>% names()
ratio <- table(subset(seu,type=='Con')$celltype_plot) %>% as.numeric()
pielabel <- paste0(names," (", round(ratio/sum(ratio)*100,2), "%)")
pie(ratio, labels=pielabel,
    radius = 1.0,clockwise=T,
    main = "Con",col = color.set[names])
dev.off()
### abundant test----
seu@meta.data$sample<-rownames(seu@meta.data)

AbundanceTest <- function(cellmeta, celltype.col, sample.col, group.col) {
  count.matrix <- table(cellmeta[[celltype.col]], cellmeta[[sample.col]])
  count.matrix <- as.matrix(count.matrix)
  count.matrix.norm <- apply(count.matrix, 2, function(xx) xx / sum(xx))
  
  sample.meta <- cellmeta[, c(sample.col, group.col)] %>%
    distinct() %>%
    as.data.frame()
  rownames(sample.meta) <- sample.meta[[sample.col]]
  sample.meta[[sample.col]] <- NULL
  sample.meta <- sample.meta[colnames(count.matrix), ]
  
  celltypes <- sort(rownames(count.matrix.norm))
  
  group.levels <- levels(cellmeta[[group.col]])
  if (is.null(group.levels)) {
    group.levels <- unique(cellmeta[[group.col]])
  }
  
  test.df <- lapply(celltypes, function(ct) {
    c1 <- count.matrix.norm[ct, sample.meta == group.levels[1]]
    c2 <- count.matrix.norm[ct, sample.meta == group.levels[2]]
    res <- wilcox.test(c1, c2)
    data.frame(
      celltype = ct,
      fold.change = mean(c2) / mean(c1),#Met/con
      p.value = res$p.value,
      mean.perc = mean(c(c1,c2))
    )
  }) %>% do.call(rbind, .)
  return(test.df)
}
da.test <- AbundanceTest(cellmeta = seu@meta.data,
                         celltype.col = "celltype_cluster",
                         sample.col = "sample",
                         group.col = "type")
(da.test %>% filter(fold.change<1) %>% dplyr::select('celltype') )[,1]
(da.test %>% filter(fold.change>1) %>% dplyr::select('celltype') )[,1]

##averageheatmap-----
Idents(seu)<-seu$celltype_plot
table(Idents(seu))
Idents(seu) <-factor(seu$celltype_plot ,
                     levels=c('Microglia','BMDM','Mature DC',
                              'Astrocytes','Oligodendrocytes',
                              'Vascular endothelial cells','Pericytes',
                              'T cells & NK cells','B cells','Ependymal cells'))


pdf('./submit/dataset3/01seu.mmc.averageheatmap.pdf',height=6,width=6)
DefaultAssay(seu)<-'RNA'
AverageHeatmap(object = seu,
               markerGene = c('Tmem119', 'P2ry12', 'Cx3cr1',
                              'Gfap', 'Aldh1l1', 'Aldoc',
                              'Gal3st1', 'Gjc2','Mbp',
                              'Lyz2', 'Plac8','Cd163',
                              'Ccr7','Flt3',
                              'Pecam1', 'Vwf' ,'Cldn5',
                              'Acta2', 'Cspg4' ,'Rgs5',
                              'Ccdc153','Rarres2',
                              'Ms4a1', 'Cd79a','Igkc',
                              'Cd3d', 'Cd2', 'Nkg7'),
               clusterAnnoName = F,
               myanCol = color.set,
               annoCol = TRUE)
pdf('./Figuresubmit/seu.mmc.averageheatmap.pdf',height=6,width=6)
DefaultAssay(seu)<-'RNA'
AverageHeatmap(object = seu,
               markerGene = c('Tmem119', 'P2ry12', 
                              'Lyz2', 'Cd163',
                              'Ccr7','Flt3',
                              'Gfap', 'Aldoc',#'Aldh1l1'
                              'Gal3st1', 'Gjc2',
                              'Pecam1', 'Cldn5',
                              'Acta2', 'Rgs5',
                              'Cd3d', 'Nkg7',
                              'Igkc','Mzb1', 
                              'Ccdc153','Rarres2'),
               clusterAnnoName = F,
               #htCol = c("#276419","#F7F7F7","#8E0152"),
               myanCol = color.set[levels(Idents(seu))],
               annoCol = TRUE,width=8,height=8)
dev.off()

FeaturePlot(seu,  
            reduction = "umap",
            features = c(
              "Clec9a",
              "Itgam",
              "Ccr7", 
              "Flt3",
              "Clec10a",
              "Itga4"), 
            order = TRUE,ncol=3)






qs::qsave(color.set,file='./submit/dataset3/01coloruse.qs')
qs::qsave(seu,file='./submit/dataset3/01seu.qs')

markers_celltype<-FindAllMarkers(seu)
markers<-markers_celltype%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05 & avg_log2FC>1.0) %>% 
  arrange(cluster,desc(avg_log2FC))%>%
  ungroup()

data.table::fwrite(markers,
                   './submit/dataset3/01seu.mmc_celltype_markers.txt',sep='\t',
                   quote = F,row.names = F)






