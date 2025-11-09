rm(list=ls())
gc()
library(Seurat);library(tidyverse);library(Nebulosa);
library(clustree);library(GEOquery);library(scCustomize);library(scRNAtoolVis)
library(UCell);library(clusterProfiler);library(harmony)
##seu.TAM-----
rm(list=ls())
gc()
seu<-qs::qread('./01seuRegular/seu.anno.qs')
seu.TAM<-subset(seu,celltype=='Myeloid cells')
seu.TAM <- NormalizeData(seu.TAM)
seu.TAM <- FindVariableFeatures(seu.TAM, selection.method = "vst", nfeatures = 2000)
seu.TAM<- ScaleData(seu.TAM, features = rownames(seu.TAM))
seu.TAM<- RunPCA(seu.TAM, features = VariableFeatures(object = seu.TAM),reduction.name = "pca")

table(seu.TAM$title)
seu.TAM <- RunHarmony(seu.TAM ,
                      reduction = "pca",
                      group.by.vars = "title",
                      reduction.save = "harmony")
ElbowPlot(seu.TAM, reduction = "harmony", ndims = 30)
xx <- cumsum(seu.TAM[["harmony"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.90) 
ndim = 30
seu.TAM  <- RunUMAP(seu.TAM , 
                    reduction = "harmony", 
                    dims = 1:30, 
                    reduction.name = "umap")
seu.TAM  <- RunUMAP(seu.TAM , 
                    reduction = "pca", 
                    dims = 1:30, 
                    reduction.name = "umap_pca")
seu.TAM  <- RunUMAP(seu.TAM , 
                    reduction = "harmony", 
                    dims = 1:30, 
                    n.neighbors = 80,# 100
                    min.dist = 0.1,  
                    metric = "correlation", 
                    reduction.name = "umap_add")

seu.TAM <- RunTSNE(seu.TAM, reduction = "harmony", 
                   dims = 1:30,
                   reduction.name = "tsne", 
                   perplexity = 150)
DimPlot(seu.TAM, reduction = "umap", split.by = "title")

resolution = seq(0.2,2,0.1)
seu.TAM <- FindNeighbors(seu.TAM, reduction = "harmony", dims = 1:30, k.param = 20)
seu.TAM <- FindClusters(seu.TAM, resolution = seq(0.2,2,0.1))
dir.create('./02TAM/')
pdf('./02TAM/seu.TAM_umap_resolution.pdf',width=18,height=18)
DimPlot(seu.TAM, reduction = "umap", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu.TAM, reduction = "tsne", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
dev.off()

(DimPlot(seu.TAM, reduction = "umap", group.by = 'Cell_subtype',label = T)) +
  (DimPlot(seu.TAM, reduction = "umap", group.by='Cell_type.refined',label = T)+NoLegend()) & ggsci::scale_color_d3("category20")



seu.TAM$seurat_clusters <- as.character(seu.TAM$RNA_snn_res.0.6)
seu.TAM$TAM_cluster<-paste0('TAM.',seu.TAM$seurat_clusters)
Idents(seu.TAM)<-'TAM_cluster'
DefaultAssay(seu.TAM)<-'RNA'
##Findmarkers------

Idents(seu.TAM)<-seu.TAM$TAM_cluster
DefaultAssay(seu.TAM)<-'RNA'
TAM_markers<-FindAllMarkers(seu.TAM)
TAM_markers <-TAM_markers%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(avg_log2FC))%>%
  ungroup()


seu.TAM$sample<-rownames(seu.TAM@meta.data)
seu$sample<-rownames(seu@meta.data)
table(seu.TAM$TAM_cluster)
seu$TAM_cluster<-NULL
seu@meta.data<-seu@meta.data %>% 
  left_join(seu.TAM@meta.data[,c('sample','TAM_cluster')]) %>% 
  mutate(TAM_cluster=ifelse(is.na(TAM_cluster),
                            'others',TAM_cluster))

rownames(seu@meta.data)<-seu@meta.data$sample
Idents(seu)<-seu$TAM_cluster
table(seu$TAM_cluster)
table(Idents(seu))

DefaultAssay(seu)<-'RNA'
TAM_markers2<-FindAllMarkers(seu)
TAM_markers2 <-TAM_markers2%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(avg_log2FC))%>%
  ungroup() %>% 
  filter(!(cluster=='others'))

data.table::fwrite(TAM_markers2,
                   './02TAM/TAM_markers_all_seurat_clusters.txt',sep='\t',
                   quote = F,row.names = F)
data.table::fwrite(TAM_markers,
                   './02TAM/TAM_markers_seurat_clusters.txt',sep='\t',
                   quote = F,row.names = F)
table(seu.TAM$TAM_cluster)
seu.TAM<-subset(seu.TAM,TAM_cluster%in%c('TAM.13','TAM.12'),invert=T)
d<-c('TAM.0'='APOE.MDM',
     'TAM.1'='APOE.MDM.DAM',
     'TAM.2'='APOE.MDM',
     'TAM.3'='CD14.Monocytes',
     'TAM.4'='APOE.MDM.DAM',
     'TAM.5'='Mac_MT1H',#Mac_MT1H
     'TAM.6'='DC',
     'TAM.7'='Recycling.MDM',
     'TAM.8'='S100A8.MDM',
     'TAM.9'='APOE.MDM', #混Microglia
     'TAM.10'='Stressed_Mac',
     'TAM.11'='DAM.Microglia'
)
#TAM.13 混tumor    TAM.12 混Tcell
seu.TAM$celltype<-d[seu.TAM$TAM_cluster]
##Dimplot-----
seu.TAM<-qs::qread('./02TAM/seu.TAM.qs')
table(seu.TAM$celltype)
celltype_pl<-d<-c('DAM.Microglia'='DAM',
                  'APOE.MDM'='APOE.Mac',
                  'APOE.MDM.DAM'='DAM-like.Mac',
                  'CD14.Monocytes'='CD14.Mono',
                  'DC'='DC',
                  'Mac_MT1H'='MT1H.Mac',
                  'S100A8.MDM'='S100A8.Mac',
                  'Recycling.MDM'='Prolif.Mac',
                  'Stressed_Mac'='Stressed.Mac')
table(seu.TAM$celltype_pl %>% is.na())
seu.TAM$celltype_pl<-celltype_pl[seu.TAM$celltype]
saveRDS(seu.TAM,'./FigureSubmit/seu.TAM.celltype.pl.rds')
color.TAM_hgc<-readRDS('~/Brain/brian/DAM/color.TAM_hgc.rds')
dir.create('./FigureSubmit/')
pdf('./FigureSubmit/TAM_Dimplot_pl.pdf',width=5.5,height=4)
DimPlot(seu.TAM,
        reduction = 'umap',
        group.by = 'celltype_pl')+theme_bw()+
  labs(color = "Cell types")+
  ggtitle(NULL) + 
  theme(panel.grid  = element_blank(),
        aspect.ratio = 1,
        legend.text = element_text(size=10),  # 您可以根据需要更改字体大小
        legend.title = element_text(size=10),
        plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
        axis.text = element_text(color = "black",size=12),
        legend.position = "right")& 
  scale_color_manual(values = color.TAM_hgc[unique(seu.TAM$celltype_pl)]) 
dev.off()

## Ucell BMDM[Bowman et.al] ------
MG_dataset<-qs::qread('/Users/huzixin/Brain/brian/Geneset/MG_dataset_forhuman.qs')
MG_BMDM<-MG_dataset %>% filter(category%in%c("Core_BMDM","Core_MG","TAM_BMDM","TAM_MG"))
MG_BMDM<-split(MG_BMDM$gene, MG_BMDM$category)
library(UCell)
DefaultAssay(seu.TAM)<-'RNA'
rownames(seu.TAM)
seu.TAM<- AddModuleScore_UCell(seu.TAM,
                               features = MG_BMDM,
                               assay = 'RNA',
                               slot='counts',
                               name = "_UCell")


marker_genes <- c( "Core_BMDM_UCell",
                   "Core_MG_UCell" ,
                   "TAM_BMDM_UCell",
                   "TAM_MG_UCell")              
data<-seu.TAM@meta.data[,marker_genes]
colnames(data)<-names(MG_BMDM)
data<-t(data)
colnames(data)<-colnames(seu.TAM)
seu.TAM[['MG_BMDM']]<-Seurat::CreateAssayObject(data = data)

DefaultAssay(seu.TAM)<-'MG_BMDM'
library(scCustomize)
FeaturePlot_scCustom(seu.TAM,reduction = 'umap',
                     features=c('Core-MG','Core-BMDM','TAM-MG','TAM-BMDM'))
FeaturePlot(seu.TAM,reduction = 'umap',
            features=c('Core-MG','Core-BMDM','TAM-MG','TAM-BMDM'))
library(scRNAtoolVis)

pdf('./02TAM/MG_BMDM_avg.Recluster.pdf',height=2.5,width=6)
Idents(seu.TAM)<-factor(seu.TAM$TAM_cluster,levels=paste0('TAM.',0:11))
Idents(seu.TAM)<-seu.TAM$celltype
AverageHeatmap(object = seu.TAM ,
               assays = 'MG_BMDM',
               markerGene =c('Core-MG','Core-BMDM','TAM-MG','TAM-BMDM'),
               clusterAnnoName = F,
               # myanCol = color.TAM2[levels(Idents(seu.TAM))],
               annoCol = F)
dev.off()
##BMDM average plotFigure3B====
pdf('./Figuresubmit/MG_BMDM_avg.Recluster.pdf',height=2.5,width=6)
Idents(seu.TAM)<-factor(seu.TAM$celltype_pl,levels=c('DAM','CD14.Mono','S100A8.Mac',
                                                     'APOE.Mac','DAM-like.Mac',
                                                     'Prolif.Mac',
                                                     'Stressed.Mac',
                                                     'MT1H.Mac',
                                                     'DC'))
DefaultAssay(seu.TAM)<-'MG_BMDM'
AverageHeatmap(object = seu.TAM ,
               assays = 'MG_BMDM',
               markerGene =c('Core-MG','Core-BMDM','TAM-MG','TAM-BMDM'),
               clusterAnnoName = F,height=2,width=6,
               myanCol = color.TAM_hgc[levels(Idents(seu.TAM))],
               annoCol = T)
dev.off()




##Ucell DAM[Gan et.al]----

library(UCell)
MG_dataset<-data.table::fread('/Users/huzixin/Brain/brian/Geneset/Gan_et_al_MGdataet1_hg.txt')
MG_dataset$Category[str_detect(MG_dataset$Category,'Homeostasis')]<-'Homeostasis'
data1<-MG_dataset[MG_dataset$Category=='Homeostasis',] %>% distinct(Gene,.keep_all = T)

MG_dataset<-rbind(data1,MG_dataset[!(MG_dataset$Category=='Homeostasis'),])

table(MG_dataset$Category)
MG_list<-split(MG_dataset$Gene, MG_dataset$Category)
DefaultAssay(seu.TAM)<-'RNA'
seu.TAM<- AddModuleScore_UCell(seu.TAM, 
                               features = MG_list,
                               assay = 'RNA',
                               slot='counts',
                               name = "_UCell")
marker_genes <- paste0(names(MG_list),'_UCell')
data<-seu.TAM@meta.data[,marker_genes]
colnames(data)<-names(MG_list)
data<-t(data)
colnames(data)<-colnames(seu.TAM)
seu.TAM[['DAM']]<-Seurat::CreateAssayObject(data = data)
DefaultAssay(seu.TAM)<-'DAM'
library(RColorBrewer)
col <- colorRampPalette(brewer.pal(11, "PiYG"))(25)



pdf('./02TAM/DAM_TAM_recluster_avgexp.pdf',height=3,width=5)
Idents(seu.TAM)<-factor(seu.TAM$TAM_cluster,levels=paste0('TAM.',0:11))
Idents(seu.TAM)<-seu.TAM$celltype
AverageHeatmap(object = seu.TAM,
               assays = "DAM",
               markerGene = c("DAM-core-program", 
                              "Stage1-DAM-upregulated",
                              #"Stage1-DAM-downregulated",
                              "Stage2-DAM","Homeostasis" ),
               annoCol = F,
               #myanCol = color.TAM2[levels(Idents(seu.TAM_TNBC))],
               clusterAnnoName = T)

dev.off()


##Ucell TAM[Katrina T et.al]------
MG_dataset<-qs::qread('/Users/huzixin/Brain/brian/Geneset/MG_dataset_forhuman.qs')

TAM.features<-MG_dataset %>% 
  filter(category%in%c('APC_Topic_Score',
                       'IFN_Response_Topic_Score',
                       'Secretory_Topic_Score'))

TAM.features<-split(TAM.features$gene, TAM.features$category)
DefaultAssay(seu.TAM)<-'RNA'
seu.TAM<- AddModuleScore_UCell(seu.TAM,features = TAM.features,name = "_UCell")

marker_genes <-c("APC_Topic_Score_UCell",
                 "IFN_Response_Topic_Score_UCell",
                 "Secretory_Topic_Score_UCell" ) 

data<-seu.TAM@meta.data[,marker_genes]
colnames(data)<-names(TAM.features)
data<-t(data)
colnames(data)<-colnames(seu.TAM)
seu.TAM[['TAM.features']]<-Seurat::CreateAssayObject(data = data)
col <- colorRampPalette(brewer.pal(11, "PiYG"))(25)
DefaultAssay(seu.TAM)<-'TAM.features'

FeaturePlot_scCustom(seurat_object = seu.TAM, 
                     features = c('APC-Topic-Score','IFN-Response-Topic-Score',
                                  'Secretory-Topic-Score'),
                     colors_use = rev(col))&
  theme(legend.title = element_blank())
pdf('./02TAM/TAM.features_avgexp.pdf',height=3,width=6)
Idents(seu.TAM)<-factor(seu.TAM$TAM_cluster,levels=paste0('TAM.',0:11))
Idents(seu.TAM)<-seu.TAM$celltype
AverageHeatmap(object = seu.TAM,
               assays = "TAM.features",
               markerGene = c('APC-Topic-Score','IFN-Response-Topic-Score',
                              'Secretory-Topic-Score'),
               #myanCol = color.TAM2[levels(Idents(seu.TAM_TNBC))],
               annoCol =F,
               clusterAnnoName = T)
dev.off()
##DAM&TAM.featuresFigure4A------
pdf('./Figuresubmit/TAM.features.pdf',height=2.5,width=6)
Idents(seu.TAM)<-factor(seu.TAM$celltype_pl,levels=c('DAM','CD14.Mono','S100A8.Mac',
                                                     'APOE.Mac','DAM-like.Mac',
                                                     'Prolif.Mac',
                                                     'Stressed.Mac',
                                                     'MT1H.Mac',
                                                     'DC'))
AverageHeatmap(object = seu.TAM ,
               assays = 'DAM',
               markerGene =c("DAM-core-program", 
                             "Stage1-DAM-upregulated",
                             "Stage2-DAM"),
               clusterAnnoName = F,
               myanCol = color.TAM_hgc[levels(Idents(seu.TAM))],
               annoCol = T,height=2,width=6)
dev.off()
pdf('./Figuresubmit/DAM.pdf',height=2.5,width=6)

AverageHeatmap(object = seu.TAM ,
               assays = 'TAM.features',
               markerGene =c('APC-Topic-Score','IFN-Response-Topic-Score',
                             'Secretory-Topic-Score'),
               clusterAnnoName = F,
               myanCol = color.TAM_hgc[levels(Idents(seu.TAM))],
               annoCol = T,height=2,width=6)
dev.off()



##提取DAM------
DefaultAssay(seu.TAM)<-'DAM'
geneExp <-FetchData(
  object = seu.TAM,
  vars =c("DAM-core-program", 
          "Stage1-DAM-upregulated",
          "Stage2-DAM"),
  slot = 'data'
)
geneExp$id <- seu.TAM$celltype_pl
PercentAbove <- utils::getFromNamespace("PercentAbove", "Seurat")
data.plot <- lapply(
  X = unique(geneExp$id),
  FUN = function(ident) {
    data.use <- geneExp[geneExp$id == ident, 1:(ncol(geneExp) - 1), drop = FALSE]
    avg.exp <- apply(
      data.use,2,
      function(x) {mean(expm1(x))}
    )
    pct.exp <- apply(data.use,2, PercentAbove, threshold = 0)
    
    res <- data.frame(id = ident, avg.exp = avg.exp, pct.exp = pct.exp * 100)
    res$gene <- rownames(res)
    return(res)
  }
) %>%
  do.call("rbind", .) %>%
  data.frame()
data.plot<-purrr::map_df(unique(data.plot$gene), function(x) {
  tmp <- data.plot %>%
    dplyr::filter(gene == x)
  avg.exp.scale <- tmp %>%
    dplyr::select(avg.exp) %>%
    scale(.) %>%
    Seurat::MinMax(., min = -2.5, max = 2.5)
  tmp$avg.exp.scaled <- as.numeric(avg.exp.scale)
  return(tmp)}) 


data.plot$gene <- factor(data.plot$gene, levels =c("DAM-core-program", 
                                                   "Stage1-DAM-upregulated",
                                                   "Stage2-DAM" ))
data.plot$id<-factor(data.plot$id, 
                     levels = rev(c('DAM','CD14.Mono','S100A8.Mac',
                                    'APOE.Mac','DAM-like.Mac',
                                    'Prolif.Mac',
                                    'Stressed.Mac',
                                    'MT1H.Mac',
                                    'DC')))
data.plot.wide<-data.plot[,c(1,4,5)] %>% 
  pivot_wider(names_from = 'id',values_from = 'avg.exp.scaled') %>% 
  column_to_rownames('gene')
saveRDS(data.plot.wide,'./Figuresubmit/DAM.wide.Lung3.rds')

Idents(seu.TAM)<-seu.TAM$celltype_pl

AverageExp<-AverageExpression(seu.TAM,
                              assays = 'DAM',
                              features=c("DAM-core-program", 
                                         "Stage1-DAM-upregulated",
                                         "Stage2-DAM" ))$DAM

p <- t(AverageExp) %>% scale() %>% t()

saveRDS(p,'./Figuresubmit/DAM.wide.Lung3_v2.rds')

AverageExp<-AverageExpression(seu.TAM,
                              assays = 'TAM.features',
                              features=c(c('APC-Topic-Score','IFN-Response-Topic-Score',
                                           'Secretory-Topic-Score')))$TAM.features

p <- t(AverageExp) %>% scale() %>% t()
saveRDS(p,'./Figuresubmit/TAM.features.wide.Lung3.rds')
##TAM markers plotFigure3C------
DefaultAssay(seu.TAM)<-'RNA'
TAM_features<-c( "TMEM119","P2RY12", "CX3CR1", 
                 "CSF1R","OLFML3","GPR34",
                 "ITGA4","PLAC8",#'CCR2','CD14',
                 'S100A8',"CXCL8",'FCN1',
                 "APOE","C1QB","SPP1",
                 "MRC1","CD163","SIGLEC1",
                 'STMN1','TOP2A',
                 'HSPA1B', 'HSPB1',
                 'MT1H','MT1M',
                 'CST3','CD1C',
                 #'CAMP', 'S100A9', 
                 "HLA-DPB1",'CD74')
TAM_features<-data.frame(gene=TAM_features,
                         cluster=c(rep('Homeo.MG',6),
                                   rep('BMDM',2),
                                   rep('Mono&S100A8.Mac',3),
                                   rep('APOE.Mac',3),
                                   rep('Mf',3),
                                   rep('Prolif.MG',2),
                                   rep('Stressed.Mac',2),
                                   rep('MT1H.Mac',2),
                                   rep('DC',2),
                                   rep('AP',2)))

geneExp <-FetchData(
  object = seu.TAM,
  vars = TAM_features$gene,
  slot = 'data'
)
geneExp$id <- seu.TAM$celltype_pl
PercentAbove <- utils::getFromNamespace("PercentAbove", "Seurat")
data.plot <- lapply(
  X = unique(geneExp$id),
  FUN = function(ident) {
    data.use <- geneExp[geneExp$id == ident, 1:(ncol(geneExp) - 1), drop = FALSE]
    avg.exp <- apply(
      data.use,2,
      function(x) {mean(expm1(x))}
    )
    pct.exp <- apply(data.use,2, PercentAbove, threshold = 0)
    
    res <- data.frame(id = ident, avg.exp = avg.exp, pct.exp = pct.exp * 100)
    res$gene <- rownames(res)
    return(res)
  }
) %>%
  do.call("rbind", .) %>%
  data.frame()
data.plot<-purrr::map_df(unique(data.plot$gene), function(x) {
  tmp <- data.plot %>%
    dplyr::filter(gene == x)
  avg.exp.scale <- tmp %>%
    dplyr::select(avg.exp) %>%
    scale(.) %>%
    Seurat::MinMax(., min = -2.5, max = 2.5)
  tmp$avg.exp.scaled <- as.numeric(avg.exp.scale)
  return(tmp)}) 

table(seu.TAM$celltype_pl)
table(data.plot$id)
data.plot$gene <- factor(data.plot$gene, levels = unique(data.plot$gene))
data.plot$id<-factor(data.plot$id, 
                     levels = rev(c('DAM','CD14.Mono','S100A8.Mac',
                                    'APOE.Mac','DAM-like.Mac',
                                    'Prolif.Mac',
                                    'Stressed.Mac',
                                    'MT1H.Mac',
                                    'DC'
                     )))

pdf('./Figuresubmit/MG_features.pdf',height=5,width=10)
ggplot(data.plot,aes(x = gene, y = id)) +
  geom_point(aes(fill = avg.exp.scaled, size = avg.exp.scaled),
             color = "black", shape = 21) +
  ggplot2::scale_size(range = c(1, 6))+
  theme_bw(base_size = 14) +
  xlab("") +ylab("") +
  coord_fixed(clip = "off") +
  theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
        axis.text = element_text(color = "black"),
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0),
        legend.position = "right")+
  guides(fill = guide_colorbar(title = "Mean expression \n in group",
                               title.position = "top",
                               title.hjust = 0.5,
                               barwidth = unit(2, "cm"),
                               frame.colour = "black",
                               frame.linewidth = 0.5,
                               ticks.colour = "black"))+
  guides(size = guide_legend(title = "Mean expression \n in group",
                             title.position = "top",
                             title.hjust = 0.5,
                             label.position = "bottom",
                             override.aes = list(color = "black",fill = "grey50"),
                             keywidth = unit(0.3, "cm")))+
  scale_fill_gradient2(low = '#0099CC',mid = "white",high = '#CC3333')
dev.off()
pdf('./Figuresubmit/MG_features_add.pdf',height=5,width=10)
jjDotPlot(object = seu.TAM,
          markerGene  = TAM_features, 
          anno = T,
          plot.margin = c(1,1,1,1), lwd = 0.5,textSize = 14,
          dot.col = c('#0099CC',"white",'#CC3333'),
          id = 'celltype_pl',ytree = F)
dev.off()

##Findmarkers======
rm(list=ls());gc()
seu.TAM<-qs::qread('./02TAM/seu.TAM.qs')
seu.TAM$celltype.new<-seu.TAM$celltype
seu<- qs::qread('./01seuRegular/seu.anno.qs')
seu$sample<-rownames(seu@meta.data)
seu.TAM$sample<-rownames(seu.TAM@meta.data)

seu@meta.data<-seu@meta.data %>% 
  left_join(seu.TAM@meta.data[,c('sample','celltype.new')]) %>% 
  mutate(celltype.new=ifelse(is.na(celltype.new),
                             'others',celltype.new))

table(seu$celltype.new)
DefaultAssay(seu.TAM)<-'RNA'
Idents(seu.TAM)<-seu.TAM$celltype.new
TAM_markers<-FindAllMarkers(seu.TAM)
TAM_markers <-TAM_markers%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(avg_log2FC))%>%
  ungroup()
Idents(seu)<-seu$celltype.new
DefaultAssay(seu)<-'RNA'
TAM_markers2<-FindAllMarkers(seu)
TAM_markers2 <-TAM_markers2%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(avg_log2FC))%>%
  ungroup() %>% 
  filter(!(cluster=='others'))

data.table::fwrite(TAM_markers,
                   './02TAM/TAMcelltype_markers.txt',sep='\t',
                   quote = F,row.names = F)
data.table::fwrite(TAM_markers2,
                   './02TAM/TAMcelltype_cluster_markers_all.txt',sep='\t',
                   quote = F,row.names = F)


##Enrichment-----
seu.TAM<-qs::qread('./02TAM/seu.TAM.qs')
TAM_markers<-data.table::fread('./02TAM/TAMcelltype_cluster_markers_all.txt')
library(clusterProfiler);library(org.Hs.eg.db)
er.list<-lapply(unique(TAM_markers$cluster), function(ct){
  data<-TAM_markers%>% 
    filter(cluster==ct) %>% 
    arrange(desc(avg_log2FC))
  if(!(nrow(data)==0)){
    geneList <- data$avg_log2FC
    names(geneList) = data$gene
    geneList = sort(geneList, decreasing = TRUE)
    ego <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 keyType = 'SYMBOL',
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
    ego<-ego@result
    if (nrow(ego)!=0) {
      ego$celltype<-ct
      return(ego)
    }
  }
})

filtered_list <- er.list[!sapply(er.list, is.null)]
er.df_TAM<- purrr::reduce(filtered_list,rbind) %>% 
  arrange(celltype,desc(NES))
data.table::fwrite(er.df_TAM,
                   './02TAM/TAMcelltype_enrichment.txt',sep='\t',
                   quote = F,row.names = F)


View(er.df_TAM[er.df_TAM$celltype=='DAM.Microglia',])
View(er.df_TAM[er.df_TAM$celltype=='APOE.MDM.DAM',])

##DAM genesHotmap Figure4E----

DAM_features<-c('APOE','TREM2','TYROBP',
                'SPP1','OLR1','CCL3','CCL4',
                'CD83','CSF1R','LILRB4', 'TMIGD3', 'GPR34','VSIG4','MS4A7','FCGR1A',
                'CD163','DAB2',
                'CTSB','CTSS','CST3',
                'PSAP','FTL',"GRN",'CAPG','FCGR2A','FCER1G',
                'MS4A4A','MS4A6A')
DefaultAssay(seu.TAM)<-'RNA'
geneExp <-FetchData(
  object = subset(seu.TAM,celltype_pl%in%c('DAM','DAM-like.Mac','APOE.Mac')),
  vars = DAM_features,
  slot = 'data'
)
geneExp$id <- subset(seu.TAM,celltype_pl%in%c('DAM','DAM-like.Mac','APOE.Mac'))$celltype_pl
PercentAbove <- utils::getFromNamespace("PercentAbove", "Seurat")
data.plot <- lapply(
  X = unique(geneExp$id),
  FUN = function(ident) {
    data.use <- geneExp[geneExp$id == ident, 1:(ncol(geneExp) - 1), drop = FALSE]
    avg.exp <- apply(
      data.use,2,
      function(x) {mean(expm1(x))}
    )
    pct.exp <- apply(data.use,2, PercentAbove, threshold = 0)
    
    res <- data.frame(id = ident, avg.exp = avg.exp, pct.exp = pct.exp * 100)
    res$gene <- rownames(res)
    return(res)
  }
) %>%
  do.call("rbind", .) %>%
  data.frame()
data.plot<-purrr::map_df(unique(data.plot$gene), function(x) {
  tmp <- data.plot %>%
    dplyr::filter(gene == x)
  avg.exp.scale <- tmp %>%
    dplyr::select(avg.exp) %>%
    scale(.) %>%
    Seurat::MinMax(., min = -2.5, max = 2.5)
  tmp$avg.exp.scaled <- as.numeric(avg.exp.scale)
  return(tmp)}) 


data.plot$gene <- factor(data.plot$gene, levels = unique(data.plot$gene))
data.plot$id<-factor(data.plot$id, 
                     levels = rev(c('DAM','APOE.Mac','DAM-like.Mac')))

dir.create('./Figuresubmit/')
pdf('./Figuresubmit/DAM_gene.pdf',height=3,width=9)
ggplot(data.plot,aes(x = gene, y = id)) +
  geom_tile(aes(fill = avg.exp.scaled),
            color = "black") +
  theme_bw(base_size = 14) +
  xlab("") +ylab("") +
  coord_fixed(clip = "off") +
  theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
        axis.text = element_text(color = "black"),
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0),
        legend.position = "right")+
  coord_fixed(clip = "off") +
  guides(fill = guide_colorbar(title = "Mean expression \n in group",
                               title.position = "top",
                               title.hjust = 0.5,
                               barwidth = unit(4, "cm"),
                               frame.colour = "black",
                               frame.linewidth = 0.5,
                               ticks.colour = "black"))+
  scale_fill_gradient2(low = '#0099CC',mid = "white",high = '#CC3333')
dev.off()

data.plot.wide<-data.plot[,c(1,4,5)] %>% 
  pivot_wider(names_from = 'id',values_from = 'avg.exp.scaled')
saveRDS(data.plot.wide,'./Figuresubmit/data.plot.wide.Lung3.rds')
saveRDS(data.plot.wide,'./Figuresubmit/data.plot.wide.Lung3_2.rds')

##DAM cellmarkersFigure4C,D ------
TAM_markers_BM<-data.table::fread('./02TAM/TAMcelltype_cluster_markers_all.txt')
table(TAM_markers_BM$cluster)
DAM_markers<-readRDS('~/Brain/brian/DAM/KEY_Gene_DAM_ALL.rds')
table(TAM_markers$cluster)
data<-TAM_markers%>% 
  filter(cluster%in%c('DAM.Microglia')) %>% 
  arrange(desc(avg_log2FC)) %>% 
  filter(avg_log2FC>=0)
data<-data%>% 
  mutate(Gene=1:nrow(data))
data <- head(data, n=200)

data.label<-data[data$gene%in%DAM_markers,]

data$pt.col <- ifelse(data$gene%in%data.label$gene, "#B294C7", "#BECEE3")
pdf('./Figuresubmit/DAMgene.pdf',height=7,width=8)
ggplot(data, aes(Gene, avg_log2FC)) +
  geom_point(size=3, color=data$pt.col) +
  ggrepel::geom_text_repel(inherit.aes = FALSE, max.overlaps = Inf,
                           data = data.label, aes(Gene,avg_log2FC, 
                                                  label=gene), size=6) +
  theme_bw(base_size = 12) +
  theme(panel.grid  = element_blank(),
        aspect.ratio = 1,
        legend.text = element_text(size=12),  
        plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
        axis.text = element_text(color = "black",size=14),
        axis.title = element_text(color = "black",size=14),
        legend.position = "right",
        legend.title = element_blank())+
  labs(color = NULL) 
dev.off()
##DAM-like.Mac Supplementary FigureS5A------
TAM_markers<-data.table::fread('./06trajectory.new/TAMcelltypetra_markers_all.txt')
DAM.Mac_Gene<-readRDS('~/Brain/brian/DAM/Key_Gene_DAM_MDM_Tra.rds')
data<-TAM_markers%>% 
  filter(cluster%in%c('DAM-like.Mac.Terminal')) %>% 
  arrange(cluster,desc(avg_log2FC)) %>% 
  filter(avg_log2FC>=0)
data<-data%>% 
  mutate(Gene=1:nrow(data))
data <- head(data, n=200)

data.label<-data[data$gene%in%DAM.Mac_Gene,] 
data$pt.col <- ifelse(data$gene%in%data.label$gene, "#98D277", "#BECEE3")
pdf('./Figuresubmit/DAM-like.Mac.gene.pdf',height=7,width=8)
ggplot(data, aes(Gene, avg_log2FC)) +
  geom_point(size=3, color=data$pt.col) +
  ggrepel::geom_text_repel(inherit.aes = FALSE, max.overlaps = Inf,
                           data = data.label, aes(Gene,avg_log2FC, 
                                                  label=gene), size=6) +
  theme_bw(base_size = 12) +
  theme(panel.grid  = element_blank(),
        aspect.ratio = 1,
        legend.text = element_text(size=12),  
        plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
        axis.text = element_text(color = "black",size=14),
        axis.title = element_text(color = "black",size=14),
        legend.position = "right",
        legend.title = element_blank())+
  labs(color = NULL) 
dev.off()

##HALLMARK_INTERFERON_GAMMA_RESPONSEFigure4F----
TAM_markers<-data.table::fread('./02TAM/TAMcelltype_cluster_markers_all.txt')
geneList <- TAM_markers$avg_log2FC[TAM_markers$cluster=='DAM.Microglia']
names(geneList) = TAM_markers$gene[TAM_markers$cluster=='DAM.Microglia']

geneList <- TAM_markers$avg_log2FC[TAM_markers$cluster=='APOE.MDM.DAM']
names(geneList) = TAM_markers$gene[TAM_markers$cluster=='APOE.MDM.DAM']

geneList <- TAM_markers$avg_log2FC[TAM_markers$cluster=='APOE.MDM']
names(geneList) = TAM_markers$gene[TAM_markers$cluster=='APOE.MDM']

geneList = sort(geneList, decreasing = TRUE)
library(clusterProfiler)
library(enrichplot)
hallmarks <- read.gmt("~/Brain/brian/Geneset/h.all.v2022.1.Hs.symbols.gmt")
y <- GSEA(geneList, TERM2GENE = hallmarks)
yd <- as.data.frame(y)

dotplot(y,showCategory=30,split=".sign") + facet_grid(~.sign)

gseaplot2(y, "HALLMARK_INTERFERON_GAMMA_RESPONSE",color = "red", pvalue_table = T)
saveRDS(y@result,'./02TAM/DAM_enrichment_GSVA_l3.rds')
saveRDS(y@result,'./02TAM/DAM.Mac_enrichment_GSVA_l3.rds')

pdf('./Figuresubmit/HALLMARK_TNFA_SIGNALING_VIA_NFKB_gsea.pdf',height=8,width=10)
gseaplot2(y, "HALLMARK_TNFA_SIGNALING_VIA_NFKB",color = "red", pvalue_table = T)
dev.off()




##add Primary TAMs and Compare with Brain Figure8B----
seu<-qs::qread('./05addPrimary/seu.anno.qs')
seu.TAM<-subset(seu,celltype=='TAM_DC')
seu.TAM <- NormalizeData(seu.TAM)
seu.TAM <- FindVariableFeatures(seu.TAM, selection.method = "vst", nfeatures = 2000)
seu.TAM<- ScaleData(seu.TAM, features = rownames(seu.TAM))
seu.TAM<- RunPCA(seu.TAM, features = VariableFeatures(object = seu.TAM),reduction.name = "pca")

table(seu.TAM$title)
seu.TAM <- RunHarmony(seu.TAM ,
                      reduction = "pca",
                      group.by.vars = "title",
                      reduction.save = "harmony")
ElbowPlot(seu.TAM, reduction = "harmony", ndims = 30)
xx <- cumsum(seu.TAM[["harmony"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.90) 
ndim = 30
seu.TAM  <- RunUMAP(seu.TAM , 
                    reduction = "harmony", 
                    dims = 1:30, 
                    reduction.name = "umap")
seu.TAM  <- RunUMAP(seu.TAM , 
                    reduction = "pca", 
                    dims = 1:30, 
                    reduction.name = "umap_pca")
seu.TAM  <- RunUMAP(seu.TAM , 
                    reduction = "harmony", 
                    dims = 1:30, 
                    n.neighbors = 80,# 100
                    min.dist = 0.1,  
                    metric = "correlation", 
                    reduction.name = "umap_add")

seu.TAM <- RunTSNE(seu.TAM, reduction = "harmony", 
                   dims = 1:30,
                   reduction.name = "tsne", 
                   perplexity = 150)
DimPlot(seu.TAM, reduction = "umap", split.by = "title")

resolution = seq(0.2,2,0.1)
seu.TAM <- FindNeighbors(seu.TAM, reduction = "harmony", dims = 1:30, k.param = 20)
seu.TAM <- FindClusters(seu.TAM, resolution = seq(0.2,2,0.1))
dir.create('./02TAM/')
pdf('./05addPrimary_TAM/seu.TAM_umap_resolution.pdf',width=18,height=18)
DimPlot(seu.TAM, reduction = "umap", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu.TAM, reduction = "tsne", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
dev.off()

seu.TAM.BM<-qs::qread('./02TAM/seu.TAM.qs')
seu.TAM.BM$celltype.brain<-seu.TAM.BM$celltype
seu.TAM$sample<-rownames(seu.TAM@meta.data)
seu.TAM$celltype %>% table()
seu.TAM@meta.data<-seu.TAM@meta.data %>% 
  left_join(seu.TAM.BM@meta.data[,c('sample','celltype.brain')]) %>% 
  mutate(celltype.brain=ifelse(is.na(celltype.brain),celltype,celltype.brain))
rownames(seu.TAM@meta.data)<-seu.TAM@meta.data$sample
(DimPlot(seu.TAM, reduction = "umap", group.by = 'celltype.brain',label = T))  & ggsci::scale_color_d3("category20")

Idents(seu.TAM)<-seu.TAM$TAM_cluster
DefaultAssay(seu.TAM)<-'RNA'
TAM_markers<-FindAllMarkers(seu.TAM)
TAM_markers <-TAM_markers%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(avg_log2FC))%>%
  ungroup()
d<-c('TAM.0'='APOE.Mac.active',
     'TAM.1'='CD14.Mono',
     'TAM.2'='Alveolar.Mac',
     'TAM.3'='APOE.Mac.active',
     'TAM.4'='APOE.Mac.active', #DAM
     'TAM.5'='APOE.Mac.secretory',#secretory
     'TAM.6'='APOE.Mac',#silent
     'TAM.7'='MT1H.Mac',
     'TAM.8'='APOE.Mac.APC',#APC
     'TAM.9'='DC',
     'TAM.10'='CXCL10.Mac',#CXCL10,IFN
     'TAM.11'='Prolif.Mac',
     'TAM.12'='APOE.Mac',
     'TAM.13'='S100A8.Mac',
     'TAM.14'='S100A8.Mac',
     'TAM.15'='Stressed.Mac',
     'TAM.16'='S100A8.Mac',
     'TAM.17'='DC',
     'TAM.18'='active.DC',
     'TAM.19'='DC',
     'TAM.20'='APOE.Mac',#silent
     'TAM.21'='APOE.Mac')#
seu.TAM$celltype<-d[seu.TAM$TAM_cluster]
seu.TAM$celltype_type<-paste(seu.TAM$Tissue_origins,seu.TAM$celltype,sep='.')
seu.TAM$seurat_clusters <- as.character(seu.TAM$RNA_snn_res.1)
seu.TAM$TAM_cluster<-paste0('TAM.',seu.TAM$seurat_clusters)
Idents(seu.TAM)<-'TAM_cluster'

Idents(seu.TAM)<-seu.TAM$celltype_type
table(seu.TAM$celltype)
DefaultAssay(seu.TAM)<-'RNA'
markers_apoe<-FindMarkers(seu.TAM,
                          ident.1 = c('Brain.APOE.Mac.active'),
                          ident.2 = c('Primary.APOE.Mac.active'))

markers_apoe<-markers_apoe %>% 
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  rownames_to_column('gene')
data.table::fwrite(markers_apoe,
                   './05addPrimary_TAM/markers_APOE.Mac.txt',sep='\t',
                   quote = F,row.names = F)
seu.TAM$resource<-seu.TAM$Tissue_origins

seu$celltype_ty
pdf('./FigureSubmit/primaryvsbrain.pdf',height=5,width=10)
FeaturePlot_scCustom(seu.TAM,features=c('SPP1','VSIG4'),split.by = 'resource')
DimPlot(seu.TAM,group.by = 'celltype')
dev.off()
seu$Tissue_origins<-ifelse(seu$Tissue_origins=='mBrain','Brain','Primary')
seu$resource<-seu$Tissue_origins
FeaturePlot_scCustom(seu,features=c('SPP1'),split.by = 'resource')
DimPlot(seu,group.by = 'celltype',label=T)
qs::qsave(seu.TAM,'./05addPrimary_TAM/seu.TAM_all.qs')
data.table::fwrite(TAM_markers,
                   './05addPrimary_TAM/Primary&Brain_TAMcelltype_markers.txt',sep='\t',
                   quote = F,row.names = F)
data.table::fwrite(TAM_markers2,
                   './05addPrimary_TAM/Primary&Brain_TAMcelltype_cluster_markers_all.txt',sep='\t',
                   quote = F,row.names = F)
seu$celltype_resource<-paste(seu$resource,seu$celltype,sep='.')
table(seu$celltype_resource)
Idents(seu)<-seu$celltype_resource
markers_TAM<-FindMarkers(seu,
                         ident.1 = 'Brain.TAM_DC',
                         ident.2 = 'Primary.TAM_DC')

markers_TAM<-markers_TAM %>% 
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  rownames_to_column('gene')
markers_TAM$p_val_adj[markers_TAM$p_val_adj == 0] <- 1e-300  # 将 0 替换为 1e-10


pdf('./05addPrimary_TAM/ApoEmarkers_brainvsprimary.pdf',height=6,width=6)
data.plot<-subset(markers_apoe,avg_log2FC>=0)
pdf('./05addPrimary_TAM/TAMmarkers_brainvsprimary.pdf',height=6,width=6)
data.plot<-subset(markers_TAM,avg_log2FC>=0)

ggplot(data.plot, 
       aes(x=avg_log2FC, 
           y=-log10(p_val_adj),
           fill=avg_log2FC)) +
  geom_jitter(size = 2, shape = 21, color = 'gray', width = 0.2, height = 0.2) +
  scale_fill_gradient2(low = '#0099CC',
                       mid = "white",
                       high = '#CC3333')+
  ggrepel::geom_text_repel(inherit.aes = F,
                           data = subset(data.plot, 
                                         gene%in%DAM_features),
                           mapping = aes(x=avg_log2FC,
                                         y=-log10(p_val_adj),
                                         label = gene),
                           nudge_x = .1, nudge_y = .01,
                           max.overlaps = Inf,
                           color='black',
                           size = 6) +
  geom_jitter(data = subset(data.plot, 
                            gene%in%DAM_features),
              aes(x=avg_log2FC,
                  y=-log10(p_val_adj),
                  fill=avg_log2FC),
              size=4,shape=21,color='black')+
  # scale_x_continuous(limits=c(-4,6)) +
  # scale_y_continuous(limits=c(0,0.2)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        
        legend.position = 'right',
        legend.text = element_text(size=12),
        panel.grid = element_blank(),
        axis.title =element_text(size=14),
        axis.text =element_text(color = "black",size=14))+
  # xlab('log2FoldChange')+
  # ylab('log10(P.val)')+
  labs(color = '', 
       x = expression("log"[2]*"FoldChange"), 
       y = expression("log"[10]*"P.Val"))

dev.off()
