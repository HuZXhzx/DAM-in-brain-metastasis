library(GEOquery)
library(tidyverse)
library(arrayQualityMetrics)
library(affy)
library(limma)
library(Seurat)
library(harmony)
library(Nebulosa)
library(scCustomize)
##seu.TAM regular-----
seu<-qs::qread('./dataset6/Lung/01merge/seu.anno.qs')
seu.TAM<-subset(seu,celltype%in%c('TAM'))
DimPlot(seu.TAM,split.by = 'geo_accession')
seu.TAM <- NormalizeData(seu.TAM)
seu.TAM <- FindVariableFeatures(seu.TAM, selection.method = "vst", nfeatures = 2000)
seu.TAM<- ScaleData(seu.TAM, features = rownames(seu.TAM))
seu.TAM<- RunPCA(seu.TAM, features = VariableFeatures(object = seu.TAM),reduction.name = "pca")
table(seu.TAM$geo_accession)
seu.TAM <- RunHarmony(seu.TAM ,
                      reduction = "pca",
                      group.by.vars = "geo_accession",
                      reduction.save = "harmony")
ElbowPlot(seu.TAM, reduction = "harmony", ndims = 30)
xx <- cumsum(seu.TAM[["harmony"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.80) 
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
DimPlot(seu.TAM, reduction = "umap", group.by = "geo_accession")+
  DimPlot(seu.TAM, reduction = "umap_pca", group.by = "geo_accession")+
  DimPlot(seu.TAM, reduction = "umap_add", group.by = "geo_accession")
DimPlot(seu.TAM, reduction = "umap_add", split.by = "geo_accession")
table(seu.TAM$age,seu.TAM$geo_accession)
resolution = seq(0.2,2,0.1)
DefaultAssay(seu.TAM)<-'RNA'
seu.TAM <- FindNeighbors(seu.TAM, reduction = "harmony", dims = 1:30, k.param = 20)
seu.TAM <- FindClusters(seu.TAM, resolution = seq(0.2,2,0.1))
dir.create('./dataset6/Lung/02TAM/')
pdf('./dataset6/Lung/02TAM/seu.TAM_umap_resolution.pdf',width=18,height=18)
DimPlot(seu.TAM, reduction = "umap", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu.TAM, reduction = "umap_add", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu.TAM, reduction = "tsne", group.by = paste0("RNA_snn_res.", resolution), ncol = 3, label = T)
dev.off()
seu.TAM$seurat_clusters <- as.character(seu.TAM$RNA_snn_res.0.7)
Idents(seu.TAM)<-seu.TAM$seurat_clusters


seu.TAM$TAM_cluster<-paste0('TAM_',as.character(seu.TAM$RNA_snn_res.0.6))
Idents(seu.TAM)<-seu.TAM$TAM_cluster
TAM_markers<-FindAllMarkers(seu.TAM)
TAM_markers <-TAM_markers%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(avg_log2FC))%>%
  ungroup()

TAM_markers$cluster<-paste('TAM',as.character(TAM_markers$cluster),sep='.')
seu.TAM$TAM_cluster<-paste('TAM',as.character(seu.TAM$seurat_clusters),sep='.')
seu.TAM$sample<-rownames(seu.TAM@meta.data)
seu$TAM_cluster<-NULL
seu@meta.data<-seu@meta.data %>% 
  left_join(seu.TAM@meta.data[,c('sample','TAM_cluster')]) %>% 
  mutate(TAM_cluster=ifelse(is.na(TAM_cluster),
                            'others',TAM_cluster))
table(seu$TAM_cluster)
rownames(seu@meta.data)<-seu@meta.data$sample
Idents(seu)<-seu$TAM_cluster
Idents(seu)<-factor(Idents(seu),levels = c(paste0('TAM.',0:10),'others'))
DefaultAssay(seu)<-'RNA'
TAM_markers2<-FindAllMarkers(seu)
TAM_markers2 <-TAM_markers2%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(avg_log2FC))%>%
  ungroup() %>% 
  filter(!(cluster=='others'))
d <- c(
  "TAM.0"="APOE.MDM",
  "TAM.1"="APOE.MDM.DAM",
  "TAM.2" = "recycling.MDM", #HSPA+
  "TAM.3" = "DC",
  "TAM.4" = "Lung derived MDM", 
  "TAM.5"="DAM.Microglia",
  "TAM.6"="APOE.MDM",
  "TAM.7"="Mac_MT1H", #stressed HSPA APOE.MDM
  "TAM.8" = "APOE.MDM.DAM", 
  "TAM.9" = "DC",
  "TAM.10" = "APOE.MDM"
)
seu.TAM$celltype.new<-d[seu.TAM$TAM_cluster]

##Findmarkers TAM celltype-----
DefaultAssay(seu.TAM)<-'RNA'
Idents(seu.TAM)<-seu.TAM$celltype.new
TAM_markers<-FindAllMarkers(seu.TAM)
TAM_markers <-TAM_markers%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(avg_log2FC))%>%
  ungroup()
seu<-qs::qread('./Lung/03addTNBC/seu.all.lung.breast.qs')

seu@meta.data<-seu@meta.data %>% 
  left_join(seu.TAM@meta.data[,c('sample','celltype.new')]) %>% 
  mutate(celltype.new=ifelse(is.na(celltype.new),
                             'others',celltype.new))
table(seu$celltype.new)
rownames(seu@meta.data)<-seu@meta.data$sample
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
                   './Lung/04clusterMDM/TAMcelltype_markers.txt',sep='\t',
                   quote = F,row.names = F)
data.table::fwrite(TAM_markers2,
                   './Lung/04clusterMDM/TAMcelltype_cluster_markers_all.txt',sep='\t',
                   quote = F,row.names = F)


##Dimplot Figure3A-----

celltype_pl<-d<-c('DAM.Microglia'='DAM',
                  'APOE.MDM'='APOE.Mac',
                  'APOE.MDM.DAM'='DAM-like.Mac',
                  'DC'='DC',
                  'recycling.MDM'='Prolif.Mac',
                  'Mac_MT1H'='MT1H.Mac',
                  'Lung derived MDM'='Lung.derived.Mac')
table(seu.TAM$celltype_pl)
seu.TAM$celltype_pl<-celltype_pl[seu.TAM$celltype.new]
color.TAM_hgc<-readRDS('../DAM/color.TAM_hgc.rds')

pdf('./Figuresubmit/TAM_Dimplot_pl.pdf',width=5.8,height=4.4)
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
saveRDS(seu.TAM,'./Figuresubmit/seu.TAM.rds')


# Ucell BMDM [Bowman et.al] Figure3B------
MG_dataset<-qs::qread('./Geneset/MG_dataset_forhuman.qs')
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


marker_genes <- colnames(seu.TAM@meta.data)[grep("_UCell",colnames(seu.TAM@meta.data))]
data<-seu.TAM@meta.data[,marker_genes]
colnames(data)<-names(MG_BMDM)
data<-t(data)
colnames(data)<-colnames(seu.TAM)
seu.TAM[['MG_BMDM']]<-Seurat::CreateAssayObject(data = data)

DefaultAssay(seu.TAM)<-'MG_BMDM'

pdf('./Figuresubmit/MG_BMDM_avg.Recluster.pdf',height=2.5,width=6)
Idents(seu.TAM)<-factor(seu.TAM$celltype_pl,levels=c('DAM','APOE.Mac','DAM-like.Mac',
                                                     'Prolif.Mac',
                                                     'MT1H.Mac','Lung.derived.Mac','DC'))
DefaultAssay(seu.TAM)<-'MG_BMDM'

AverageHeatmap(object = seu.TAM ,
               assays = 'MG_BMDM',
               markerGene =c('Core-MG','Core-BMDM','TAM-MG','TAM-BMDM'),
               clusterAnnoName = F,height=2,width=6,
               myanCol = color.TAM_hgc[levels(Idents(seu.TAM))],
               annoCol = T)
dev.off()
# Ucell DAM [Gan et.al]----
library(UCell)
MG_dataset<-data.table::fread('./Geneset/Gan_et_al_MGdataet1_hg.txt')
MG_list<-split(MG_dataset$Gene, MG_dataset$Category)
seu.TAM<- AddModuleScore_UCell(seu.TAM, 
                               features = MG_list,
                               assay = 'RNA',
                               slot='counts',
                               name = "_UCell")
marker_genes <- colnames(seu.TAM@meta.data)[grep("_UCell",colnames(seu.TAM@meta.data))][5:11]
data<-seu.TAM@meta.data[,marker_genes]
colnames(data)<-names(MG_list)
data<-t(data)
colnames(data)<-colnames(seu.TAM)
seu.TAM[['DAM']]<-Seurat::CreateAssayObject(data = data)
# Ucell TAM [Katrina T et.al]----
MG_dataset<-qs::qread('./Geneset/MG_dataset_forhuman.qs')

TAM.features<-MG_dataset %>% 
  filter(category%in%c('APC_Topic_Score',
                       'IFN_Response_Topic_Score',
                       'Secretory_Topic_Score'))

TAM.features<-split(TAM.features$gene, TAM.features$category)
DefaultAssay(seu.TAM)<-'RNA'
seu.TAM<- AddModuleScore_UCell(seu.TAM, 
                               features = TAM.features,
                               name = "_UCell")

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

saveRDS(seu.TAM,'./Lung/04clusterMDM/seu.TAM.rds')

##DAM&TAM features Figure4A------
pdf('./Figuresubmit/TAM.features.pdf',height=2.5,width=6)
Idents(seu.TAM)<-factor(seu.TAM$celltype_pl,levels=c('DAM','APOE.Mac','DAM-like.Mac',
                                                     'Prolif.Mac',
                                                     'MT1H.Mac','Lung.derived.Mac','DC'))

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
                     levels = rev(c('DAM','APOE.Mac','DAM-like.Mac',
                                    'Prolif.Mac',
                                    'MT1H.Mac','Lung.derived.Mac','DC')))
data.plot.wide<-data.plot[,c(1,4,5)] %>% 
  pivot_wider(names_from = 'id',values_from = 'avg.exp.scaled') %>% 
  column_to_rownames('gene')
saveRDS(data.plot.wide,'./Figuresubmit/DAM.wide.dataset6.rds')
Idents(seu.TAM)<-seu.TAM$celltype_pl

AverageExp<-AverageExpression(seu.TAM,
                              assays = 'DAM',
                              features=c("DAM-core-program", 
                                         "Stage1-DAM-upregulated",
                                         "Stage2-DAM" ))$DAM

p <- t(AverageExp) %>% scale() %>% t()

saveRDS(p,'./Figuresubmit/DAM.wide.dataset6_v2.rds')

AverageExp<-AverageExpression(seu.TAM,
                              assays = 'TAM.features',
                              features=c(c('APC-Topic-Score','IFN-Response-Topic-Score',
                                           'Secretory-Topic-Score')))$TAM.features

p <- t(AverageExp) %>% scale() %>% t()
saveRDS(p,'./Figuresubmit/TAM.features.wide.dataset6.rds')
##TAM featuresplot  Figure3C------
DefaultAssay(seu.TAM)<-'RNA'
TAM_features<-c( "TMEM119","P2RY12", "CX3CR1", 
                 "CSF1R","OLFML3","GPR34",
                 "ITGA4","PLAC8",#'CCR2','CD14',
                 'S100A8',"CXCL8",'FCN1',
                 "APOE","C1QB","SPP1",
                 "MRC1","CD163","SIGLEC1",
                 #'CAMP', 'S100A9', 
                 'STMN1','TOP2A',
                 #'HSPA1B', 'HSPB1',
                 'MT1H','MT1M',
                 'SCGB3A1','SFTPA2',
                 'CST3','CD1C',
                 "HLA-DPB1",'CD74')
TAM_features<-data.frame(gene=TAM_features,
                         cluster=c(rep('Homeo.MG',6),
                                   rep('BMDM',2),
                                   rep('Mono&S100A8.Mac',3),
                                   rep('APOE.Mac',3),
                                   rep('Mf',3),
                                   rep('Prolif.MG',2),
                                   rep('MT1H.Mac',2),
                                   rep('Lung.Mac',2),
                                   rep('DC',2),
                                   rep('AP',2)))
DefaultAssay(seu.TAM)<-'RNA'
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

gene<-setdiff(TAM_features$gene,unique(data.plot$gene))
celltype<-unique(data.plot$id)
colnames(data.plot)
ma_impue<-data.frame(id=rep(celltype,length(gene)),
                     avg.exp=rep(0,length(rep(celltype,length(gene)))),
                     pct.exp=rep(0,length(rep(celltype,length(gene)))),
                     gene=rep(gene,each=length(celltype)),
                     avg.exp.scaled=rep(0,length(rep(celltype,length(gene)))))

data.plot<-rbind(data.plot,ma_impue) %>% as.data.frame()
data.plot$gene <- factor(data.plot$gene, levels = TAM_features$gene)

data.plot$id<-factor(data.plot$id, 
                     levels = rev(c('DAM','APOE.Mac','DAM-like.Mac',
                                    'Prolif.Mac',
                                    'MT1H.Mac','Lung.derived.Mac','DC')))

dir.create('./Figuresubmit/')
pdf('./Figuresubmit/MG_features.pdf',height=5,width=10)
ggplot(data.plot,aes(x = gene, y = id)) +
  geom_point(aes(fill = avg.exp.scaled, size =avg.exp.scaled),
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

##DAM genes Hotmap Figure4E----

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
saveRDS(data.plot.wide,'./Figuresubmit/data.plot.wide.dataset6.rds')

##DAM cellmarkersFigure4C,D ------
TAM_markers<-data.table::fread('./Lung/04clusterMDM/TAMcelltype_cluster_markers_all.txt')
DAM_markers<-readRDS('../DAM/KEY_Gene_DAM_ALL.rds')
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
TAM_markers<-data.table::fread('./Lung/06trajectory.new/TAMcelltypetra_markers_all.txt')
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
##HALLMARK_TNFA_SIGNALING_VIA_NFKB----
TAM_markers<-data.table::fread('./Lung/04clusterMDM/TAMcelltype_cluster_markers_all.txt')
table(TAM_markers$cluster)
geneList <- TAM_markers$avg_log2FC[TAM_markers$cluster=='DAM.Microglia']
names(geneList) = TAM_markers$gene[TAM_markers$cluster=='DAM.Microglia']

geneList = sort(geneList, decreasing = TRUE)
library(clusterProfiler)
library(enrichplot)
hallmarks <- read.gmt("../Geneset/h.all.v2022.1.Hs.symbols.gmt")
y <- GSEA(geneList, TERM2GENE = hallmarks)
yd <- as.data.frame(y)
dotplot(y,showCategory=30,split=".sign") + facet_grid(~.sign)

pdf('./Figuresubmit/HALLMARK_TNFA_SIGNALING_VIA_NFKB_gsea.pdf',height=8,width=10)
gseaplot2(y, "HALLMARK_TNFA_SIGNALING_VIA_NFKB",color = "red", pvalue_table = T)
dev.off()

