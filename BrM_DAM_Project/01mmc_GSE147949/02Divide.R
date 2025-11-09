## seuTAM Baseline------
rm(list=ls())
library(Seurat);library(tidyverse);library(Nebulosa);library(clustree);library(GEOquery)
seu<-qs::qread('dataset3/tmp/03-1.mmc.seurat.aucell.qs')
class(seu$seurat_clusters)

seu.TAM<-subset(seu,(seurat_clusters %in% c("15","18"))|(celltype_cluster=='microglia'))
DimPlot(seu.TAM,label=T,reduction='umap',group.by='seurat_clusters')
DimPlot(seu,label=T,reduction='umap',group.by='seurat_clusters')

DefaultAssay(seu.TAM)<-'RNA'
seu.list<-SplitObject(seu.TAM, split.by = "geo_accession")
for (i in 1:length(seu.list)) {
  seu.list[[i]] <- NormalizeData(seu.list[[i]], verbose = FALSE)
  seu.list[[i]] <- FindVariableFeatures(seu.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}
names(seu.list)
seu.anchors <- FindIntegrationAnchors(object.list = seu.list, dims = 1:20,
                                      anchor.features = 2000)

seu.TAM.integrated <- IntegrateData(anchorset = seu.anchors, dims = 1:20)
DefaultAssay(seu.TAM.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
seu.TAM.integrated <- ScaleData(seu.TAM.integrated, verbose = FALSE)
seu.TAM.integrated <- RunPCA(seu.TAM.integrated, npcs = 80, verbose = FALSE)
DimHeatmap(seu.TAM.integrated, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(seu.TAM.integrated, dims = 16:30, cells = 500, balanced = TRUE)

ElbowPlot(seu.TAM.integrated, reduction = "pca", ndims = 30)
xx <- cumsum(seu.TAM.integrated[["pca"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.9) # 20 PCs解释了95%的方差，假设5%的方差来自于噪声
ndim = 20

## try a larger perplexity!
seu.TAM.integrated <- RunTSNE(seu.TAM.integrated, reduction = "pca", dims = 1:20, perplexity = 150)
seu.TAM.integrated <- RunUMAP(seu.TAM.integrated, reduction = "pca", dims = 1:20)
table(seu.TAM.integrated$type)

DimPlot(seu.TAM.integrated, reduction = "umap", group.by = "geo_accession")

resolution = seq(0.2,1.2,0.1)
seu.TAM.integrated <- FindNeighbors(seu.TAM.integrated, reduction = "pca", dims = 1:20, k.param = 20)
seu.TAM.integrated <- FindClusters(seu.TAM.integrated, resolution = seq(0.2,1.2,0.1))
dir.create('./dataset3/TAM/')
pdf('./dataset3/TAM/seu_TAM_umap_resolution2.pdf',width=12,height=12)
DimPlot(seu.TAM.integrated, reduction = "umap", group.by = paste0("integrated_snn_res.", resolution), ncol = 3, label = T)
DimPlot(seu.TAM.integrated, reduction = "tsne", group.by = paste0("integrated_snn_res.", resolution), ncol = 3, label = T)
dev.off()

seu.TAM.integrated@meta.data$seurat_clusters <- seu.TAM.integrated@meta.data$integrated_snn_res.0.5
class(seu.TAM.integrated@meta.data$seurat_clusters)
seu.TAM.integrated@meta.data$seurat_clusters<-as.character(seu.TAM.integrated@meta.data$seurat_clusters)
seu.TAM.integrated@meta.data$integrated_snn_res.1.2<-as.character(seu.TAM.integrated@meta.data$integrated_snn_res.1.2)

seu.TAM.integrated@meta.data$seurat_clusters[seu.TAM.integrated@meta.data$integrated_snn_res.1.2==18]<-'15' #'Apoe+ BMDM'
table(seu.TAM.integrated@meta.data$seurat_clusters)
Idents(seu.TAM.integrated)<-'seurat_clusters'
seu.TAM.integrated@meta.data$seurat_clusters<-as.factor(seu.TAM.integrated@meta.data$seurat_clusters)
seu.TAM$seurat_clusters[as.character(seu.TAM$integrated_snn_res.0.7)=='13' & seu.TAM$seurat_clusters=='6'] <- '6_2'
seu.TAM$seurat_clusters[as.character(seu.TAM$integrated_snn_res.0.7)=='9' & seu.TAM$seurat_clusters=='6' ] <-'6_1'
seu.TAM$seurat_clusters[seu.TAM$seurat_clusters=='6' ] <-'6_1'
seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('0','1','2','4','5')] <-'c0'
seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('3')] <-'c1'
seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('7')] <-'c2'
seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('8')] <-'c3'
seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('6_1')] <-'c4'
seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('6_2')] <-'c5'
seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('9')] <-'c6'
seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('11')] <-'c7'
seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('14')] <-'c8'
seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('10')] <-'c9'
seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('15')] <-'c10'

seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('12')] <-'c12'
seu.TAM@meta.data$TAM_cluster[seu.TAM$seurat_clusters%in%c('13')] <-'c11'
seu.TAM@meta.data$TAM_cluster[seu.TAM@meta.data$TAM_cluster=='c11'] <-'c9'
seu.TAM@meta.data$TAM_cluster[seu.TAM@meta.data$TAM_cluster=='c12'] <-'c11'
Idents(seu.TAM)<-seu.TAM$TAM_cluster
DimPlot(seu.TAM,reduction = 'umap',label=T)



seu.TAM$celltype_cluster %>% table()
seu.TAM$celltype_cluster[seu.TAM$seurat_clusters=='15']<-'MDM'
seu.TAM$celltype_cluster[seu.TAM$seurat_clusters=='10']<-'Monocytes'
seu.TAM$celltype_cluster[seu.TAM$TAM_cluster=='c9']<-'Monocytes'

table(subset(seu.TAM,celltype_cluster=='myeloid_cells_DC ')$seurat_clusters)
seu.TAM$celltype_cluster[seu.TAM$seurat_clusters%in%c('0','1','2','4','5','3','7','8','6_1','6_2','11','9','14')]<-'Microglia'

Idents(seu.TAM)<-seu.TAM$celltype_cluster

saveRDS(color.TAM,'../DAM/color.TAM_hcluster.mgi.rds')
color.TAM<-c( 'MDM'="#B49D99", 
              'Monocytes'="#D8AC60",
              'matureDC'= "#C48244", 
              'Microglia'="#FB9A99")
color.TAM_cluster<-c( c0="#969D62",                               
                      c3="#815AA8", 
                      c1="#E93A3B", 
                      c6="#B294C7", 
                      c10="#B49D99", 
                      c2="#B15928",
                      c5="#FDAE53",        
                      c7="#92CF72",        
                      c4="#84B8D7",
                      c9="#D8AC60",    
                      c8="#F47A79",
                      c11="#C48244" )


Idents(seu.TAM)<-factor(Idents(seu.TAM),
                        levels=paste0('c',0:11))

pdf('./Figuresubmit/TAM_Dimplot.pdf',width=8,height=8)
clusterCornerAxes(object = seu.TAM,reduction = 'umap',
                  noSplit = T,
                  #groupFacet = 'type',
                  clusterCol = 'TAM_cluster',
                  cellLabel = T,
                  keySize = 8,aspect.ratio = 1,
                  relLength = 0.5,themebg = 'bwCorner')& 
  scale_color_manual(values = color.TAM_cluster) 
dev.off()


DimPlot(seu.TAM,group.by = 'TAM_cluster',split.by = 'type',label=T)
pdf('./Figuresubmit/TAM_Dimplot2.pdf',width=8,height=8)
clusterCornerAxes(object = seu.TAM,
                  reduction = 'umap',
                  addCircle = T,
                  clusterCol = 'celltype',
                  noSplit = T,cellLabel = T,
                  keySize = 8,aspect.ratio = 1,
                  relLength = 0.5,themebg = 'bwCorner')& 
  scale_color_manual(values = c(`Br-resident MG`="#FB9A99",BMDM="#A6CEE3") )&
  scale_fill_manual(values = c(`Br-resident MG`="#FB9A99",BMDM="#A6CEE3") )
dev.off()


## MG BMDM features ["Core_BMDM","Core_MG","TAM_BMDM","TAM_MG"]-----

MG_dataset<-qs::qread('./Geneset/MG_dataset_formu.qs')
MG_BMDM<-MG_dataset %>% filter(category%in%c("Core_BMDM","Core_MG","TAM_BMDM","TAM_MG"))
MG_BMDM<-split(MG_BMDM$gene, MG_BMDM$category)
library(UCell)
DefaultAssay(seu.TAM)<-'RNA'
seu.TAM<- AddModuleScore_UCell(seu.TAM, 
                               features = MG_BMDM,name = "_UCell")

marker_genes <- colnames(seu.TAM@meta.data)[grep("_UCell",colnames(seu.TAM@meta.data))]
data<-seu.TAM@meta.data[,marker_genes]
colnames(data)<-names(MG_BMDM)
data<-t(data)
colnames(data)<-colnames(seu.TAM)
seu.TAM[['MG_BMDM']]<-Seurat::CreateAssayObject(data = data)

seu.TAM@assays$MG_BMDM@data %>% View()

seu.TAM@meta.data[,marker_genes] %>% head()
seu.TAM@meta.data[,marker_genes]<-NULL
DefaultAssay(seu.TAM)<-'MG_BMDM'

FeaturePlot(seu.TAM,features=c('Core-MG','Core-BMDM'),cols=col)
FeaturePlot(seu.TAM,features=c('TAM-MG','TAM-BMDM'))
library(scCustomize)
library(RColorBrewer)
library(scRNAtoolVis)
col <- colorRampPalette(brewer.pal(11, "PiYG"))(25)
#col<-c("#8E0152","#DF7CB1","#FAD5E9","#C7E79E","#6EAE36","#276419")
colorRampPalette(c('#0099CC',"white",'#CC3333'))(40)
DefaultAssay(seu.TAM)<-'MG_BMDM'
pdf('./Figuresubmit/MGBMDM2.pdf',width=6,height=5)
FeaturePlot_scCustom(seurat_object = seu.TAM, 
                     raster.dpi = c(50,50),
                     features = c('Core-MG','Core-BMDM',
                                  'TAM-MG','TAM-BMDM'), 
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),)
dev.off()

## MG features ['APC_Topic_Score','IFN_Response_Topic_Score','Secretory_Topic_Score';DAM related]-----
library(UCell)
MG_dataset<-qs::qread('./Geneset/MG_dataset_formu.qs')
TAM.features<-MG_dataset %>% 
  filter(category%in%c('APC_Topic_Score',
                       'IFN_Response_Topic_Score',
                       'Secretory_Topic_Score'))
TAM.features<-split(TAM.features$gene, TAM.features$category)
seu.TAM<- AddModuleScore_UCell(seu.TAM, features = TAM.features,name = "_UCell")
MG_dataset<-data.table::fread("./Geneset/Gan_et_al_MGdataet1_mmc.txt")
MG_list<-split(MG_dataset$Gene, MG_dataset$Category)
DefaultAssay(seu.TAM)<-'RNA'
seu.TAM<- AddModuleScore_UCell(seu.TAM, 
                               features = MG_list,
                               name = "_UCell")

seu.MG<-subset(seu.TAM,TAM_cluster%in%paste0('c',0:8))
seu.MG$TAMcluster_type<-paste0(seu.MG$TAM_cluster,'(',seu.MG$type,')')
Idents(seu.MG)<-factor(seu.MG$TAMcluster_type,
                       levels=c('c0(Con)','c0(Met)',
                                'c1(Con)','c1(Met)',
                                'c2(Con)','c2(Met)',
                                'c3(Con)','c3(Met)',
                                'c4(Con)','c4(Met)',
                                'c5(Con)','c5(Met)',
                                'c6(Con)','c6(Met)',
                                'c7(Con)','c7(Met)',
                                'c8(Con)','c8(Met)'))
### Visualization=====
table(seu.MG$TAMcluster_type)
color.mg.type<-rep(color.TAM_cluster[paste0('c',0:8)],each=2)
names(color.mg.type)<-c('c0(Con)','c0(Met)',
                        'c1(Con)','c1(Met)',
                        'c2(Con)','c2(Met)',
                        'c3(Con)','c3(Met)',
                        'c4(Con)','c4(Met)',
                        'c5(Con)','c5(Met)',
                        'c6(Con)','c6(Met)',
                        'c7(Con)','c7(Met)',
                        'c8(Con)','c8(Met)')
pdf('./Figuresubmit/TAM.features.averageplot.pdf',height = 2,width=10)
AverageHeatmap(object = seu.MG,
               assays = "TAM.features",
               markerGene = c('APC-Topic-Score','IFN-Response-Topic-Score',
                              'Secretory-Topic-Score'),
               myanCol = color.mg.type,
               annoCol =T,
               clusterAnnoName = T,height=2,width=12)

dev.off()
pdf('./Figuresubmit/DAM.averageplot.pdf',height = 2,width=10)
AverageHeatmap(object = seu.MG,
               assays = "DAM",height=2,width=12,
               markerGene = c("DAM-core-program", 
                              "Stage1-DAM-upregulated",
                              "Stage2-DAM",
                              "Homeostasis"),
               myanCol = color.mg.type,
               annoCol =T,
               clusterAnnoName = T)

dev.off()


## MG feature markers------
TAM_features<-c("Tmem119", "P2ry12", 
                "Cx3cr1","Csf1r",
                "Egr1","Jun",
                "Rps20","Rps24",
                "Cd74", "H2-Aa",
                "Isg15","Ifit3",
                "Spp1","Ccl4",
                "Top2a","Stmn1",
                "Pf4","Cd38")

TAM_features<-data.frame(gene=TAM_features,
                         cluster=c(rep('M0.MG',4),
                                   rep('Diff.MG',2),
                                   rep('Trans.MG',2),
                                   rep('AP.MG',2),
                                   rep('Def.MG',2),
                                   rep('Sec.MG',2),
                                   rep('Prolif.MG',2),
                                   rep('BAM',2)))
geneExp <-FetchData(
  object = subset(seu.TAM,TAM_cluster%in%paste0('c',0:8)),
  vars = TAM_features$gene,
  slot = 'data'
)
geneExp$TAM_cluster <- subset(seu.TAM,TAM_cluster%in%paste0('c',0:8))@meta.data[['TAM_cluster']]
geneExp$type<- subset(seu.TAM,TAM_cluster%in%paste0('c',0:8))@meta.data[['type']]
geneExp$id<-paste0(geneExp$TAM_cluster,'(',geneExp$type,')')
geneExp$id<-geneExp$TAM_cluster
PercentAbove <- utils::getFromNamespace("PercentAbove", "Seurat")
data.plot <- lapply(
  X = unique(geneExp$id),
  FUN = function(ident) {
    data.use <- geneExp[geneExp$id == ident, 1:(ncol(geneExp) - 3), drop = FALSE]
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

# no scale or rescale using log1p
#data.plot$avg.exp.scaled <- log1p(data.plot$avg.exp)
data.plot.res <- data.plot
data.plot.res$gene <- factor(data.plot.res$gene, levels = unique(data.plot.res$gene))


data.plot.res$id<-factor(data.plot.res$id, levels = rev(paste0('c',0:8)))

pdf('./Figuresubmit/MG_features.pdf',height=4,width=8)
ggplot(data.plot.res,aes(x = gene, y = id)) +
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
  guides(size = guide_legend(title = "Fraction of cells \n in group (%)",
                             title.position = "top",
                             title.hjust = 0.5,
                             label.position = "bottom",
                             override.aes = list(color = "black",fill = "grey50"),
                             keywidth = unit(0.3, "cm")))+
  scale_fill_gradient2(low = '#0099CC',mid = "white",high = '#CC3333')
dev.off()
pdf('./Figuresubmit/MG_features_add.pdf',height=4,width=8)

jjDotPlot(object = subset(seu.TAM,TAM_cluster%in%paste0('c',0:8)),
          markerGene  = TAM_features, 
          anno = T,
          plot.margin = c(1,1,1,1), lwd = 0.5,textSize = 14,
          dot.col = c('#0099CC',"white",'#CC3333'),
          id = 'TAM_cluster',ytree = F)
dev.off()
##cell abundance-====
df<-data.frame(table(seu.TAM$TAM_cluster,seu.TAM$type)) 
df<- df %>% group_by(Var2) %>% 
  mutate(percent=Freq/sum(Freq)) %>% 
  rename(Celltype=Var1,
         Group=Var2) %>% 
  mutate(Celltype=factor(Celltype,levels=paste0('c',0:11)),
         Group=factor(Group,levels=c('Con','Met'))) %>% 
  arrange(Group,Celltype) 


df$Subject<-rep(1:12,2)
library(ggalluvial)
pdf('./Figuresubmit/Cellabundance_alluvial.pdf',height=6,width=4)
ggplot(df, aes(x = Group, stratum = Celltype, alluvium = Subject, y = percent , 
               fill = Celltype,label = Celltype, position = "fill")) +
  scale_x_discrete(expand = c(.1, .2)) +
  #scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  geom_flow(alpha = 0.5, width = 1/2,) +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        axis.text.y = element_text(lineheight=.8,  
                                   hjust=0.5, size =14),
        axis.text.x = element_text(size =12),
        axis.title.x = element_text(lineheight=.8,
                                    hjust=0.5, size =14),
        legend.title = element_blank())+
  scale_fill_manual(values = color.TAM_cluster)+
  geom_stratum(alpha = 1, width = 1/2) +
  #geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "right") + ylab("Ratio") + xlab("") 
dev.off()

## TAM markers-----
table(seu.TAM$TAMcluster_type)
Idents(seu.TAM) <-seu.TAM$TAMcluster_type
DefaultAssay(seu.TAM)<-'RNA'
TAM_markers.type <- FindAllMarkers(seu.TAM, only.pos = TRUE)

TAM_markers.type <-TAM_markers.type %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(abs(avg_log2FC)))%>% 
  ungroup()
colnames(TAM_markers.type)[6]<-'TAM_cluster_type'

data.table::fwrite(TAM_markers.type,
                   './V2/02markers_TAM_clusters_type.txt',sep='\t',
                   quote = F,row.names = F)

## DAM(Met.c6) Met vs Con-----
data<-TAM_markers.type %>% 
  filter(TAM_cluster_type%in%c('Met.c6')) %>% 
  arrange(TAM_cluster_type,desc(avg_log2FC)) %>% 
  filter(avg_log2FC>=0)
data<-data%>% 
  mutate(Gene=1:nrow(data))
data <- head(data, n=500)
DAM_gene<-readRDS('../DAM/KEY_Gene_DAM_ALL_MGI.rds') %>% 
  distinct(Gene.name,.keep_all = T)
de.gene<-data.table::fread('./V2/markers_convsmet.txt') %>% filter(cluster=='c6')
DAM_gene_up<-DAM_gene$MGI.symbol[DAM_gene$MGI.symbol%in%de.gene$gene[de.gene$avg_log2FC>=0]]
DAM_gene_down<-DAM_gene$MGI.symbol[DAM_gene$MGI.symbol%in%de.gene$gene[de.gene$avg_log2FC<-0.5]]


data$regulate[data$gene%in%DAM_gene$MGI.symbol]<-'DAM key gene'
data$regulate[data$gene%in%DAM_gene_up]<-'Up in brain metastasis'
data$regulate[data$gene%in%DAM_gene_down]<-'Down in brain metastasis'
data$regulate[!(data$gene%in%DAM_gene$MGI.symbol)]<-'NoSig'
data.label<-subset(data,gene%in%DAM_gene$MGI.symbol)
data.label.up<-subset(data.label,gene%in%DAM_gene_up)
data.label.down<-subset(data.label,gene%in%DAM_gene_down)

d<-c("Up in brain metastasis"='#CC3333',
     "Down in brain metastasis"='#0099CC',
     "NoSig"="#BECEE3",
     "DAM key gene"="#B294C7")
#data.label<-data[data$gene%in%DAM_gene$MGI.symbol]
#data$pt.col <- ifelse(data$gene %in%data.label$gene, "#B294C7", "#BECEE3")
pdf('./Figuresubmit/DAMgene.pdf',height=7,width=8)
ggplot(data, aes(Gene, avg_log2FC)) +
  geom_point(size=3,aes(color=regulate)) +
  ggrepel::geom_text_repel(inherit.aes = FALSE, max.overlaps = Inf,
                           data = data.label, 
                           aes(Gene, avg_log2FC,label=gene,
                               color=regulate), size=6, show.legend = FALSE) +
  scale_color_manual(values=c(`Up in brain metastasis`='#CC3333',
                              `Down in brain metastasis`='#0099CC',
                              NoSig="#BECEE3",
                              `DAM key gene`="#B294C7"))+
  #ggtitle("Metc6") + 
  ylab("Specificity score") +
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


## markers con vs met------

table(seu.TAM$TAM_cluster)
seu.TAM$TAMcluster_type<-paste0(seu.TAM$type,seu.TAM$TAM_cluster)
Idents(seu.TAM) <-'TAMcluster_type'
celltypes <- unique(seu.TAM$TAM_cluster)
groups <- unique(seu.TAM$type)


DefaultAssay(seu.TAM)<-'RNA'
de.list <- lapply(celltypes, function(ct) {
  message(glue::glue("processing {ct} ..."))
  ct1 <- paste0(groups[1],ct)
  ct2 <- paste0(groups[2],ct)
  n1<-length(colnames(subset(seu.TAM,TAMcluster_type==ct1)))
  n2<-length(colnames(subset(seu.TAM,TAMcluster_type==ct2)))
  if(n1>10&n2>10){
    de <- FindMarkers(seu.TAM, ident.1 = ct2, 
                      ident.2 = ct1, 
                      test.use = "wilcox", 
                      fc.name = "avg_diff", 
                      logfc.threshold = 0)
    de$change <- ifelse(de$avg_diff > 0,
                        paste("up in", groups[2]),
                        paste("up in", groups[1]))
    de$cluster <- ct
    return(de)
  }
})

names(de.list) <- celltypes

filtered_list <- de.list[!sapply(de.list, is.null)]
filtered_list<-lapply(seq_along(filtered_list),function(xx){
  filtered_list[[xx]]<-filtered_list[[xx]] %>% 
    filter(p_val<0.05) %>% 
    rownames_to_column('gene')
  return(filtered_list[[xx]])
})
de.gene <- purrr::reduce(filtered_list,rbind) 
de.gene <-de.gene %>% arrange(cluster,desc(avg_diff)) 
colnames(de.gene)[3]<-'avg_log2FC'
data.table::fwrite(de.gene ,
                   './V2/markers_convsmet.txt',sep='\t',
                   quote = F,row.names = F)


colnames(de.gene)[8]<-'cluster'

table(seu.TAM$TAM_cluster)
pdf('./Figuresubmit/c6_metvscon.pdf',height=4,width=5)
ggplot(markers[markers$avg_log2FC>=0,], 
       aes(x = avg_log2FC, 
           y = -log10(p_val), 
           fill = avg_log2FC)) +
  geom_jitter(size=2,shape=21,color='gray') +
  scale_fill_gradient2(low = '#0099CC',
                       mid = "white",
                       high = '#CC3333')+
  ggrepel::geom_text_repel(inherit.aes = F,
                           data = subset(markers[markers$avg_log2FC>=0,] , 
                                         gene%in%DAM_features),
                           mapping = aes(x = avg_log2FC,
                                         y = -log10(p_val),
                                         label = gene),
                           nudge_x = .1, nudge_y = .01,
                           max.overlaps = Inf,
                           color='black',
                           size = 6) +
  geom_jitter(data = subset(markers[markers$avg_log2FC>=0,], 
                            gene%in%DAM_features),
              
              aes(x = avg_log2FC,
                  y = -log10(p_val),
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

pdf('./FigureSubmit/controlvsmet.pdf',height=10,width=10)
DefaultAssay(seu.TAM)<-'RNA'

FeaturePlot_scCustom(seurat_object = seu.TAM, reduction='umap',
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),
                     features = c("Spp1"),split.by = 'type')&
  ( theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
          aspect.ratio = 1,
          legend.text = element_text(size=10),  # 您可以根据需要更改字体大小
          legend.title = element_text(size=10),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text = element_text(color = "black",size=12)))

dev.off()
##Enrichment-----
seu@meta.data<-seu@meta.data %>% 
  left_join(seu.TAM@meta.data[,c('sample','TAM_cluster')]) %>% 
  mutate(TAM_cluster=ifelse(is.na(TAM_cluster),'non_TAM',TAM_cluster))
rownames(seu@meta.data)<-seu$sample
DimPlot(seu,label=T,group.by = 'TAM_cluster')
DefaultAssay(seu)<-'RNA'

Idents(seu)<-seu$TAM_cluster
table(seu$TAM_cluster)
markers_all<-FindAllMarkers(seu)

markers_all_view<-markers_all%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(avg_log2FC))%>% 
  ungroup() %>% 
  filter(!(cluster=='non_TAM') & abs(avg_log2FC)>1)
cluster<-unique(markers_all_view$cluster)
library(clusterProfiler);library(org.Mm.eg.db)

markers_all_view$cluster<-as.character(markers_all_view$cluster)
er.list<-lapply(cluster, function(ct){
  de<-markers_all_view %>% 
    filter(cluster==ct) 
  gene <- markers_all_view$gene
  gene <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  gene_df <- data.frame(avg_diff=de$avg_log2FC,
                        SYMBOL = de$gene)
  gene_df <- merge(gene_df,gene,by="SYMBOL")  
  geneList <- gene_df$avg_diff
  names(geneList) = gene_df$ENTREZID
  geneList = sort(geneList, decreasing = TRUE)
  
  ego <- gseGO(geneList     = geneList,
               OrgDb        = org.Mm.eg.db,
               
               pvalueCutoff = 0.05,
               verbose      = FALSE)
  ego<-ego@result
  if (nrow(ego)!=0) {
    ego$celltype<-ct
    #ekegg<-ekegg@result
    return(ego)
    
  }
  
  
})

filtered_list <- er.list[!sapply(er.list, is.null)]
er.df_TAM<- purrr::reduce(filtered_list,rbind) %>% 
  arrange(celltype,desc(NES))
data<-er.df_TAM[,c('Description','NES','pvalue','p.adjust','celltype')] %>% 
  #mutate(celltype=str_extract(celltype,'\\d+')) %>% 
  filter(celltype%in%c('c2','c3','c4','c5','c6'))%>% 
  group_by(celltype) %>% 
  arrange(desc(NES)) %>% 
  #dplyr::slice(1:5) %>% 
  dplyr::slice(c(1:5, tail(row_number(), 5))) %>%
  ungroup() %>% 
  mutate(change=ifelse(NES>0,'up','down')) %>% 
  mutate(ID=paste0(Description,celltype),
         change=factor(change,levels=c('up','down'))) %>% 
  arrange(change,celltype,desc(abs(NES)))

data$ID<-factor(data$ID,levels = rev(data$ID))

pdf('./submit/dataset3/03inflamedTAM_enrichment.pdf',height=10,width=15)
ggplot(data,aes(x = ID, y = NES)) +
  coord_flip() +
  geom_segment(aes(x =  ID, xend =  ID,
                   y = 0, yend = NES, 
                   linetype = change,
                   color = celltype), 
               size = 1) + 
  geom_point(aes(size = 2, 
                 color = celltype,
                 shape=change)) + #添加散点/气泡
  scale_shape_manual(values =c(19,15))+
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  labs(x = NULL, y = 'NES',
       title = "GO Pathway Enrichment",
       color = "TAM cluster",
       shape= "change") +
  scale_y_continuous(expand = expansion(add = c(0.1, 0.1)),
                     limits = c(-8,8), breaks = seq(-8,8,2))+
  scale_color_manual(values=color.TAM_cluster[unique(data$celltype)])+ theme_bw() + 
  geom_text(data = data[data$change=='up',],
            aes(x = ID,y = -0.1,label = Description),
            hjust = 1,color = 'black',size = 6) +
  geom_text(data = data[data$change=='down',],
            aes(x = ID,y =0.1,label = Description),
            hjust = 0,color = 'black',size = 6) +
  scale_x_discrete(labels = NULL) 
dev.off()