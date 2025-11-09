
rm(list=ls())
library(Seurat);library(tidyverse);library(Nebulosa);library(clustree);library(GEOquery)
library(scCustomize)
source("~/Brain/brian/Rscript/DimReduction.R")
source("~/Brain/brian/Rscript/Lineage.R")
##Regular----
seu.TAM<-qs::qread('./02TAM/seu.TAM.qs')
seu.MDM<-subset(seu.TAM,celltype%in%c('APOE.MDM.DAM','APOE.MDM'))
seu.MDM<-subset(seu.TAM,celltype%in%c('APOE.MDM.DAM','APOE.MDM'))
DefaultAssay(seu.MDM)<-'RNA'
seu.MDM <- NormalizeData(seu.MDM)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
seu.MDM <- CellCycleScoring(seu.MDM, s.features = s.genes, g2m.features = g2m.genes)
seu.MDM <- FindVariableFeatures(seu.MDM, nfeatures = 2000)
seu.MDM <- RunPCA(seu.MDM, npcs = 300) 
seu.MDM <- RunHarmony(seu.MDM, group.by.vars = c("title"), reduction.use = "pca", dims.use = 1:300)
set.seed(13) 
seu.MDM  <- RunFDG(seu.MDM , reduction = "harmony", dims = 1:200)

color.MDM<-c( "#B49D99",#APOE.MDM
              "#98D277",#APOE.MDM.DAM
              "#2D5D8A" )
names(color.MDM)<-c('APOE.Mac','DAM-like.Mac','int')
table(seu.MDM$TAM_cluster)
table(seu.MDM$TAM_cluster,seu.MDM$celltype)
seu.MDM<-subset(seu.MDM,TAM_cluster=='TAM.0',invert=T)
d<-c('TAM.2'='APOE.Mac',
     'TAM.1'='int',
     'TAM.4'='DAM-like.Mac',
     'TAM.9'='int')
seu.MDM$MDMtype<-d[seu.MDM$TAM_cluster]
qs::qsave(seu.MDM, "./06trajectory.new/Brain.seu.MDM.trajectory.qs")

seu <- seu.MDM[rownames(seu.MDM[["RNA"]]@scale.data), ]
sceasy::convertFormat(seu, from = "seurat", to = "anndata", main_layer = "scale.data", 
                      outFile =  "./06trajectory.new/Brain.seu.MDM.scaled.h5ad")

##Import palantir results-----
palantir.res<-data.table::fread('./06trajectory.new/palantir_results.csv') 
seu.MDM@meta.data$ [,c('V1','pseudotime','DP','APOE.MDM.DAM','APOE.MDM.inflamed')]<-NULL
meta.data<-seu.MDM@meta.data %>% 
  rownames_to_column('V1') %>% 
  left_join(palantir.res[,c('V1','pseudotime','DP','APOE.MDM.DAM')])
rownames(meta.data)<-meta.data$V1
seu.MDM@meta.data<-meta.data
colnames(palantir.res)
pdf('./06trajectory.new/MDM.pseudotime.pdf',height=4,width=8)
p1<-FeaturePlot_scCustom(seu.MDM,reduction='fr',features=c('pseudotime'))
p2<-(DimPlot(seu.MDM,group.by = 'MDMtype',reduction = 'fr',label=T)+NoLegend())&ggsci::scale_color_d3()
p1|p2
dev.off()
##Lineage Figure6C-----
data.use <- FetchData(seu.MDM, vars = c("MDMtype", "FDG_1", "FDG_2"))
data.use <- data.use %>%
  mutate(cellID = rownames(.)) %>%
  group_by(MDMtype) %>%
  mutate(x = abs(FDG_1 - median(FDG_1)),
         y = abs(FDG_2 - median(FDG_2)),
         d = x + y) %>%
  arrange(MDMtype, d) %>%
  slice_head(n=1)
cell.embeddings <- Embeddings(seu.MDM, reduction = "fr")[data.use$cellID, ]
cell.pair.dis <- dist(cell.embeddings)
g <- igraph::graph_from_adjacency_matrix(adjmatrix = as.matrix(cell.pair.dis),
                                         mode = "undirected",
                                         weighted = TRUE,
                                         add.colnames = "cellID")

nodes2cellID <- igraph::get.vertex.attribute(g)$cellID
cellID2nodes <- setNames(1:length(nodes2cellID), nodes2cellID)

mst <- igraph::minimum.spanning.tree(g)
edges <- igraph::ends(mst, igraph::E(mst))
edges <- as.data.frame(edges)
edges$from <- nodes2cellID[edges$V1]
edges$to <- nodes2cellID[edges$V2]
edges$from.x <- cell.embeddings[edges$from, 1]
edges$from.y <- cell.embeddings[edges$from, 2]
edges$to.x <- cell.embeddings[edges$to, 1]
edges$to.y <- cell.embeddings[edges$to, 2]


root.cells <- "GGGCATCTCTCCCTGA_NS_19"
terminal.cells <- c(
  "APOE.MDM.DAM" = "CCGTTCATCCCAAGAT_NS_02"
)

nodes2cluster <- seu.MDM$MDMtype[nodes2cellID]
cluster2nodes <- setNames(names(nodes2cluster), nodes2cluster)

terminal.nodes <- cellID2nodes[cluster2nodes[seu.MDM$MDMtype[terminal.cells]]]
names(terminal.nodes) <- names(terminal.cells)

root.nodes <- cellID2nodes[cluster2nodes[seu.MDM$MDMtype[root.cells]]]

result <- igraph::shortest_paths(mst, from = root.nodes, to = terminal.nodes) ## modify me
lineages <- result$vpath
names(lineages) <- names(terminal.nodes)
seu.MDM$lineage.APOE.MDM.DAM <- seu.MDM$MDMtype %in% nodes2cluster[lineages$APOE.MDM.DAM]
## fit the principle curve
p.curve.APOE.MDM.DAM <- fit_pc(seu.MDM, lineage = "lineage.APOE.MDM.DAM", reduction = "fr", sample.n = 100) # 3
## plot on FR
data.point <- FetchData(seu.MDM, vars = c(paste0("FDG_", 1:2), "MDMtype"))
data.path.APOE.MDM.DAM <- get_path(p.curve.APOE.MDM.DAM, df=6)
data.arrow.APOE.MDM.DAM <- get_arrow(data.path.APOE.MDM.DAM, reverse = F)
p1<-ggplot() +
  geom_point(data = data.point, aes(FDG_1, FDG_2, color = MDMtype), size = .2) +
  geom_path(data = data.path.APOE.MDM.DAM, aes(X,Y), size = .8) +
  geom_segment(data = data.arrow.APOE.MDM.DAM, aes(x = X, xend = Xend, y = Y, yend = Yend),
               arrow = arrow(length = unit(0.1, "in"), angle = 30, type = "closed"), size = .5) +
  theme_classic(base_size = 15)+
  scale_color_manual(values=color.MDM)+
  theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
        aspect.ratio = 1,
        legend.text = element_text(size=10),  # 您可以根据需要更改字体大小
        legend.title = element_text(size=10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text = element_text(color = "black",size=12))+
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  ggtitle('Trajectory')
p1
table(is.na(seu.MDM$pseudotime))
p2<-FeaturePlot_scCustom(seu.MDM,
                         reduction='fr',
                         colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),
                         features=c('pseudotime'))+
  theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
        aspect.ratio = 1,
        legend.title = element_blank(),
        legend.text = element_text(size=10), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text = element_text(color = "black",size=12))+
  ggtitle('Pseudotime')
p2

# cells around the curve
# seu.MDM <- cells_on_lineage(seu.MDM, lineage = "APOE.MDM.DAM", reduction = "fr", data.path = data.path.APOE.MDM.DAM, delta = 0.1)
# p2<-DimPlot(seu.MDM, reduction = "fr", group.by = "lineage.finetune.APOE.MDM.DAM") +
#   scale_color_manual(values = c("grey", "red"))+NoLegend()
pdf('./06trajectory.new/MDM.pseudo.trajectory.pdf',height = 4,width=10)
p1|p2
dev.off()

##identityFigure6D-F-----
pdf('./06trajectory.new/maco_pheno.pdf',height=5,width=10)
DefaultAssay(seu.MDM)<-'maco_pheno'
FeaturePlot_scCustom(seurat_object = seu.MDM, reduction='fr',
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),
                     features = c("Angiogenesis","Phagocytosis"))&
  ( theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
          aspect.ratio = 1,
          legend.text = element_text(size=10),  # 您可以根据需要更改字体大小
          legend.title = element_text(size=10),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text = element_text(color = "black",size=12)))
dev.off()
pdf('./06trajectory.new/TAM.features.pdf',height=5,width=15)
DefaultAssay(seu.MDM)<-'TAM.features'
FeaturePlot_scCustom(seurat_object = seu.MDM, reduction='fr',
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),
                     features = c('APC-Topic-Score','IFN-Response-Topic-Score',
                                  'Secretory-Topic-Score'), num_columns = 3)&
  ( theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
          aspect.ratio = 1,
          legend.text = element_text(size=10),  # 您可以根据需要更改字体大小
          legend.title = element_text(size=10),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text = element_text(color = "black",size=12)))
dev.off()
pdf('./06trajectory.new/DAM.features.pdf',height=5,width=15)
DefaultAssay(seu.MDM)<-'DAM'
FeaturePlot_scCustom(seu.MDM,features=c("DAM-core-program", 
                                        "Stage1-DAM-upregulated",
                                        "Stage2-DAM"),reduction='fr', 
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100), num_columns = 3)&
  ( theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
          aspect.ratio = 1,
          legend.text = element_text(size=10),  # 您可以根据需要更改字体大小
          legend.title = element_text(size=10),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text = element_text(color = "black",size=12)))
dev.off()
pdf('./06trajectory.new/DAM.Gene.pdf',height=12,width=12)

DefaultAssay(seu.MDM)<-'RNA'
FeaturePlot_scCustom(seu.MDM, reduction = "fr",  num_columns = 3,
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),
                     features = c('TREM2','TYROBP','APOE',
                                  'SPP1','OLR1','CCL3',
                                  'CCL4','CD83','CSF1R'))
FeaturePlot_scCustom(seu.MDM, reduction = "fr",  num_columns = 3,
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),
                     features = c('MS4A7','LILRB4', 'VSIG4',
                                  'FCER1G','FCGR2A',
                                  'MS4A4A','MS4A6A','CCL18','CAPG'))
FeaturePlot_scCustom(seu.MDM, reduction = "fr",  num_columns = 3,
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),
                     features = c('CTSD','CTSB','CTSL',
                                  'CTSS','PSAP','DAB2',
                                  'FTL',"GRN",'CD74','HLA-DPB1'))

dev.off()
pdf('./06trajectory.new/DAM.Gene.aveg.pdf',height=12,width=12)
Idents(seu.MDM)<-factor(seu.MDM$MDMtype,levels = c('APOE.Mac','int','DAM-like.Mac'))
AverageHeatmap(object = seu.MDM ,
               assays = 'RNA',
               markerGene =c('APOE','TREM2','TYROBP',
                             'SPP1','OLR1','CCL3','CCL4',
                             'CD83','CSF1R','MS4A7',
                             'LILRB4', 'VSIG4',
                             'FCER1G','FCGR2A',
                             'MS4A4A','MS4A6A','CCL18','CAPG',
                             'CTSD','CTSB','CTSL','CTSS',
                             'PSAP','DAB2','FTL',"GRN"),
               clusterAnnoName = F,
               myanCol = color.MDM[levels(Idents(seu.MDM))],
               annoCol = T,height=10,width=2)

dev.off()
##MDM_markers Figure6G------
Idents(seu.MDM)<-seu.MDM$TAM_cluster
markers_DAM.MDM<-FindMarkers(seu.MDM,ident.1 = 'TAM.4',ident.2 = 'TAM.2')
markers_DAM.MDM<-markers_DAM.MDM %>% 
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(desc(avg_log2FC))%>%
  rownames_to_column('gene')


data.table::fwrite(markers_DAM.MDM,
                   './06trajectory/DAM_vs_Orig.txt',sep='\t',
                   quote = F,row.names = F)

MDM_key<-readRDS('~/Brain/brian/DAM/DAM.Mac.Gene.rds')
MDM_markers<-data.table::fread('06trajectory/DAM_vs_Orig.txt')
pdf('./Figuresubmit/DAM.Macmarkers_Trajectory.pdf',height=6,width=6)
MDM_markers$p_val_adj[MDM_markers$p_val_adj==0]<-NA
MDM_markers$p_val_adj[is.na(MDM_markers$p_val_adj)]<-min(MDM_markers$p_val_adj,na.rm=T)

ggplot(MDM_markers[MDM_markers$avg_log2FC>=0,], 
       aes(x=avg_log2FC, 
           y=-log10(p_val_adj),
           fill=avg_log2FC)) +
  geom_jitter(size = 2, shape = 21, color = 'gray', width = 0.2, height = 0.2) +
  scale_fill_gradient2(low = '#0099CC',
                       mid = "white",
                       high = '#CC3333')+
  ggrepel::geom_text_repel(inherit.aes = F,
                           data = subset(MDM_markers, 
                                         gene%in%MDM_key),
                           mapping = aes(x=avg_log2FC,
                                         y=-log10(p_val_adj),
                                         label = gene),
                           nudge_x = .1, nudge_y = .01,
                           max.overlaps = Inf,
                           color='black',
                           size = 6) +
  geom_jitter(data = subset(MDM_markers, 
                            gene%in%MDM_key),
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


##Key regulonFigure6H-K------
seu.MDM<-readRDS('./06trajectory.new/seu.MDM.addScenic.qs')
DAM.Mac.Regulon<-readRDS('~/Brain/brian/DAM/DAM.MDM.Terminal.regulon.rds')
pdf('./06trajectory.new/Regulon.pdf',height=10,width=10)
DefaultAssay(seu.MDM)<-'AUCell'
FeaturePlot_scCustom(seurat_object = seu.MDM, reduction='fr',
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),
                     features = paste0(DAM.Mac.Regulon,'(+)')[1:4])&
  ( theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
          aspect.ratio = 1,
          legend.text = element_text(size=10),  # 您可以根据需要更改字体大小
          legend.title = element_text(size=10),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text = element_text(color = "black",size=12)))
FeaturePlot_scCustom(seurat_object = seu.MDM, reduction='fr',
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),
                     features = paste0(DAM.Mac.Regulon,'(+)')[5:7])&
  ( theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
          aspect.ratio = 1,
          legend.text = element_text(size=10),  # 您可以根据需要更改字体大小
          legend.title = element_text(size=10),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text = element_text(color = "black",size=12)))
dev.off()

pdf('./06trajectory.new/Regulon_RNA.pdf',height=10,width=10)
DefaultAssay(seu.MDM)<-'RNA'
FeaturePlot_scCustom(seurat_object = seu.MDM, reduction='fr',
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),
                     features = DAM.Mac.Regulon[1:4])&
  ( theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
          aspect.ratio = 1,
          legend.text = element_text(size=10),  # 您可以根据需要更改字体大小
          legend.title = element_text(size=10),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text = element_text(color = "black",size=12)))
FeaturePlot_scCustom(seurat_object = seu.MDM, reduction='fr',
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),
                     features = DAM.Mac.Regulon[5:7])&
  ( theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
          aspect.ratio = 1,
          legend.text = element_text(size=10),  # 您可以根据需要更改字体大小
          legend.title = element_text(size=10),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text = element_text(color = "black",size=12)))
dev.off()
saveRDS(seu.MDM,'./06trajectory.new/seu.MDM.addScenic.qs')

##Findmarkers new----
rm(list=ls());gc();
seu.TAM<-readRDS("./Figuresubmit/seu.TAM.celltype.pl.rds")
seu<-qs::qread('./01seuRegular/seu.anno.qs')
table(seu.TAM$celltype_pl,seu.TAM$TAM_cluster)
seu.TAM$celltype_tra<-seu.TAM$celltype_pl
seu.TAM$celltype_tra[seu.TAM$TAM_cluster=='TAM.4']<-'DAM-like.Mac.Terminal'
seu.TAM$celltype_tra[seu.TAM$TAM_cluster=='TAM.2']<-'APOE.Mac.Original'
table(seu.TAM$celltype_tra,seu.TAM$TAM_cluster)
DefaultAssay(seu.TAM)<-'RNA'
Idents(seu.TAM)<-seu.TAM$celltype_tra
TAM_markers<-FindAllMarkers(seu.TAM)
TAM_markers <-TAM_markers%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(avg_log2FC))%>%
  ungroup()
seu$sample<-rownames(seu@meta.data)
seu@meta.data<-seu@meta.data %>% 
  left_join(seu.TAM@meta.data[,c('sample','celltype_tra')]) %>% 
  mutate(celltype_tra=ifelse(is.na(celltype_tra),
                             'others',celltype_tra))
table(seu$celltype_tra)
rownames(seu@meta.data)<-seu@meta.data$sample
Idents(seu)<-seu$celltype_tra
DefaultAssay(seu)<-'RNA'
TAM_markers2<-FindAllMarkers(seu)
TAM_markers2 <-TAM_markers2%>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(cluster,desc(avg_log2FC))%>%
  ungroup() %>% 
  filter(!(cluster=='others'))
data.table::fwrite(TAM_markers,
                   './06trajectory.new/TAMcelltypetra_markers.txt',sep='\t',
                   quote = F,row.names = F)
data.table::fwrite(TAM_markers2,
                   './06trajectory.new/TAMcelltypetra_markers_all.txt',sep='\t',
                   quote = F,row.names = F)
qs::qsave(seu.TAM,'./06trajectory.new/seu.TAM.tra.qs')

##enrichmentFigure6B------
geneList <- TAM_markers2$avg_log2FC[TAM_markers2$cluster=='DAM-like.Mac.Terminal']
names(geneList) = TAM_markers2$gene[TAM_markers2$cluster=='DAM-like.Mac.Terminal']

geneList = sort(geneList, decreasing = TRUE)
library(clusterProfiler)
library(enrichplot)
hallmarks <- read.gmt("~/Brain/brian/Geneset/h.all.v2022.1.Hs.symbols.gmt")
y <- GSEA(geneList, TERM2GENE = hallmarks)
yd <- as.data.frame(y)
dotplot(y,showCategory=30,split=".sign") + facet_grid(~.sign)

saveRDS(y@result,'./06trajectory.new/DAM.Mac_enrichment_GSVA.rds')
saveRDS(y@result,'./06trajectory.new/APOE.orig_enrichment_GSVA.rds')
e.regulon<-readRDS('./06trajectory.new/DAM.Mac_enrichment_GSVA.rds')
DAM.enrich<-readRDS('~/Brain/brian/DAM/DAM.mac.tra.enrichent.rds')
enrichment_data<-e.regulon %>%
  filter(ID%in%DAM.enrich) %>% 
  arrange(desc(NES))


enrichment_data$ID<-factor(enrichment_data$ID,levels=rev(DAM.enrich))
pdf('./Figuresubmit/DAM.mac.enrich.pdf',height=3,width=8)
ggplot(enrichment_data,
       aes(x = ID, y= NES, fill =-log10(p.adjust)) ) +
  geom_bar(stat = "identity")+  # 用 size 映射 p 值
  #scale_size(range = c(3, 10), name = "-log10(p.adjust)") +  
  scale_fill_gradient(low = '#0099CC',high = '#CC3333')+
  labs(x = "", y = "NES", title = "") +
  theme_bw() +
  coord_flip()+
  theme(
    #aspect.ratio = 1,
    legend.text = element_text(size=12),  
    plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
    axis.text = element_text(color = "black",size=14),
    axis.title = element_text(color = "black",size=14),
    legend.position = "right",
    legend.title = element_blank())
dev.off()


