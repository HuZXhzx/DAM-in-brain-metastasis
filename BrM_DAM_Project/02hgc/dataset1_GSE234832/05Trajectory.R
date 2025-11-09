rm(list=ls())
library(Seurat);library(tidyverse);library(Nebulosa);library(clustree);library(GEOquery)
library(scCustomize)
##TrajectoryRegular-----
source("~/Brain/brian/Rscript/DimReduction.R")
source("~/Brain/brian/Rscript/Lineage.R")
seu.TAM<-qs::qread('./dataset1/Lung/04clusterMDM/seu.TAM.qs')
color.MDM<-c( "#B49D99",#APOE.MDM
              "#98D277",#APOE.MDM.DAM
              "#2D5D8A" )
names(color.MDM)<-c('APOE.Mac','DAM-like.Mac','int')
table(seu.MDM$TAM_cluster)
d<-c('TAM.1'='APOE.Mac',
     'TAM.4'='int',
     'TAM.0'='DAM-like.Mac',
     'TAM.3'='DAM-like.Mac')
seu.MDM$MDMtype<-d[seu.MDM$TAM_cluster]
seu.MDM<-subset(seu.MDM,TAM_cluster%in%c('TAM.0','TAM.1','TAM.4'))
Idents(seu.MDM)<-factor(seu.MDM$MDMtype,levels = c('APOE.Mac','int','DAM-like.Mac'))

table(seu.MDM$celltype.new)
DefaultAssay(seu.MDM)<-'RNA'
seu.MDM <- NormalizeData(seu.MDM)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
seu.MDM <- CellCycleScoring(seu.MDM, s.features = s.genes, g2m.features = g2m.genes)
seu.MDM <- FindVariableFeatures(seu.MDM, nfeatures = 2000)

seu.MDM <- ScaleData(seu.MDM, vars.to.regress = c("S.Score", "G2M.Score"))


seu.MDM <- RunPCA(seu.MDM, npcs = 300) 
seu.MDM <- RunHarmony(seu.MDM, group.by.vars = c("title"), reduction.use = "pca", dims.use = 1:300)
set.seed(13) 
seu.MDM  <- RunFDG(seu.MDM , reduction = "harmony", dims = 1:200)


seu <- seu.MDM[rownames(seu.MDM[["RNA"]]@scale.data), ]
sceasy::convertFormat(seu, from = "seurat", to = "anndata", main_layer = "scale.data", 
                      outFile =  "./Lung/06trajectory.new/Brain.seu.MDM.scaled.h5ad")


##Import palantir results-----
palantir.res<-data.table::fread('./Lung/06trajectory.new/palantir_results.csv') 
seu.MDM@meta.data[,c('V1','pseudotime','DP','APOE.MDM.DAM')]<-NULL
meta.data<-seu.MDM@meta.data %>% 
  rownames_to_column('V1') %>% 
  left_join(palantir.res[,c('V1','pseudotime','DP','APOE.MDM.DAM')])
rownames(meta.data)<-meta.data$V1
seu.MDM@meta.data<-meta.data
colnames(palantir.res)

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

root.cells <- "GSM7475325_GTTCGCTGTGACTATC-1"
terminal.cells <- c(
  "APOE.MDM.DAM"="GSM7475326_TTTCCTCCACATTCTT-1"
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
data.path.APOE.MDM.DAM <- get_path(p.curve.APOE.MDM.DAM, df=5)
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
p3<-(DimPlot(seu.MDM,group.by = 'MDMtype',reduction = 'fr',label=T)+NoLegend())&ggsci::scale_color_d3()

# cells around the curve
# seu.MDM <- cells_on_lineage(seu.MDM, lineage = "APOE.MDM.DAM", reduction = "fr", data.path = data.path.APOE.MDM.DAM, delta = 0.1)
# p2<-DimPlot(seu.MDM, reduction = "fr", group.by = "lineage.finetune.APOE.MDM.DAM") +
#   scale_color_manual(values = c("grey", "red"))+NoLegend()
pdf('./Lung/06trajectory.new/MDM.pseudo.trajectory.pdf',height = 4,width=10)
p1|p2
dev.off()

##identity Figure6D-F-----
pdf('./Lung/06trajectory.new/maco_pheno.pdf',height=5,width=10)
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
pdf('./Lung/06trajectory.new/TAM.features.pdf',height=5,width=15)
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
pdf('./Lung/06trajectory.new/DAM.features.pdf',height=5,width=15)
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




##Key regulonFigure6H-K------
DAM.Mac.Regulon<-readRDS('~/Brain/brian/DAM/DAM.Mac.Trajectory.Regulon.rds')
pdf('./Lung/06trajectory.new/Regulon.pdf',height=10,width=10)
DefaultAssay(seu.MDM)<-'AUCell'
FeaturePlot_scCustom(seurat_object = seu.MDM, reduction='fr',
                     colors_use = colorRampPalette(c('#0099CC',"white",'#CC3333'))(100),
                     features = DAM.Mac.Regulon)&
  ( theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
          aspect.ratio = 1,
          legend.text = element_text(size=10),  # 您可以根据需要更改字体大小
          legend.title = element_text(size=10),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text = element_text(color = "black",size=12)))

dev.off()




##enrichmentFigure6B------
geneList <- TAM_markers2$avg_log2FC[TAM_markers2$cluster=='DAM-like.Mac.Terminal']
names(geneList) = TAM_markers2$gene[TAM_markers2$cluster=='DAM-like.Mac.Terminal']

geneList = sort(geneList, decreasing = TRUE)
library(clusterProfiler)
library(enrichplot)
hallmarks <- read.gmt("../Geneset/h.all.v2022.1.Hs.symbols.gmt")
y <- GSEA(geneList, TERM2GENE = hallmarks)
yd <- as.data.frame(y)
dotplot(y,showCategory=30,split=".sign") + facet_grid(~.sign)

saveRDS(y@result,'./Lung/06trajectory.new/DAM.Mac_enrichment_GSVA.rds')

e.regulon<-readRDS('./dataset1/Lung/06trajectory.new/DAM.Mac_enrichment_GSVA.rds')
DAM.enrich<-readRDS('./DAM/DAM.mac.tra.enrichent.rds')
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



