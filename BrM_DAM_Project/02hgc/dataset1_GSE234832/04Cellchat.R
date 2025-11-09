rm(list=ls());gc()
library(GEOquery);library(tidyverse);library(arrayQualityMetrics);library(limma)
library(Seurat);library(harmony);library(scCustomize);library(scRNAtoolVis);library(clusterProfiler);library(UCell)
library(scCustomize);library(scRNAtoolVis);library(CellChat)
seu.TAM<-readRDS('./Figuresubmit/seu.TAM.rds')
##Cellchat between tumor and TAMs----
seu.mtc<-qs::qread('./Lung/05cellchat/01TumorCluster/seu.mtc.qs')
table(seu.mtc$MTC_cluster)

genes<-intersect(rownames(seu.mtc@assays$RNA) ,rownames(seu.TAM@assays$RNA))

xcounts<-cbind(seu.mtc@assays$RNA@counts[genes,],
               seu.TAM@assays$RNA@counts[genes,])

meta.data.mtc<-seu.mtc@meta.data[,c("sample","orig.ident","title","group","type","MTC_cluster")]
meta.data.TAM<-seu.TAM@meta.data[,c("sample","orig.ident","title","group","type","celltype_pl")]
colnames(meta.data.mtc)[6]<-c('cluster')
colnames(meta.data.TAM)[6]<-c('cluster')
meta.data<-rbind(meta.data.mtc,meta.data.TAM)
rownames(meta.data)
table(meta.data$cluster)
seu <- CreateSeuratObject(counts = xcounts,
                          meta.data = meta.data)
seu <- NormalizeData(seu)

data.input = seu[["RNA"]]@data # normalized data matrix
Idents(seu)<-seu$cluster
labels <- Idents(seu)
library(CellChat)
cc.obj <- createCellChat(object = seu, 
                         group.by = "ident", 
                         assay = "RNA")

cc.obj  <- addMeta(cc.obj , meta = meta.data)
cc.obj <- setIdent(cc.obj, ident.use = "cluster") 

levels(cc.obj@idents) 
groupSize <- as.numeric(table(cc.obj@idents))
CellChatDB <- CellChatDB.human
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB 
cc.obj@DB <- CellChatDB.use
cc.obj<- subsetData(cc.obj)
future::plan("multisession", workers = 4) # do parallel
cc.obj <- identifyOverExpressedGenes(cc.obj)
cc.obj <- identifyOverExpressedInteractions(cc.obj)
cc.obj <- computeCommunProb(cc.obj, type = "triMean")
cc.obj <- filterCommunication(cc.obj, min.cells = 10)
cc.obj<- computeCommunProbPathway(cc.obj)
cc.obj<- aggregateNet(cc.obj)
df.net <- subsetCommunication(cc.obj) 

dir.create('./Lung/05cellchat/new/')
qs::qsave(cc.obj,'./Lung/05cellchat/new/cc.obj_mtc.new.qs')
data.table::fwrite(df.net,
                   './Lung/05cellchat/new/df_net_DAM_MTC.new.txt',sep='\t',
                   quote = F,row.names = F)

##chord plot Figure7B Tumor->TAMs----
cc.obj<-qs::qread('./Lung/05cellchat/new/cc.obj_mtc.new.qs')
df.net<-data.table::fread('./Lung/05cellchat/new/df_net_DAM_MTC.new.txt')

prob <- slot(cc.obj, "net")$prob
pval <- slot(cc.obj, "net")$pval
prob[pval > 0.05] <- 0
net <- reshape2::melt(prob, value.name = "prob")
colnames(net)[1:3] <- c("source","target","interaction_name")
cols.default <- c("interaction_name_2", "pathway_name",  "ligand",  "receptor" ,"annotation","evidence")
cols.common <- intersect(cols.default,colnames(cc.obj@LR$LRsig))
pairLR = dplyr::select(cc.obj@LR$LRsig, cols.common)
idx <- match(net$interaction_name, rownames(pairLR))
temp <- pairLR[idx,]
net <- cbind(net, temp)
#signaling<-c("ApoE","APP")
DAM.MTC.pathway<-c("APP - CD74",
                   "APP - (TREM2+TYROBP)",
                   "MDK - LRP1", 
                   "MDK - NCL",
                   "CD99 - CD99",
                   "CD99 - PILRA",
                   "ICAM1 - (ITGAM+ITGB2)",
                   "LAMA5 - CD44",
                   "COL1A1 - CD44")


net <- subset(net, interaction_name_2 %in% DAM.MTC.pathway)
table(cc.obj@idents)
targets.use = c('DAM','DAM-like.Mac'); 
sources.use = c('EMT_high',"REP_high",'EMT_low')
net <- subset(net, source %in% sources.use)
net <- subset(net, target %in% targets.use)
net$source<-as.character(net$source)
net$source[net$source%in%c('EMT_high',"REP_high",'EMT_low')]<-'Tumor'
sources.use = c('Tumor')
levels(cc.obj@idents)[levels(cc.obj@idents)%in%c('EMT_high',"REP_high",'EMT_low')]<-'Tumor'

df <- subset(net, prob > 0)
df$id<-paste(df$source,df$target,df$interaction_name,sep='.')
df<-df %>% 
  arrange(id,desc(prob)) %>% 
  distinct(id,.keep_all = T)
df$id <- 1:nrow(df)
df$receptor[df$receptor=='CD99']<-'CD99-R'

ligand.uni <- unique(df$ligand)
for (i in 1:length(ligand.uni)) {
  df.i <- df[df$ligand == ligand.uni[i], ]
  source.uni <- unique(df.i$source)
  for (j in 1:length(source.uni)) {
    df.i.j <- df.i[df.i$source == source.uni[j], ]
    df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
    df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
  }
}

receptor.uni <- unique(df$receptor)
for (i in 1:length(receptor.uni)) {
  df.i <- df[df$receptor == receptor.uni[i], ]
  target.uni <- unique(df.i$target)
  for (j in 1:length(target.uni)) {
    df.i.j <- df.i[df.i$target == target.uni[j], ]
    df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
    df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
  }
}

cell.order.sources <- levels(cc.obj@idents)[levels(cc.obj@idents) %in% sources.use]
cell.order.targets <- levels(cc.obj@idents)[levels(cc.obj@idents) %in% targets.use]

df$source <- factor(df$source, levels = cell.order.sources)
df$target <- factor(df$target, levels = cell.order.targets)

df.ordered.source <- df[with(df, order(source, -prob)), ]
df.ordered.target <- df[with(df, order(target, -prob)), ]

order.source <- unique(df.ordered.source[ ,c('ligand','source')])
order.target <- unique(df.ordered.target[ ,c('receptor','target')])
order.sector <- c(order.source$ligand, order.target$receptor)
color.TAM_hgc<-readRDS('../DAM/color.TAM_hgc.rds')

color.use <- c("Tumor"="#6A3D9A",color.TAM_hgc)

# define edge color
edge.color <- color.use[as.character(df.ordered.source$source)]
names(edge.color) <- as.character(df.ordered.source$source)

grid.col.ligand <- color.use[as.character(order.source$source)]
names(grid.col.ligand) <- as.character(order.source$source)
grid.col.receptor <- color.use[as.character(order.target$target)]
names(grid.col.receptor) <- as.character(order.target$target)
grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
names(grid.col) <- order.sector
df.plot <- df.ordered.source[ ,c('ligand','receptor','prob')]
circos.clear()

pdf('./Figuresubmit/TAM&MTC2.pdf',height=10,width=10)  
chordDiagram(df.plot,
             order = order.sector,
             col = edge.color,
             grid.col = grid.col,
             transparency = 0.4,
             link.border = NA,
             directional = 1,
             direction.type = c("diffHeight","arrows"),
             link.arr.type = "big.arrow",
             annotationTrack = "grid",
             annotationTrackHeight = 0.03,
             preAllocateTracks = list(track.height = max(strwidth(order.sector))),
             small.gap = 1,
             big.gap = 10,
             link.visible = T,
             scale = F,
             link.target.prop = T,
             reduce = -1)

circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = 0.8)
}, bg.border = NA)

lgd <- ComplexHeatmap::Legend(at = names(color.use), type = "grid", legend_gp = grid::gpar(fill = color.use), title = "Cell State")
ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(20, "mm"), y = unit(20, "mm"), just = c("right", "bottom"))

circos.clear()
gg <- recordPlot()
gg
dev.off()


##chord plot Figure7B  TAMs->Tumor-----
MTC.DAM.pathway<-c(
  "SPP1 - CD44",
  "SPP1 - (ITGAV+ITGB1)",
  "LGALS9 - P4HB",
  "LGALS9 - CD44",
  "TNF - TNFRSF1A",
  "TNFSF12 - TNFRSF12A",
  "GRN - SORT1")
cc.obj<-qs::qread('./Lung/05cellchat/new/cc.obj_mtc.new.qs')
df.net<-data.table::fread('./Lung/05cellchat/new/df_net_DAM_MTC.new.txt')

prob <- slot(cc.obj, "net")$prob
pval <- slot(cc.obj, "net")$pval
prob[pval > 0.05] <- 0
net <- reshape2::melt(prob, value.name = "prob")
colnames(net)[1:3] <- c("source","target","interaction_name")
cols.default <- c("interaction_name_2", "pathway_name",  "ligand",  "receptor" ,"annotation","evidence")
cols.common <- intersect(cols.default,colnames(cc.obj@LR$LRsig))
pairLR = dplyr::select(cc.obj@LR$LRsig, cols.common)
idx <- match(net$interaction_name, rownames(pairLR))
temp <- pairLR[idx,]
net <- cbind(net, temp)

net <- subset(net, interaction_name_2 %in% MTC.DAM.pathway)
table(cc.obj@idents)
sources.use = c('DAM','DAM-like.Mac'); 
targets.use = c('EMT_high',"REP_high",'EMT_low')
net <- subset(net, source %in% sources.use)
net <- subset(net, target %in% targets.use)
net$target<-as.character(net$target)
net$target[net$target%in%c('EMT_high',"REP_high",'EMT_low')]<-'Tumor'
targets.use = c('Tumor')
levels(cc.obj@idents)[levels(cc.obj@idents)%in%c('EMT_high',"REP_high",'EMT_low')]<-'Tumor'

df <- subset(net, prob > 0)
df$id<-paste(df$source,df$target,df$interaction_name,sep='.')
df<-df %>% 
  arrange(id,desc(prob)) %>% 
  distinct(id,.keep_all = T)
df$id <- 1:nrow(df)
df$receptor[df$receptor=='CD99']<-'CD99-R'

ligand.uni <- unique(df$ligand)
for (i in 1:length(ligand.uni)) {
  df.i <- df[df$ligand == ligand.uni[i], ]
  source.uni <- unique(df.i$source)
  for (j in 1:length(source.uni)) {
    df.i.j <- df.i[df.i$source == source.uni[j], ]
    df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
    df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
  }
}

receptor.uni <- unique(df$receptor)
for (i in 1:length(receptor.uni)) {
  df.i <- df[df$receptor == receptor.uni[i], ]
  target.uni <- unique(df.i$target)
  for (j in 1:length(target.uni)) {
    df.i.j <- df.i[df.i$target == target.uni[j], ]
    df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
    df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
  }
}

cell.order.sources <- levels(cc.obj@idents)[levels(cc.obj@idents) %in% sources.use]
cell.order.targets <- levels(cc.obj@idents)[levels(cc.obj@idents) %in% targets.use]

df$source <- factor(df$source, levels = cell.order.sources)
df$target <- factor(df$target, levels = cell.order.targets)

df.ordered.source <- df[with(df, order(source, -prob)), ]
df.ordered.target <- df[with(df, order(target, -prob)), ]

order.source <- unique(df.ordered.source[ ,c('ligand','source')])
order.target <- unique(df.ordered.target[ ,c('receptor','target')])
order.sector <- c(order.source$ligand, order.target$receptor)
color.TAM_hgc<-readRDS('../DAM/color.TAM_hgc.rds')

color.use <- c("Tumor"="#6A3D9A",color.TAM_hgc[c('DAM','DAM-like.Mac')])

# define edge color
edge.color <- color.use[as.character(df.ordered.source$source)]
names(edge.color) <- as.character(df.ordered.source$source)

grid.col.ligand <- color.use[as.character(order.source$source)]
names(grid.col.ligand) <- as.character(order.source$source)
grid.col.receptor <- color.use[as.character(order.target$target)]
names(grid.col.receptor) <- as.character(order.target$target)
grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
names(grid.col) <- order.sector
df.plot <- df.ordered.source[ ,c('ligand','receptor','prob')]
circos.clear()

pdf('./Figuresubmit/TAM&MTC2.pdf',height=10,width=10)  
chordDiagram(df.plot,
             order = order.sector,
             col = edge.color,
             grid.col = grid.col,
             transparency = 0.4,
             link.border = NA,
             directional = 1,
             direction.type = c("diffHeight","arrows"),
             link.arr.type = "big.arrow",
             annotationTrack = "grid",
             annotationTrackHeight = 0.03,
             preAllocateTracks = list(track.height = max(strwidth(order.sector))),
             small.gap = 1,
             big.gap = 10,
             link.visible = T,
             scale = F,
             link.target.prop = T,
             reduce = -1)

circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = 0.8)
}, bg.border = NA)

lgd <- ComplexHeatmap::Legend(at = names(color.use), type = "grid", legend_gp = grid::gpar(fill = color.use), title = "Cell State")
ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(20, "mm"), y = unit(20, "mm"), just = c("right", "bottom"))

circos.clear()
gg <- recordPlot()
gg
dev.off()
##chord plot Figure7C TAM----
TAM.pathway<-c("APOE - TREM2","SPP1 - CD44",
               "LGALS9 - CD44",
               "LGALS9 - P4HB",
               "LGALS9 - HAVCR2","LGALS9 - CD45")

net <- reshape2::melt(prob, value.name = "prob")
colnames(net)[1:3] <- c("source","target","interaction_name")
cols.default <- c("interaction_name_2", "pathway_name",  "ligand",  "receptor" ,"annotation","evidence")
cols.common <- intersect(cols.default,colnames(cc.obj@LR$LRsig))
pairLR = dplyr::select(cc.obj@LR$LRsig, cols.common)
idx <- match(net$interaction_name, rownames(pairLR))
temp <- pairLR[idx,]
net <- cbind(net, temp)
net <- subset(net, interaction_name_2 %in% TAM.pathway)

targets.use = c('DAM','DAM-like.Mac'); 
sources.use = c('DAM','DAM-like.Mac','APOE.Mac')
net <- subset(net, source %in% sources.use)
net <- subset(net, target %in% targets.use)
net$source<-as.character(net$source)

df <- subset(net, prob > 0)
df$id<-paste(df$source,df$target,df$interaction_name,sep='.')
df<-df %>% 
  arrange(id,desc(prob)) %>% 
  distinct(id,.keep_all = T)
df$id <- 1:nrow(df)
df$receptor[df$receptor=='CD99']<-'CD99-R'

ligand.uni <- unique(df$ligand)
for (i in 1:length(ligand.uni)) {
  df.i <- df[df$ligand == ligand.uni[i], ]
  source.uni <- unique(df.i$source)
  for (j in 1:length(source.uni)) {
    df.i.j <- df.i[df.i$source == source.uni[j], ]
    df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
    df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
  }
}

receptor.uni <- unique(df$receptor)
for (i in 1:length(receptor.uni)) {
  df.i <- df[df$receptor == receptor.uni[i], ]
  target.uni <- unique(df.i$target)
  for (j in 1:length(target.uni)) {
    df.i.j <- df.i[df.i$target == target.uni[j], ]
    df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
    df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
  }
}

cell.order.sources <- levels(cc.obj@idents)[levels(cc.obj@idents) %in% sources.use]
cell.order.targets <- levels(cc.obj@idents)[levels(cc.obj@idents) %in% targets.use]

df$source <- factor(df$source, levels = cell.order.sources)
df$target <- factor(df$target, levels = cell.order.targets)

df.ordered.source <- df[with(df, order(source, -prob)), ]
df.ordered.target <- df[with(df, order(target, -prob)), ]

order.source <- unique(df.ordered.source[ ,c('ligand','source')])
order.target <- unique(df.ordered.target[ ,c('receptor','target')])
order.sector <- c(order.source$ligand, order.target$receptor)

color.use <- color.TAM_hgc[c('DAM','DAM-like.Mac','APOE.Mac')]

# define edge color
edge.color <- color.use[as.character(df.ordered.source$source)]
names(edge.color) <- as.character(df.ordered.source$source)

grid.col.ligand <- color.use[as.character(order.source$source)]
names(grid.col.ligand) <- as.character(order.source$source)
grid.col.receptor <- color.use[as.character(order.target$target)]
names(grid.col.receptor) <- as.character(order.target$target)
grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
names(grid.col) <- order.sector
df.plot <- df.ordered.source[ ,c('ligand','receptor','prob')]
circos.clear()

pdf('./Figuresubmit/TAM&TAM.pdf',height=10,width=10)  
chordDiagram(df.plot,
             order = order.sector,
             col = edge.color,
             grid.col = grid.col,
             transparency = 0.4,
             link.border = NA,
             directional = 1,
             direction.type = c("diffHeight","arrows"),
             link.arr.type = "big.arrow",
             annotationTrack = "grid",
             annotationTrackHeight = 0.03,
             preAllocateTracks = list(track.height = max(strwidth(order.sector))),
             small.gap = 1,
             big.gap = 10,
             link.visible = T,
             scale = F,
             link.target.prop = T,
             reduce = -1)

circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = 0.8)
}, bg.border = NA)

lgd <- ComplexHeatmap::Legend(at = names(color.use), type = "grid", legend_gp = grid::gpar(fill = color.use), title = "Cell State")
ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(20, "mm"), y = unit(20, "mm"), just = c("right", "bottom"))

circos.clear()
gg <- recordPlot()
gg
dev.off()


###Prepare DIY dotplot----
rm(list=ls())
TAM.pathway<-c("APOE - TREM2","SPP1 - CD44",
               "LGALS9 - CD44",
               "LGALS9 - P4HB",
               "LGALS9 - HAVCR2","LGALS9 - CD45",
               "CD99 - CD99","CD99 - PILRA",
               "ANXA1 - FPR1","LAIR1 - LILRB4")

DAM.MTC.pathway<-c("APP - CD74",
                   "APP - (TREM2+TYROBP)",
                   "MDK - LRP1", 
                   "MDK - NCL",
                   "CD99 - CD99",
                   "CD99 - PILRA",
                   "ICAM1 - (ITGAM+ITGB2)",
                   "LAMA5 - CD44",
                   "COL1A1 - CD44",
                   "SPP1 - CD44",
                   "SPP1 - (ITGAV+ITGB1)",
                   "LGALS9 - P4HB",
                   "LGALS9 - CD44",
                   "TNF - TNFRSF1A",
                   "TNFSF12 - TNFRSF12A",
                   "GRN - SORT1")
cc.obj<-qs::qread('./dataset1/Lung/05cellchat/new/cc.obj_mtc.new.qs')
#pairLR <- extractEnrichedLR(cc.obj, signaling = c("ApoE","APP"), geneLR.return = FALSE)
targets.use = c('DAM','DAM-like.Mac','EMT_high',"REP_high",'EMT_low'); 
sources.use = c('DAM','DAM-like.Mac','APOE.Mac','EMT_high',"REP_high",'EMT_low')
df.net <- subsetCommunication(cc.obj, 
                              slot.name = "net",
                              sources.use = sources.use, 
                              targets.use = targets.use,
                              thresh = 0.05)
df.net<-subset(df.net,interaction_name_2%in%c(DAM.MTC.pathway,TAM.pathway))

df.net$source<-as.character(df.net$source)
df.net$target<-as.character(df.net$target)
df.net$source[df.net$source%in%c('EMT_high',"REP_high",'EMT_low')]<-'Tumor'
df.net$target[df.net$target%in%c('EMT_high',"REP_high",'EMT_low')]<-'Tumor'

df.net$source.target <- paste(df.net$source, 
                              df.net$target, sep = " -> ")
df.net$id<-paste(df.net$source.target,df.net$interaction_name_2)
df.net<-df.net %>% 
  arrange(id,desc(prob)) %>% 
  distinct(id,.keep_all = T)

df.net$pval[df.net$pval > 0.05] = 1
df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
df.net$pval[df.net$pval <= 0.01] = 3
df.net$prob[df.net$prob == 0] <- NA

table(df.net$source)
source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
df.net$target <- factor(df.net$target, levels = c('DAM','DAM-like.Mac','Tumor'))
df.net$source <- factor(df.net$source , levels = c('Tumor','DAM','DAM-like.Mac','APOE.Mac'))
df.net<- with(df.net, df.net[order(target, source),])
source.target.order <- unique(as.character(df.net$source.target))
df.net$source.target <- factor(df.net$source.target, 
                               levels = source.target.order)
df.net$dataset<-'Dataset1'
values <- c(1,2,3); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
df.net$interaction_name_2<-factor(df.net$interaction_name_2,
                                  levels=unique(c(DAM.MTC.pathway,TAM.pathway)))
ggplot(df.net[df.net$target=='Tumor',], 
       aes(x = source.target, 
           y = interaction_name_2, 
           color = prob, 
           size = pval)) +
  geom_point(pch = 16) +
  theme_bw() + 
  # theme(panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust= 1, 
                                   vjust = 1,),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(position = "bottom")+
  scale_radius(range = c(min(df.net[df.net$target=='Tumor',]$pval), max(df.net[df.net$target=='Tumor',]$pval)),
               breaks = sort(unique(df.net[df.net$target=='Tumor',]$pval)),
               labels = names(values)[values %in% sort(unique(df.net[df.net$target=='Tumor',]$pval))], 
               name = "p-value")+
  scale_fill_gradient2(low = '#0099CC',mid = "white",high = '#CC3333')+
  guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))+
  theme(text = element_text(size = 12),
        plot.title = element_text(size=12)) +
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 6))
qs::qsave(df.net,'./Lung/05cellchat/new/TAM_MTC_d1.qs')

##Cellchat between T cells and TAMs FigureS7======
rm(list=ls())
seu.T<-qs::qread('./Lung/05cellchat/02Tcluster/seu.T.qs')
table(seu.T$celltype)
genes<-intersect(rownames(seu.T@assays$RNA) ,rownames(seu.TAM@assays$RNA))

xcounts<-cbind(seu.T@assays$RNA@counts[genes,],
               seu.TAM@assays$RNA@counts[genes,])
table(seu.TAM$celltype_pl %>% is.na())
table(seu.T$celltype%>% is.na())
meta.data.T<-seu.T@meta.data[,c("sample","orig.ident","title","type","celltype")]
meta.data.TAM<-seu.TAM@meta.data[,c("sample","orig.ident","title","type","celltype_pl")]
colnames(meta.data.T)[5]<-c('cluster')
colnames(meta.data.TAM)[5]<-c('cluster')

meta.data<-rbind(meta.data.T,meta.data.TAM)
table(is.na(meta.data$cluster))
rownames(meta.data)
seu <- CreateSeuratObject(counts = xcounts,
                          meta.data = meta.data)
seu <- NormalizeData(seu)

data.input = seu[["RNA"]]@data # normalized data matrix
Idents(seu)<-seu$cluster
labels <- Idents(seu)
cc.obj <- createCellChat(object = seu, 
                         group.by = "ident", 
                         assay = "RNA")

cc.obj  <- addMeta(cc.obj , meta = meta.data)
cc.obj <- setIdent(cc.obj, ident.use = "cluster") 

levels(cc.obj@idents) 
groupSize <- as.numeric(table(cc.obj@idents))
CellChatDB <- CellChatDB.human
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB 

cc.obj@DB <- CellChatDB.use
cc.obj<- subsetData(cc.obj)
future::plan("multisession", workers = 4) # do parallel
cc.obj <- identifyOverExpressedGenes(cc.obj)
cc.obj <- identifyOverExpressedInteractions(cc.obj)
cc.obj <- computeCommunProb(cc.obj, type = "triMean")
cc.obj <- filterCommunication(cc.obj, min.cells = 10)
cc.obj<- computeCommunProbPathway(cc.obj)
cc.obj<- aggregateNet(cc.obj)
df.net <- subsetCommunication(cc.obj) 
qs::qsave(cc.obj,'./Lung/05cellchat/new/cc.obj_T.new.qs')
data.table::fwrite(df.net,
                   './Lung/05cellchat/new/df_net_DAM_T.new.txt',sep='\t',
                   quote = F,row.names = F)
df.net<-data.table::fread('./Lung/05cellchat/new/df_net_DAM_T.new.txt')
cc.obj <- netAnalysis_computeCentrality(cc.obj, slot.name = "netP") 
cc.obj <- computeNetSimilarity(cc.obj, type = "functional")
cc.obj <- netEmbedding(cc.obj, type = "functional")
cc.obj <- netClustering(cc.obj, type = "functional")
cc.obj <- computeNetSimilarity(cc.obj, type = "structural")
cc.obj <- netEmbedding(cc.obj, type = "structural")
cc.obj <- netClustering(cc.obj, type = "structural")



###chord plot T cell <-> TAMs----

#DAM_T_pathway<-readRDS('../DAM/DAM_T_pathway.rds')
DAM_T_pathway<-c('LGALS9 - CD45',
                 'LGALS9 - HAVCR2',
                 'SPP1 - CD44',
                 'CD99 - CD99')
cc.obj<-qs::qread('./Lung/05cellchat/new/cc.obj_T.new.qs')
df.net<-data.table::fread('./Lung/05cellchat/new/df_net_DAM_T.new.txt')

prob <- slot(cc.obj, "net")$prob
pval <- slot(cc.obj, "net")$pval
prob[pval > 0.05] <- 0
net <- reshape2::melt(prob, value.name = "prob")
colnames(net)[1:3] <- c("source","target","interaction_name")
cols.default <- c("interaction_name_2", "pathway_name",  "ligand",  "receptor" ,"annotation","evidence")
cols.common <- intersect(cols.default,colnames(cc.obj@LR$LRsig))
pairLR = dplyr::select(cc.obj@LR$LRsig, cols.common)
idx <- match(net$interaction_name, rownames(pairLR))
temp <- pairLR[idx,]
net <- cbind(net, temp)


net <- subset(net, interaction_name_2 %in% DAM_T_pathway)
sources.use = c('DAM','DAM-like.Mac'); 
#'DAM','DAM-like.Mac','APOE.Mac',
targets.use = c('CD4.T','CD8.T','Treg','DC')
table(cc.obj@idents)

net <- subset(net, source %in% sources.use)
net <- subset(net, target %in% targets.use)
net$source<-as.character(net$source)

df <- subset(net, prob > 0)
df$id<-paste(df$source,df$target,df$interaction_name,sep='.')
df<-df %>% 
  arrange(id,desc(prob)) %>% 
  distinct(id,.keep_all = T)
df$id <- 1:nrow(df)
df$receptor[df$receptor=='CD99']<-'CD99-R'

ligand.uni <- unique(df$ligand)
for (i in 1:length(ligand.uni)) {
  df.i <- df[df$ligand == ligand.uni[i], ]
  source.uni <- unique(df.i$source)
  for (j in 1:length(source.uni)) {
    df.i.j <- df.i[df.i$source == source.uni[j], ]
    df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
    df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
  }
}

receptor.uni <- unique(df$receptor)
for (i in 1:length(receptor.uni)) {
  df.i <- df[df$receptor == receptor.uni[i], ]
  target.uni <- unique(df.i$target)
  for (j in 1:length(target.uni)) {
    df.i.j <- df.i[df.i$target == target.uni[j], ]
    df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
    df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
  }
}

cell.order.sources <- levels(cc.obj@idents)[levels(cc.obj@idents) %in% sources.use]
cell.order.targets <- levels(cc.obj@idents)[levels(cc.obj@idents) %in% targets.use]

df$source <- factor(df$source, levels = cell.order.sources)
df$target <- factor(df$target, levels = cell.order.targets)

df.ordered.source <- df[with(df, order(source, -prob)), ]
df.ordered.target <- df[with(df, order(target, -prob)), ]

order.source <- unique(df.ordered.source[ ,c('ligand','source')])
order.target <- unique(df.ordered.target[ ,c('receptor','target')])
order.sector <- c(order.source$ligand, order.target$receptor)
color.TAM_hgc<-readRDS('../DAM/color.TAM_hgc.rds')

color.use <- c("CD4.T"="#56B4E9","CD8.T"="#B55D8D",
               "Treg"="#00BFC4",
               color.TAM_hgc[c('DAM','DAM-like.Mac','DC')])

# define edge color
edge.color <- color.use[as.character(df.ordered.source$source)]
names(edge.color) <- as.character(df.ordered.source$source)

grid.col.ligand <- color.use[as.character(order.source$source)]
names(grid.col.ligand) <- as.character(order.source$source)
grid.col.receptor <- color.use[as.character(order.target$target)]
names(grid.col.receptor) <- as.character(order.target$target)
grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
names(grid.col) <- order.sector
df.plot <- df.ordered.source[ ,c('ligand','receptor','prob')]
circos.clear()
library(circlize);
pdf('./Figuresubmit/TAM&Tcell.pdf',height=10,width=10)  
chordDiagram(df.plot,
             order = order.sector,
             col = edge.color,
             grid.col = grid.col,
             transparency = 0.4,
             link.border = NA,
             directional = 1,
             direction.type = c("diffHeight","arrows"),
             link.arr.type = "big.arrow",
             annotationTrack = "grid",
             annotationTrackHeight = 0.03,
             preAllocateTracks = list(track.height = max(strwidth(order.sector))),
             small.gap = 1,
             big.gap = 10,
             link.visible = T,
             scale = F,
             link.target.prop = T,
             reduce = -1)

circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = 0.8)
}, bg.border = NA)

lgd <- ComplexHeatmap::Legend(at = names(color.use), type = "grid", legend_gp = grid::gpar(fill = color.use), title = "Cell State")
ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(20, "mm"), y = unit(20, "mm"), just = c("right", "bottom"))

circos.clear()
gg <- recordPlot()
gg
dev.off()



###Prepare DIY dotplot----
rm(list=ls())
cc.obj<-qs::qread('./Lung/05cellchat/new/cc.obj_T.new.qs')

#pairLR <- extractEnrichedLR(cc.obj, signaling = c("ApoE","APP"), geneLR.return = FALSE)
sources.use = c('DAM','DAM-like.Mac'); 
targets.use = c('CD4.T','CD8.T','Treg','DC')
DAM_T_pathway<-readRDS('~/Brain/brian/DAM/DAM_T_pathway.rds')

df.net <- subsetCommunication(cc.obj, 
                              slot.name = "net",
                              sources.use = sources.use, 
                              targets.use = targets.use,
                              thresh = 0.05)
df.net<-subset(df.net,interaction_name_2%in%DAM_T_pathway)

df.net$target<-as.character(df.net$target)

df.net$source.target <- paste(df.net$source, 
                              df.net$target, sep = " -> ")
df.net$id<-paste(df.net$source.target,df.net$interaction_name_2)
df.net<-df.net %>% 
  arrange(id,desc(prob)) %>% 
  distinct(id,.keep_all = T)

df.net$pval[df.net$pval > 0.05] = 1
df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
df.net$pval[df.net$pval <= 0.01] = 3
df.net$prob[df.net$prob == 0] <- NA

table(df.net$target)
source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
df.net$target <- factor(df.net$target, levels = c('CD8.T','CD4.T','Treg','DC'))

df.net$source <- factor(df.net$source , levels = c('DAM','DAM-like.Mac'))
df.net<- with(df.net, df.net[order(target, source),])
source.target.order <- unique(as.character(df.net$source.target))
df.net$source.target <- factor(df.net$source.target, 
                               levels = source.target.order)
df.net$dataset<-'Dataset1'
values <- c(1,2,3); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
df.net$interaction_name_2<-factor(df.net$interaction_name_2,
                                  levels=DAM_T_pathway)
intersect(df.net$interaction_name_2,DAM_T_pathway)
ggplot(df.net, aes(x = source.target, y = interaction_name_2,
                   color = prob, size = pval)) +
  geom_point(pch = 16) +
  theme_bw() + 
  # theme(panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust= 1, 
                                   vjust = 1,),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(position = "bottom")+
  scale_radius(range = c(min(df.net$pval), max(df.net$pval)),
               breaks = sort(unique(df.net$pval)),
               labels = names(values)[values %in% sort(unique(df.net$pval))], 
               name = "p-value")+
  scale_fill_gradient2(low = '#0099CC',mid = "white",high = '#CC3333')+
  guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))+
  theme(text = element_text(size = 12),
        plot.title = element_text(size=12)) +
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 6))
qs::qsave(df.net,'./Lung/05cellchat/new/TAM_T_d1.qs')


