
color.use<-readRDS('./Atlas_ProjectTILs/panC.colSet.list.rds')

library(ComplexHeatmap)
##DAM_gene------
mat1<-readRDS('../dataset1/Figuresubmit/data.plot.wide.dataset1.rds')
mat2<-readRDS('../dataset2/Figuresubmit/data.plot.wide.dataset2.rds')
mat3<-readRDS('../dataset6/Figuresubmit/data.plot.wide.dataset6.rds')
mat4<-readRDS('../Lung_to_brain/Lung1_GSE123902/Figuresubmit/data.plot.wide.Lung1.rds')
mat5<-readRDS('../Lung_to_brain/Lung2_GSE202371/Figuresubmit/data.plot.wide.Lung2.rds')
mat6<-readRDS('../Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907//Figuresubmit/data.plot.wide.Lung3.rds')
mat_all<-cbind(mat1[,c('gene','DAM','APOE.Mac','DAM-like.Mac')],
               mat2[,c('DAM','APOE.Mac','DAM-like.Mac')],
               mat3[,c('DAM','APOE.Mac','DAM-like.Mac')],
               mat4[,c('DAM','APOE.Mac','DAM-like.Mac')],
               mat5[,c('DAM','APOE.Mac','DAM-like.Mac')],
               mat6[,c('DAM','APOE.Mac','DAM-like.Mac')]) %>% 
  as.data.frame() %>% 
  column_to_rownames('gene') %>% as.matrix()
dim(mat_all)
split = rep(1:6, each = 3)
col_dataset<-color.use$dataset %>% tail(6)
names(col_dataset)<-c('GSE234832','GSE186344','GSE143423','GSE123902','GSE202371','GSE131907')

col_dataset <- c("#1C7F93", "#D85FF7", "#FF84A0", "#E74C3C", "#8E44AD", "#3498DB")  # 自定义颜色
names(col_dataset) <- c('GSE234832', 'GSE186344', 'GSE143423', 'GSE123902', 'GSE202371', 'GSE131907')

# 构造热图注释
ha <- HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = col_dataset), 
                   labels = c('GSE234832', 'GSE186344', 'GSE143423', 'GSE123902', 'GSE202371', 'GSE131907')),
  col = list(foo = col_dataset)
)

library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), 
                     c("#0099CC", "white", '#CC3333'))
col_fun= colorRamp2(seq(min(mat_all), max(mat_all), length = 3), c("#0099CC", "white", '#CC3333'))
col_fun= colorRamp2(seq(min(mat_all), max(mat_all), length = 3), c("#0099CC", "white", '#CC3333'), space = "RGB")

dir.create('../DAM/FigureIntegrate/')
pdf('../DAM/FigureIntegrate/DAM_Gene.pdf',height=7,width=6)
Heatmap(mat_all, name = "mat", col = col_fun,
        column_split = split, 
        column_gap = unit(0.5, "cm"),  # 设置column之间的距离为2cm
       # row_gap = unit(3, "cm"),  
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_rows = F,cluster_columns = F,
        top_annotation = ha, 
        column_title = NULL)
dev.off()
##DAM score----
mat1<-readRDS('./dataset1/Figuresubmit/DAM.wide.dataset1_v2.rds')
mat2<-readRDS('./dataset2/Figuresubmit/DAM.wide.dataset2_v2.rds')
mat3<-readRDS('./dataset6/Figuresubmit/DAM.wide.dataset6_v2.rds')

mat4<-readRDS('./Lung_to_brain/Lung1_GSE123902/Figuresubmit/DAM.wide.Lung1_v2.rds')
mat5<-readRDS('./Lung_to_brain/Lung2_GSE202371//Figuresubmit/DAM.wide.Lung2_v2.rds')
mat6<-readRDS('./Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/Figuresubmit/DAM.wide.Lung3_v2.rds')
mat_all<-cbind(mat1[,c('DAM','CD14.Mono','S100A8.Mac','APOE.Mac','DAM-like.Mac',
                       'Prolif.Mac',
                       'Stressed.Mac',
                       'DC')],
               mat2[,c('DAM','MG','CD14.Mono','S100A8.Mac','APOE.Mac','DAM-like.Mac',
                       'Prolif.Mac','DC')],
               mat3[,c('DAM','APOE.Mac','DAM-like.Mac',
                       'Prolif.Mac',
                       'MT1H.Mac','Lung.derived.Mac','DC')],
               mat4[,c('DAM','CD14.Mono','S100A8.Mac',
                       'APOE.Mac','DAM-like.Mac',
                       'Prolif.Mac',
                       'MT1H.Mac',
                       'DC')],
               mat5[,c('DAM','CD14.Mono','S100A8.Mac',
                       'APOE.Mac','DAM-like.Mac',
                       'MT1H.Mac',
                       'DC',
                       'Neu')],
               mat6[,c('DAM','CD14.Mono','S100A8.Mac',
                       'APOE.Mac','DAM-like.Mac',
                       'Prolif.Mac',
                       'Stressed.Mac',
                       'MT1H.Mac',
                       'DC')]) %>% as.matrix()

split = c(rep(1:2, each=8),rep(3,7),rep(4:5,each=8),rep(6,9))




col_dataset <- c("#1C7F93", "#D85FF7", "#FF84A0", "#E74C3C", "#8E44AD", "#3498DB")  # 自定义颜色
names(col_dataset) <- c('GSE234832', 'GSE186344', 'GSE143423', 'GSE123902', 'GSE202371', 'GSE131907')
ha <- HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = col_dataset), 
                   labels = c('GSE234832', 'GSE186344', 'GSE143423', 'GSE123902', 'GSE202371', 'GSE131907')),
  col = list(foo = col_dataset)
)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), 
                     c("#0099CC", "white", '#CC3333'))
col_fun= colorRamp2(seq(min(mat_all), max(mat_all), length = 3), c("#0099CC", "white", '#CC3333'))
col_fun= colorRamp2(seq(min(mat_all), max(mat_all), length = 3), c("#0099CC", "white", '#CC3333'), space = "RGB")


pdf('./DAM/FigureIntegrate/DAM_Core_Programe.pdf',height=2.5,width=12)
Heatmap(mat_all, name = "mat", col = col_fun,
        column_split = split, 
        column_gap = unit(0.5, "cm"),  # 设置column之间的距离为2cm
        # row_gap = unit(3, "cm"),  
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_rows = F,cluster_columns = F,
        top_annotation = ha, 
        column_title = NULL)
dev.off()
##TAM.features-----
mat1<-readRDS('./dataset1/Figuresubmit/TAM.features.wide.dataset1.rds')
mat2<-readRDS('./dataset2/Figuresubmit/TAM.features.wide.dataset2.rds')
mat3<-readRDS('./dataset6/Figuresubmit/TAM.features.wide.dataset6.rds')

mat4<-readRDS('./Lung_to_brain/Lung1_GSE123902/Figuresubmit/TAM.features.wide.Lung1.rds')
mat5<-readRDS('./Lung_to_brain/Lung2_GSE202371//Figuresubmit/TAM.features.wide.Lung2.rds')
mat6<-readRDS('./Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/Figuresubmit/TAM.features.wide.Lung3.rds')
mat_all<-cbind(mat1[,c('DAM','CD14.Mono','S100A8.Mac','APOE.Mac','DAM-like.Mac',
                       'Prolif.Mac',
                       'Stressed.Mac',
                       'DC')],
               mat2[,c('DAM','MG','CD14.Mono','S100A8.Mac','APOE.Mac','DAM-like.Mac',
                       'Prolif.Mac','DC')],
               mat3[,c('DAM','APOE.Mac','DAM-like.Mac',
                       'Prolif.Mac',
                       'MT1H.Mac','Lung.derived.Mac','DC')],
               mat4[,c('DAM','CD14.Mono','S100A8.Mac',
                       'APOE.Mac','DAM-like.Mac',
                       'Prolif.Mac',
                       'MT1H.Mac',
                       'DC')],
               mat5[,c('DAM','CD14.Mono','S100A8.Mac',
                       'APOE.Mac','DAM-like.Mac',
                       'MT1H.Mac',
                       'DC',
                       'Neu')],
               mat6[,c('DAM','CD14.Mono','S100A8.Mac',
                       'APOE.Mac','DAM-like.Mac',
                       'Prolif.Mac',
                       'Stressed.Mac',
                       'MT1H.Mac',
                       'DC')]) %>% as.matrix()

split = c(rep(1:2, each=8),rep(3,7),rep(4:5,each=8),rep(6,9))




col_dataset <- c("#1C7F93", "#D85FF7", "#FF84A0", "#E74C3C", "#8E44AD", "#3498DB")  # 自定义颜色
names(col_dataset) <- c('GSE234832', 'GSE186344', 'GSE143423', 'GSE123902', 'GSE202371', 'GSE131907')
ha <- HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = col_dataset), 
                   labels = c('GSE234832', 'GSE186344', 'GSE143423', 'GSE123902', 'GSE202371', 'GSE131907')),
  col = list(foo = col_dataset)
)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), 
                     c("#0099CC", "white", '#CC3333'))

col_fun= colorRamp2(seq(min(mat_all), max(mat_all), length = 3), c("#0099CC", "white", '#CC3333'))
col_fun= colorRamp2(seq(min(mat_all), max(mat_all), length = 3), c("#0099CC", "white", '#CC3333'), space = "RGB")
rownames(mat_all)<-c("DAM-core-program","Stage1-DAM-upregulated","Stage2-DAM"   )
pdf('./DAM/FigureIntegrate/TAM.features2.pdf',height=2.5,width=12)
Heatmap(mat_all, name = "mat", col = col_fun,
        column_split = split, 
        column_gap = unit(0.5, "cm"),  # 设置column之间的距离为2cm
        # row_gap = unit(3, "cm"),  
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_rows = F,cluster_columns = F,
        top_annotation = ha, 
        column_title = NULL)
dev.off()
test <- system.file("extdata", "pbmc.markers.csv", package = "scRNAtoolVis") 
markers <- read.csv(test)

##cellchat MTC-----
TAM.pathway<-readRDS('./TAM.pathway.rds')
DAM.MTC.pathway<-readRDS('./DAM.MTC.pathway.rds')
DAM.MTC.pathway<-readRDS('./DAM.MTC.pathway.add.rds')
cc.mtc.d1<-qs::qread('../dataset1/Lung/05cellchat/new/TAM_MTC_d1.qs')
cc.mtc.d2<-qs::qread('../dataset2/Lung/05cellchat/new/TAM_MTC_d2.qs')
cc.mtc.d6<-qs::qread('../dataset6/Lung/05cellchat/new/TAM_MTC_d6.qs')
cc.mtc.l1<-qs::qread('../Lung_to_brain/Lung1_GSE123902/06cellchat/new/TAM_MTC_l1.qs')
cc.mtc.l2<-qs::qread('../Lung_to_brain/Lung2_GSE202371/03cellchat/new/TAM_MTC_l2.qs')
cc.mtc.l3<-qs::qread('../Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/03cellchat/new/TAM_MTC_l3.qs')
mat1.d1<-cc.mtc.d1[,c("interaction_name_2","source.target","prob")] %>% 
  pivot_wider(names_from = "source.target",values_from = 'prob') %>% 
  column_to_rownames('interaction_name_2')

mat1.d2<-cc.mtc.d2[,c("interaction_name_2","source.target","prob")] %>% 
  pivot_wider(names_from = "source.target",values_from = 'prob') %>% 
  column_to_rownames('interaction_name_2')
mat1.d6<-cc.mtc.d6[,c("interaction_name_2","source.target","prob")] %>% 
  pivot_wider(names_from = "source.target",values_from = 'prob') %>% 
  column_to_rownames('interaction_name_2')
mat1.l1<-cc.mtc.l1[,c("interaction_name_2","source.target","prob")] %>% 
  pivot_wider(names_from = "source.target",values_from = 'prob') %>% 
  column_to_rownames('interaction_name_2')
mat1.l2<-cc.mtc.l2[,c("interaction_name_2","source.target","prob")] %>% 
  pivot_wider(names_from = "source.target",values_from = 'prob') %>% 
  column_to_rownames('interaction_name_2')
mat1.l3<-cc.mtc.l3[,c("interaction_name_2","source.target","prob")] %>% 
  pivot_wider(names_from = "source.target",values_from = 'prob') %>% 
  column_to_rownames('interaction_name_2')
###mat mtc TAM Figure7A-----
colnames(mat1.d1)
mat1<-cbind(mat1.d1[DAM.MTC.pathway,c("Tumor -> DAM","Tumor -> DAM-like.Mac","DAM -> Tumor","DAM-like.Mac -> Tumor")],
            mat1.d2[DAM.MTC.pathway,c("Tumor -> DAM","Tumor -> DAM-like.Mac","DAM -> Tumor","DAM-like.Mac -> Tumor")],
            mat1.d6[DAM.MTC.pathway,c("Tumor -> DAM","Tumor -> DAM-like.Mac","DAM -> Tumor","DAM-like.Mac -> Tumor")],
            mat1.l1[DAM.MTC.pathway,c("Tumor -> DAM","Tumor -> DAM-like.Mac","DAM -> Tumor","DAM-like.Mac -> Tumor")],
            mat1.l2[DAM.MTC.pathway,c("Tumor -> DAM","Tumor -> DAM-like.Mac","DAM -> Tumor","DAM-like.Mac -> Tumor")],
            mat1.l3[DAM.MTC.pathway,c("Tumor -> DAM","Tumor -> DAM-like.Mac","DAM -> Tumor","DAM-like.Mac -> Tumor")]) %>% as.matrix()
mat1[is.na(mat1)]<-0
mat2<-apply(mat1,2,scale)
rownames(mat2)<-rownames(mat1)
colnames(mat2)[1:4]<-paste0('d1_',colnames(mat2)[1:4])
colnames(mat2)[5:8]<-paste0('d2_',colnames(mat2)[5:8])
colnames(mat2)[9:12]<-paste0('d6_',colnames(mat2)[9:12])
colnames(mat2)[13:16]<-paste0('l1_',colnames(mat2)[13:16])
colnames(mat2)[17:20]<-paste0('l2_',colnames(mat2)[17:20])
colnames(mat2)[21:24]<-paste0('l3_',colnames(mat2)[21:24])



mat2_long<-mat2 %>% as.data.frame() %>% 
  rownames_to_column('Pathway') %>% 
  pivot_longer(cols=2:24)
mat2_long$dataset<-sapply(strsplit(mat2_long$name,'_'),'[',1)
mat2_long$cell<-sapply(strsplit(mat2_long$name,'_'),'[',2)
mat2_long$Pathway<-factor(mat2_long$Pathway,levels = rev(DAM.MTC.pathway))
d<-c('GSE234832', 'GSE186344', 'GSE143423', 'GSE123902', 'GSE202371', 'GSE131907')
names(d)<-c('d1','d2','d6','l1','l2','l3')
mat2_long$dataset<-d[mat2_long$dataset]
mat2_long$dataset<-factor(mat2_long$dataset,levels=c('GSE234832', 'GSE186344', 'GSE143423', 'GSE123902', 'GSE202371', 'GSE131907'))
mat2_long$cell %>% table()
mat2_long$cell <-factor(mat2_long$cell ,levels=c("Tumor -> DAM","Tumor -> DAM-like.Mac","DAM -> Tumor","DAM-like.Mac -> Tumor"))
pdf('./cc.mtc.DAM.pdf',height=7,width=12)
ggplot(mat2_long,aes(x = cell, y = Pathway)) +
  geom_tile(aes(fill = value),
            color = "black") +
  theme_bw(base_size = 14) +
  xlab("") +ylab("") +
  coord_fixed(clip = "off") +
  theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
        panel.grid.major = element_blank(), # 去掉大背景网格线
        panel.grid.minor = element_blank(), # 去掉小背景网格线
        panel.background = element_blank(),
        axis.text = element_text(color = "black",size=12),
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0),
        legend.position = "right",
        strip.text = element_text(size = 12, angle = 0,
                                  face='bold'),  # 设置facet标题的字体大小和角度
        strip.background = element_blank(),                  # 去掉facet框的背景
        panel.border = element_blank())+
  guides(fill = guide_colorbar(title = "Probability",
                               title.position = "top",
                               title.hjust = 0.5,
                               barwidth = unit(0.7, "cm"),
                               direction = "vertical",
                               frame.colour = "black",
                               frame.linewidth = 0.5,
                               ticks.colour = "black"))+
  scale_fill_gradient2(low = '#0099CC',mid = "white",high = '#CC3333')+
  facet_wrap(~dataset,ncol = 6)
dev.off()


###mat TAM TAM FigureS7A-----
colnames(mat1.d1)
mat1<-cbind(mat1.d1[TAM.pathway,c('DAM -> DAM',
                                  'DAM-like.Mac -> DAM',
                                  'APOE.Mac -> DAM',
                                  'DAM -> DAM-like.Mac',
                                  'DAM-like.Mac -> DAM-like.Mac',
                                  'APOE.Mac -> DAM-like.Mac')],
            mat1.d2[TAM.pathway,c('DAM -> DAM',
                                  'DAM-like.Mac -> DAM',
                                  'APOE.Mac -> DAM',
                                  'DAM -> DAM-like.Mac',
                                  'DAM-like.Mac -> DAM-like.Mac',
                                  'APOE.Mac -> DAM-like.Mac')],
            mat1.d6[TAM.pathway,c('DAM -> DAM',
                                  'DAM-like.Mac -> DAM',
                                  'APOE.Mac -> DAM',
                                  'DAM -> DAM-like.Mac',
                                  'DAM-like.Mac -> DAM-like.Mac',
                                  'APOE.Mac -> DAM-like.Mac')],
            mat1.l1[TAM.pathway,c('DAM -> DAM',
                                  'DAM-like.Mac -> DAM',
                                  'APOE.Mac -> DAM',
                                  'DAM -> DAM-like.Mac',
                                  'DAM-like.Mac -> DAM-like.Mac',
                                  'APOE.Mac -> DAM-like.Mac')],
            mat1.l2[TAM.pathway,c('DAM -> DAM',
                                  'DAM-like.Mac -> DAM',
                                  'APOE.Mac -> DAM',
                                  'DAM -> DAM-like.Mac',
                                  'DAM-like.Mac -> DAM-like.Mac',
                                  'APOE.Mac -> DAM-like.Mac')],
            mat1.l3[TAM.pathway,c('DAM -> DAM',
                                  'DAM-like.Mac -> DAM',
                                  'APOE.Mac -> DAM',
                                  'DAM -> DAM-like.Mac',
                                  'DAM-like.Mac -> DAM-like.Mac',
                                  'APOE.Mac -> DAM-like.Mac')]) %>% as.matrix()
mat1[is.na(mat1)]<-0
mat2<-apply(mat1,2,scale)
rownames(mat2)<-rownames(mat1)
colnames(mat2)[1:6]<-paste0('d1_',colnames(mat2)[1:6])
colnames(mat2)[7:12]<-paste0('d2_',colnames(mat2)[7:12])
colnames(mat2)[13:18]<-paste0('d6_',colnames(mat2)[13:18])
colnames(mat2)[19:24]<-paste0('l1_',colnames(mat2)[19:24])
colnames(mat2)[25:30]<-paste0('l2_',colnames(mat2)[25:30])
colnames(mat2)[31:36]<-paste0('l3_',colnames(mat2)[31:36])



mat2_long<-mat2 %>% as.data.frame() %>% 
  rownames_to_column('Pathway') %>% 
  pivot_longer(cols=2:37)
mat2_long$dataset<-sapply(strsplit(mat2_long$name,'_'),'[',1)
mat2_long$cell<-sapply(strsplit(mat2_long$name,'_'),'[',2)
mat2_long$Pathway<-factor(mat2_long$Pathway,levels = rev(TAM.pathway))
d<-c('GSE234832', 'GSE186344', 'GSE143423', 'GSE123902', 'GSE202371', 'GSE131907')
names(d)<-c('d1','d2','d6','l1','l2','l3')
mat2_long$dataset<-d[mat2_long$dataset]
table(mat2_long$cell)
mat2_long$cell<-factor(mat2_long$cell,
                       levels=c('DAM -> DAM',
                                'DAM-like.Mac -> DAM',
                                'APOE.Mac -> DAM',
                                'DAM -> DAM-like.Mac',
                                'DAM-like.Mac -> DAM-like.Mac',
                                'APOE.Mac -> DAM-like.Mac'))

mat2_long$dataset<-factor(mat2_long$dataset,levels=c('GSE234832', 'GSE186344', 'GSE143423', 'GSE123902', 'GSE202371', 'GSE131907'))


pdf('./cc.TAM.DAM.pdf',height=12,width=14)
ggplot(mat2_long,aes(x = cell, y = Pathway)) +
  geom_tile(aes(fill = value),
            color = "black") +
  theme_bw(base_size = 14) +
  xlab("") +ylab("") +
  coord_fixed(clip = "off") +
  theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
        axis.text = element_text(color = "black",size=12),
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0),
        legend.position = "right",
        strip.text = element_text(size = 12, angle = 0,
                                  face='bold'),  # 设置facet标题的字体大小和角度
        strip.background = element_blank(),                  # 去掉facet框的背景
        panel.border = element_blank())+
  guides(fill = guide_colorbar(title = "Probability",
                               title.position = "top",
                               title.hjust = 0.5,
                               barwidth = unit(0.7, "cm"),
                               direction = "vertical",
                               frame.colour = "black",
                               frame.linewidth = 0.5,
                               ticks.colour = "black"))+
  scale_fill_gradient2(low = '#0099CC',mid = "white",high = '#CC3333')+
  facet_wrap(~dataset,ncol = 6)
dev.off()


###mat TAM T FigureS7B------
DAM_T_pathway<-readRDS('~/Brain/brian/DAM/DAM_T_pathway.rds')
cc.T.d1<-qs::qread('../dataset1/Lung/05cellchat/new/TAM_T_d1.qs')
cc.T.d2<-qs::qread('../dataset2/Lung/05cellchat/new/TAM_T_d2.qs')
cc.T.d6<-qs::qread('../dataset6/Lung/05cellchat/new/TAM_T_d6.qs')
cc.T.l1<-qs::qread('../Lung_to_brain/Lung1_GSE123902/06cellchat/new/TAM_T_l1.qs')
cc.T.l2<-qs::qread('../Lung_to_brain/Lung2_GSE202371/03cellchat/new/TAM_T_l2.qs')
cc.T.l3<-qs::qread('../Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/03cellchat/new/TAM_T_l3.qs')
mat1.d1<-cc.T.d1[,c("interaction_name_2","source.target","prob")] %>% 
  pivot_wider(names_from = "source.target",values_from = 'prob') %>% 
  column_to_rownames('interaction_name_2')

mat1.d2<-cc.T.d2[,c("interaction_name_2","source.target","prob")] %>% 
  pivot_wider(names_from = "source.target",values_from = 'prob') %>% 
  column_to_rownames('interaction_name_2')
mat1.d6<-cc.T.d6[,c("interaction_name_2","source.target","prob")] %>% 
  pivot_wider(names_from = "source.target",values_from = 'prob') %>% 
  column_to_rownames('interaction_name_2')

mat1.d6<-rbind(mat1.d6,rep(NA,6))
rownames(mat1.d6)[8]<-setdiff(DAM_T_pathway,rownames(mat1.d6))

mat1.l1<-cc.T.l1[,c("interaction_name_2","source.target","prob")] %>% 
  pivot_wider(names_from = "source.target",values_from = 'prob') %>% 
  column_to_rownames('interaction_name_2')
mat1.l1<-rbind(mat1.l1,rep(NA,8),rep(NA,8),rep(NA,8),rep(NA,8))
rownames(mat1.l1)[5:8]<-setdiff(DAM_T_pathway,rownames(mat1.l1))
mat1.l2<-cc.T.l2[,c("interaction_name_2","source.target","prob")] %>% 
  pivot_wider(names_from = "source.target",values_from = 'prob') %>% 
  column_to_rownames('interaction_name_2')
mat1.l3<-cc.T.l3[,c("interaction_name_2","source.target","prob")] %>% 
  pivot_wider(names_from = "source.target",values_from = 'prob') %>% 
  column_to_rownames('interaction_name_2')

colnames(mat1.d1)
mat1<-cbind(mat1.d1[DAM_T_pathway,],
            mat1.d2[DAM_T_pathway,],
            mat1.d6[DAM_T_pathway,],
            mat1.l1[DAM_T_pathway,],
            mat1.l2[DAM_T_pathway,],
            mat1.l3[DAM_T_pathway,]) %>% as.matrix()
mat1[is.na(mat1)]<-0

mat2<-apply(mat1,2,scale)
rownames(mat2)<-rownames(mat1)
colnames(mat2)[1:8]<-paste0('d1_',colnames(mat2)[1:8])
colnames(mat2)[9:16]<-paste0('d2_',colnames(mat2)[9:16])
colnames(mat2)[17:22]<-paste0('d6_',colnames(mat2)[17:22])
colnames(mat2)[23:30]<-paste0('l1_',colnames(mat2)[23:30])
colnames(mat2)[31:36]<-paste0('l2_',colnames(mat2)[31:36])
colnames(mat2)[37:44]<-paste0('l3_',colnames(mat2)[37:44])



mat2_long<-mat2 %>% as.data.frame() %>% 
  rownames_to_column('Pathway') %>% 
  pivot_longer(cols=2:45)
mat2_long$dataset<-sapply(strsplit(mat2_long$name,'_'),'[',1)
mat2_long$cell<-sapply(strsplit(mat2_long$name,'_'),'[',2)
mat2_long$Pathway<-factor(mat2_long$Pathway,
                          levels = rev(c('LGALS9 - CD45',
                                         'LGALS9 - HAVCR2',
                                         'SPP1 - CD44',
                                         'CD99 - CD99',
                                         'GAS6 - AXL',
                                         'TNF - TNFRSF1B',
                                         'PGE2-PTGES3 - PTGER4',
                                         'ADORA3 - ENTPD1')))
d<-c('GSE234832', 'GSE186344', 'GSE143423', 'GSE123902', 'GSE202371', 'GSE131907')
names(d)<-c('d1','d2','d6','l1','l2','l3')
mat2_long$dataset<-d[mat2_long$dataset]
pdf('./cc.DAM.T.pdf',height=14,width=14)
ggplot(mat2_long,aes(x = cell, y = Pathway)) +
  geom_tile(aes(fill = value),
            color = "black") +
  theme_bw(base_size = 14) +
  xlab("") +ylab("") +
  coord_fixed(clip = "off") +
  theme(plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
        panel.grid.major = element_blank(), # 去掉大背景网格线
        panel.grid.minor = element_blank(), # 去掉小背景网格线
        panel.background = element_blank(),  # 去掉背景填充
        #axis.line = element_line(colour = "black") ,
        axis.text = element_text(color = "black",size=12),
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0),
        legend.position = "right", Figure
        strip.text = element_text(size = 12, angle = 0,
                                  face='bold'),  # 设置facet标题的字体大小和角度
        strip.background = element_blank(),                  # 去掉facet框的背景
        panel.border = element_blank())+
  guides(fill = guide_colorbar(title = "Probability",
                               title.position = "top",
                               title.hjust = 0.5,
                               barwidth = unit(0.7, "cm"),
                               direction = "vertical",
                               frame.colour = "black",
                               frame.linewidth = 0.5,
                               ticks.colour = "black"))+
  scale_fill_gradient2(low = '#0099CC',mid = "white",high = '#CC3333')+
  facet_wrap(~dataset,ncol = 6)
dev.off()
