rm(list=ls())
library(GEOquery)
library(tidyverse)
library(arrayQualityMetrics)
library(affy)
library(limma)
DAM_features<-c('APOE','TREM2','TYROBP',
                'SPP1','OLR1','CCL3','CCL4',
                'CD83','CSF1R','MS4A7','FCGR1A',
                'LILRB4', 'TMIGD3', 'GPR34','VSIG4',
                'FCGR2A',
                'MS4A4A','MS4A6A','CCL18',
                'CD163','MRC1','DAB2',
                'CTSD','CTSB','CTSL','CTSS',
                'PSAP','FTL',"GRN",'CAPG')
#GSE14690------
data<-getGEO('GSE14690',file ='~/Brain/download/GSE14690/GSE14690_series_matrix.txt.gz',getGPL = F )
data@protocolData
expSet<-exprs(data)
pdata<-pData(data)
gpl<-data.table::fread('./GPL8128.txt')
gpl$ID<-as.character(gpl$ID)


colnames(data)
pdata<-pdata[,c("title","geo_accession","source_name_ch1")]
pdata$set<-sapply(strsplit(pdata$title,': '),'[',1)           
pdata$set<-sapply(strsplit(pdata$set,', '),'[',2)           
pdata$title<-sapply(strsplit(pdata$title,': '),'[',2)    
colnames(pdata)[3]<-'group'
table(pdata$group)
pdata<-pdata %>% 
  filter(!(str_detect(group,'mice')))
table(pdata$group)
pdata<-pdata %>% 
  filter(str_detect(group,'breast cancer') | group=='brain metastasis')
pdata$group<-ifelse(pdata$group=='primary breast cancer','primary','Brmt')
table(pdata$set,pdata$group)
group<-as.factor(pdata1$group)
design <- model.matrix(~group)
colnames(design)<-levels(group)
fit<-lmFit(edata1,design)
contrast.matrix<-makeContrasts(Brmt-primary,levels=design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
topTable(fit2,coef=1,adjust='BH')
results <- topTable(fit2,adjust='fdr',coef=1,number=Inf)
results <- results %>% 
  rownames_to_column("ID") %>%
  inner_join(gpl[,c('ID','SYMBOL')],by="ID") %>%
  dplyr::select(-ID) %>%
  dplyr::select(SYMBOL,everything()) %>%
  arrange(desc(logFC)) %>% 
  filter(P.Value<0.05&logFC>=1) %>% 
  distinct(SYMBOL,.keep_all = T)



##绘图火山图 Figure8C-----
pdf('./GSE14690/result_set1_volcano.pdf',height=6,width=6)
ggplot(result_set1, 
      aes(x=logFC, 
          y=-log10(adj.P.Val),
          fill=logFC)) +
  geom_jitter(size=2,shape=21,color='gray') +
  scale_fill_gradient2(low = '#0099CC',
                       mid = "white",
                       high = '#CC3333')+
  ggrepel::geom_text_repel(inherit.aes = F,
                           data = subset(result_set1,
                                         SYMBOL%in%DAM_features),
                           mapping = aes(x=logFC, 
                                         y=-log10(adj.P.Val),
                                         label = SYMBOL),
                           nudge_x = .1, nudge_y = .01,
                           max.overlaps = Inf,
                           color='blue',
                           size = 6) +
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
