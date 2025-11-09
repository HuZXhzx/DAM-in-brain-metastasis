rm(list=ls())
library(GEOquery)
library(tidyverse)
library(arrayQualityMetrics)
library(affy)
library(limma)
##compare-----
data<-getGEO('GSE248830',file ='~/Brain/download/GSE248830/GSE248830_series_matrix.txt.gz',getGPL = F )
data@protocolData
pdata<-pData(data)
colnames(pdata)
pdata<-pdata[,c(1,2,8,33:37)]
colnames(pdata)<-c('title','geo_accession','resource','age','histology','gender','smoking','treat')
expSet<-data.table::fread('~/Brain/download/GSE248830/GSE248830_Raw_data.csv.gz') %>% 
  column_to_rownames('Probe Name')

boxplot(expSet,outline=FALSE,notch=T, las=2)
expSet=normalizeBetweenArrays(expSet)
boxplot(expSet,outline=FALSE, notch=T, las=2)
expSet <- as.data.frame(expSet)
data_colnames<-data.frame(id=colnames(expSet))
data_colnames$title[str_detect(data_colnames$id,'Primary breast')]<-paste0('primary_breast_cancer_',str_extract(data_colnames$id[str_detect(data_colnames$id,'Primary breast')],'\\d+'))
data_colnames$title[str_detect(data_colnames$id,'BCBM')]<-paste0('breast_cancer_matched_paired_brain_metastasis_',str_extract(data_colnames$id[str_detect(data_colnames$id,'BCBM')],'\\d+'))
data_colnames$title[str_detect(data_colnames$id,'Primary LUAD')]<-paste0('primary_LUAD_',str_extract(data_colnames$id[str_detect(data_colnames$id,'Primary LUAD')],'\\d+'))
data_colnames$title[str_detect(data_colnames$id,'BM-LUAD')]<-paste0('lung_ adenocarcinoma_matched_paired_brain_metastasis_',str_extract(data_colnames$id[str_detect(data_colnames$id,'BM-LUAD')],'\\d+'))
d<-pdata$geo_accession
names(d)<-pdata$title
data_colnames$geo_accession<-d[data_colnames$title]
colnames(expSet)<-data_colnames$geo_accession
expSet <-expSet[,pdata$geo_accession]
table(pdata$resource)
d<-c('breast_matched_paired_brain_metastasis'='BrBM',
     'lung_matched_paired_brain_metastasis'='LuBM',
     'primary breast'='BC',
     'primary LUAD'='LUAD')
pdata$resource<-d[pdata$resource]
pdata$orig<-ifelse(pdata$resource%in%c('LuBM','LUAD'),'Lung','Breast')
table(pdata$orig)
pdata$pair<-paste(pdata$orig,str_extract(pdata$title,'\\d+'),sep='.')
pdata$group<-ifelse(pdata$resource%in%c('BC','LUAD'),'primary','brain')


ex <-expSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

## 开始判断
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  ## 取log2
  exprSet <- log2(ex)
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}
group<-factor(pdata$group,levels=c('brain','primary'))
block<-as.factor(pdata$pair)
design<-model.matrix(~0+group+block)
colnames(design)[1:2]<-levels(group)
contrast.matrix<-makeContrasts(brain-primary,levels=design)
fit<-lmFit(exprSet,design)

fit<-contrasts.fit(fit,contrast.matrix)

fit<-eBayes(fit)
#results <- topTable(fit,adjust='fdr',coef="groupBrM",number=Inf)
results <- topTable(fit,adjust='fdr',coef=1,number=Inf)
results <- results %>% 
  rownames_to_column("symbol") %>%
  dplyr::select(symbol,everything()) %>%
  arrange(desc(logFC)) %>% 
  #distinct(symbol,.keep_all = T) %>% 
  filter(P.Value<0.05)

data.table::fwrite(results,file ='./results_GSE248830.txt',
                   sep='\t',quote = F,row.names = F )
results <- results %>% 
  filter(logFC>=0)
##绘图火山图Figure8C-----
DAM_features<-c('APOE','TREM2','TYROBP',
                'SPP1','OLR1','CCL3','CCL4',
                'CD83','CSF1R','MS4A7','FCGR1A',
                'LILRB4', 'TMIGD3', 'GPR34','VSIG4',
                'FCGR2A',
                'MS4A4A','MS4A6A','CCL18',
                'CD163','MRC1','DAB2',
                'CTSD','CTSB','CTSL','CTSS',
                'PSAP','FTL',"GRN",'CAPG')

pdf('./result_GSE248830_volcano.pdf',height=6,width=6)
ggplot(results, 
       aes(x=logFC, 
           y=-log10(adj.P.Val),
           fill=logFC)) +
  geom_jitter(size=2,shape=21,color='gray') +
  scale_fill_gradient2(low = '#0099CC',
                       mid = "white",
                       high = '#CC3333')+
  ggrepel::geom_text_repel(inherit.aes = F,
                           data = subset(results,
                                         symbol%in%DAM_features),
                           mapping = aes(x=logFC, 
                                         y=-log10(adj.P.Val),
                                         label = symbol),
                           nudge_x = .1, nudge_y = .01,
                           max.overlaps = Inf,
                           color='black',
                           size = 6) +
  
  geom_jitter(data = subset(results,
                            symbol%in%DAM_features),
              aes(x=logFC, y=-log10(adj.P.Val),
                  fill=logFC),
              size=4,shape=21,color='black')+
  
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

