rm(list=ls())
library(GEOquery)
library(tidyverse)
library(arrayQualityMetrics)
library(affy)
library(limma)
# GSE125989---------
# Expression profiling by array
# Matching primary breast cancers and brain metastases from 16 patients.
# GPL571 	[HG-U133A_2] Affymetrix Human Genome U133A 2.0 Array

file<-list.files('~/Brain/download/GSE125989/GSE125989_RAW/')
cel_files <- affy::list.celfiles("~/Brain/download/GSE125989/GSE125989_RAW", 
                                 full.name = F)
head(cel_files)
data<-getGEO('GSE125989',file ='~/Brain/download/GSE125989/GSE125989_series_matrix.txt.gz',getGPL = F )
data<-pData(data)
data<-data[,c("title","geo_accession",
              "source_name_ch1","er-by-gene:ch1",
              "her2-by-gene:ch1","characteristics_ch1.5")]
data[1,]
colnames(data)<-c('id','geo_accession','histology','ER','HER2','pair')
data$id<-sapply(strsplit(data$id,': '),'[',1)

data$ER<-str_extract(data$ER,'negative|positive')
data$HER2<-str_extract(data$HER2,'negative|positive')
data$pair<-sapply(strsplit(data$pair,': '),'[',2)
table(data$ER)
table(data$HER2)
table(data$id)
table(data$pair)
cel_files <- data.frame(filename=cel_files[str_match(cel_files,"GSM\\d+")%in%data$geo_accession],
                        geo_accession=str_extract(cel_files,'GSM\\d+'))
data<-data%>%
  inner_join(cel_files)
table(data$pair)
data$histology<-ifelse(str_detect(data$histology,'brain'),'BrM','BC')
data$type<-paste0(data$ER,'_',data$HER2)
table(data$type)
data$type[data$type=='negative_negative']<-'TNBC'
data$type[data$type=='negative_positive']<-'HER2'
data$type[data$type=='positive_negative']<-'LuminalA'
data$type[data$type=='positive_positive']<-'LuminalB'
data<-data%>% arrange(pair)
target<-data
data<-lapply(1:16,function(x){
  p1<-data[data$pair==x,]
  p1$consent<-ifelse(p1$type[p1$histology=='BC']==p1$type[p1$histology=='BrM'],
                     'T','F')
  return(p1)  
  
}) %>% do.call(rbind,.)
#data<-data %>% filter(consent=='T')
rownames(data)<-data$filename
cel <- ReadAffy(filenames = data$filename,
                celfile.path = "/Users/huzixin/Brain/download/GSE125989/GSE125989_RAW",
                phenoData = data)
boxplot(cel, las = 2)
# eset <- expresso(cel, normalize.method="qspline",
#                  bgcorrect.method="rma",pmcorrect.method="pmonly",
#                  summary.method="liwong")
eset<-rma(cel)
boxplot(eset, las = 2)
exprSet<-exprs(eset)
target<-pData(eset)
colnames(exprSet)<-eset$geo_accession
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  ## 取log2
  exprSet <- log2(ex)
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}

boxplot(eset,outline=FALSE,notch=T, las=2)

library(hgu133a2.db)

probeid<-as.data.frame(toTable(hgu133a2SYMBOL))
exprSet2 <- as.data.frame(exprSet) %>%
  rownames_to_column("probe_id") %>%
  inner_join(probeid,by="probe_id") %>%
  dplyr::select(-probe_id) %>%
  dplyr::select(symbol,everything()) %>%
  mutate(rowMean =rowMeans(.[,-1])) %>%
  arrange(desc(rowMean)) %>%
  distinct(symbol,.keep_all = T) %>%
  dplyr::select(-rowMean) %>%
  column_to_rownames("symbol")
#差异分析------
group<-as.factor(target$histology)
#Levels: BC BrM
#group <- factor(group,levels=c())
design<-model.matrix(~1+group)
colnames(design)<-levels(group)
fit<-lmFit(eset,design)
contrast.matrix<-makeContrasts(BrM-BC,levels=design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
topTable(fit2,coef=1,adjust='BH')
results <- topTable(fit2,adjust='fdr',coef=1,number=Inf)
group<-as.factor(data$histology)
block<-factor(data$pair)

design<-model.matrix(~0+group+block)
colnames(design)[1:2]<-levels(group)
fit<-lmFit(eset,design)
contrast.matrix<-makeContrasts(BrM-BC,levels=design)
fit<-contrasts.fit(fit,contrast.matrix)


fit<-eBayes(fit)
results <- topTable(fit,adjust='fdr',coef=1,number=Inf)

results <- results %>% 
  rownames_to_column("probe_id") %>%
  inner_join(probeid,by="probe_id") %>%
  dplyr::select(-probe_id) %>%
  dplyr::select(symbol,everything()) %>%
  arrange(symbol,desc(logFC)) %>% 
  distinct(symbol,.keep_all = T) %>% 
  filter(P.Value<0.05)

results_GSE125989<-results
data.table::fwrite(results_GSE125989,file ='./results_GSE125989.txt',
                   sep='\t',quote = F,row.names = F )
data.table::fwrite(results_GSE125989,file ='./Transcriptome/GSE125989/results_GSE125989.txt',
                   sep='\t',quote = F,row.names = F )


qs::qsave(eset,'./GSE125989/GSE125989_eset_rma.qs')
data.table::fwrite(target,'./GSE125989/GSE125989_target_all.txt',
                   sep='\t',quote = F,row.names = F )
data.table::fwrite(data,'./GSE125989/GSE125989_target.txt',
                   sep='\t',quote = F,row.names = F )



##绘图火山图-----
DAM_features<-c('APOE','TREM2','TYROBP',
                'SPP1','OLR1','CCL3','CCL4',
                'CD83','CSF1R','MS4A7','FCGR1A',
                'LILRB4', 'TMIGD3', 'GPR34','VSIG4',
                'FCGR2A',
                'MS4A4A','MS4A6A','CCL18',
                'CD163','MRC1','DAB2',
                'CTSD','CTSB','CTSL','CTSS',
                'PSAP','FTL',"GRN",'CAPG')
eset<-qs::qread('./GSE125989_eset_rma.qs')
results<-data.table::fread(file ='./results_GSE125989.txt') %>% 
  filter(P.Value<0.05&logFC>=0)
data<-data.table::fread('./GSE125989_target.txt',)
target<-data.table::fread('./GSE125989_target_all.txt')
table(data$histology)

pdf('./result_GSE125989_volcano.pdf',height=6,width=6)
ggplot(results, 
       aes(x=logFC, 
           y=-log10(P.Value),
           fill=logFC)) +
  geom_jitter(size=2,shape=21,color='gray') +
  scale_fill_gradient2(low = '#0099CC',
                       mid = "white",
                       high = '#CC3333')+
  ggrepel::geom_text_repel(inherit.aes = F,
                           data = subset(results,
                                         symbol%in%DAM_features),
                           mapping = aes(x=logFC, 
                                         y=-log10(P.Value),
                                         label = symbol),
                           nudge_x = .1, nudge_y = .01,
                           max.overlaps = Inf,
                           color='black',
                           size = 6) +
  geom_jitter(data = subset(results,
                            symbol%in%DAM_features),
              aes(x=logFC, 
                  y=-log10(P.Value),
                  fill=logFC),
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

dev.off()Figure8C

