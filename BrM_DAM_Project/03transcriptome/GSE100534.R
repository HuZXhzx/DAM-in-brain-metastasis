#GSE100534-------
rm(list=ls())
library(GEOquery)
library(tidyverse)
library(arrayQualityMetrics)
library(affy)
library(limma)
data<-getGEO('GSE100534',file ='~/Brain/download/GSE100534/GSE100534_series_matrix.txt.gz',getGPL = F )
data<-pData(data)
DAM_features<-c('APOE','TREM2','TYROBP',
                'SPP1','OLR1','CCL3','CCL4',
                'CD83','CSF1R','MS4A7','FCGR1A',
                'LILRB4', 'TMIGD3', 'GPR34','VSIG4',
                'FCGR2A',
                'MS4A4A','MS4A6A','CCL18',
                'CD163','MRC1','DAB2',
                'CTSD','CTSB','CTSL','CTSS',
                'PSAP','FTL',"GRN",'CAPG')
#GPL6244 	
#[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]
colnames(data)
data<-data[,c(1:2,38,40:44)]
colnames(data)<-c('title','geo_accession','relation','age','gender','grade','Her2','Histology')

data$title<-sapply(strsplit(data$title,': '),'[',1)
data$relation<-sapply(strsplit(data$relation,': '),'[',2)
data<-data %>% filter(!(title=='meningioma sample'))

colnames(data_add)
data_add<-getGEO('GSE36295',file ='~/Brain/download/GSE100534/GSE36295_series_matrix.txt.gz',getGPL = F )
data_add<-pData(data_add)
data_add<-data_add[,c(1:2,38:44,46)] %>% 
  filter(!(str_detect(title,'Normal')))
data_add[1,]
colnames(data_add)<-c('title_add','relation','age_add','ER','grade_add','HER2_add','Histology_add','PR','Site','TNBC')

data<-data %>% left_join(data_add)
data$age<-str_extract(data$age,'\\d+')


data$Her2<-sapply(strsplit(data$Her2,' '),'[',2)

data$grade<-sapply(strsplit(data$grade,'grade '),'[',2)
data$Histology[1:16]<-data$Histology_add[1:16]
data$grade[1:16]<-data$grade_add[1:16]

colnames(data)
data<-data[,-c(9,10,12,14,16)]
data$TNBC<-as.numeric(data$TNBC)
data$ER<-as.numeric(data$ER)
data$PR<-as.numeric(data$PR)
data$HER2_add<-as.numeric(data$HER2_add)

file<-list.files('~/Brain/download/GSE100534/GSE100534_RAW/')
cel_files <- affy::list.celfiles("~/Brain/download/GSE100534/GSE100534_RAW/", 
                                 full.name = F)

cel_files <- cel_files [str_extract(cel_files,'GSM\\d+')%in%data$geo_accession]
cel_files <- data.frame(filename=cel_files,
                        geo_accession=str_extract(cel_files,'GSM\\d+'))

data<-data %>% inner_join(cel_files)
data$title[str_detect(data$title,'BC sample')]<-'BC'
data$title[str_detect(data$title,'BCBM sample')]<-'BrM'
data2<-data %>% filter((!(is.na(Her2))) & Histology=='IDC' & (!(Her2==0)))
data2<-data %>% 
  filter(geo_accession%in%c('GSM2686112','GSM2686120',
                            'GSM2686113','GSM2686119',
                            'GSM2686102','GSM2686118')) 

rownames(data2)<-data2$geo_accession
data2<-data2[c('GSM2686112','GSM2686120',
               'GSM2686113','GSM2686119',
               'GSM2686102','GSM2686118'),]
data2$pair<-rep(1:3,each=2)




cel <- ReadAffy(filenames = data2$filename,
                celfile.path = "/Users/huzixin/Brain/download/GSE100534/GSE100534_RAW/",
                phenoData = data2)



library(hugene10stv1cdf)
library(hugene10sttranscriptcluster.db)

eset<-rma(cel)
boxplot(eset,outline=FALSE,notch=T, las=2)

eset$geo_accession
probeid<-as.data.frame(toTable(hugene10sttranscriptclusterSYMBOL))
group<-as.factor(data2$title)
block<-as.factor(data2$pair)
design<-model.matrix(~0+group+block)
colnames(design)[1:2]<-levels(group)
fit<-lmFit(eset,design)
contrast.matrix<-makeContrasts(BrM-BC,levels=design)
fit<-contrasts.fit(fit,contrast.matrix)

fit<-eBayes(fit)
#results <- topTable(fit,adjust='fdr',coef="groupBrM",number=Inf)
results <- topTable(fit,adjust='fdr',coef=1,number=Inf)
results <- results %>% 
  rownames_to_column("probe_id") %>%
  inner_join(probeid,by="probe_id") %>%
  dplyr::select(-probe_id) %>%
  dplyr::select(symbol,everything()) %>%
  arrange(symbol,desc(logFC)) %>% 
  distinct(symbol,.keep_all = T) %>% 
  filter(P.Value<0.05)

dir.create('./Transcriptome/GSE100534')
data.table::fwrite(results,file ='./results_GSE100534.txt',
                   sep='\t',quote = F,row.names = F )
results <- results %>% 
  filter(logFC>=0)

qs::qsave(eset,'./GSE100534/GSE100534_eset_rma.qs')
data.table::fwrite(data2,'./Transcriptome/GSE100534/GSE100534_target.txt',
                   sep='\t',quote = F,row.names = F )
data.table::fwrite(data,'./Transcriptome/GSE100534/GSE100534_target_all.txt',
                   sep='\t',quote = F,row.names = F )

results<-data.table::fread('./results_GSE100534.txt')


##绘图火山图Figure8C-----
pdf('./result_GSE100534_volcano.pdf',height=6,width=6)
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

dev.off()





