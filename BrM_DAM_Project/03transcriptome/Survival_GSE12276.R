rm(list=ls())
library(GEOquery)
library(tidyverse)
library(arrayQualityMetrics)
library(affy)
library(limma)


#GSE12276--------
#GPL570
data1<-getGEO('GSE12276',file ='~/Brain/download/GSE12276/GSE12276_series_matrix.txt',getGPL = F )
data1<-pData(data1)

colnames(data1)
data1<-data1[,c(1,2,10,11,19)]
data1[1,]
colnames(data1)<-c('Primary_breast_cancer','geo_accession',
                   'Relapse','Survival_time','Description')
data1$Primary_breast_cancer<-sapply(strsplit(data1$Primary_breast_cancer,', '),'[',2)
data1$Primary_breast_cancer<-str_replace_all(data1$Primary_breast_cancer,' ','_')
data1$Relapse<-sapply(strsplit(data1$Relapse,': '),'[',2)
data1$Survival_time<-sapply(strsplit(data1$Survival_time,': '),'[',2)
table(data1$Relapse)

cel_files <- affy::list.celfiles("~/Brain/download/GSE12276/GSE12276_RAW/", 
                                 full.name = F)

cel_files <- cel_files [str_extract(cel_files,'GSM\\d+')%in%data1$geo_accession]
cel_files <- data.frame(filename=cel_files,
                        geo_accession=str_extract(cel_files,'GSM\\d+'))

data1<-data1 %>% inner_join(cel_files)

cel <- ReadAffy(filenames = data1$filename,
                celfile.path = "/Users/huzixin/Brain/download/GSE12276/GSE12276_RAW/",
                phenoData = data1)



eset<-rma(cel)
data1$Relapse_site<-ifelse(str_detect(data1$Relapse,'brain'),'brain','others')
qs::qsave(eset,'./GSE12276/eset204.qs')
qs::qsave(data1,'./GSE12276/pdata204.qs')
eset<-qs::qread('./eset204.qs')
target<-pData(eset)
data1<-qs::qread('./pdata204.qs')
exprSet<-exprs(eset)
colnames(exprSet)<-eset$geo_accession
boxplot(exprSet,outline=FALSE, notch=T, las=2)

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

table(data1$Relapse)
library(hgu133plus2.db)
probeid<-as.data.frame(toTable(hgu133plus2SYMBOL))

group<-as.factor(data1$Relapse_site)
design<-model.matrix(~0+group)
colnames(design)<-levels(group)
fit<-lmFit(eset,design)
contrast.matrix<-makeContrasts(brain-others,levels=design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
topTable(fit2,coef=1,adjust='BH')
results <- topTable(fit2,adjust='fdr',coef=1,number=Inf)
results <- results %>% 
  rownames_to_column("probe_id") %>%
  inner_join(probeid,by="probe_id") %>%
  dplyr::select(-probe_id) %>%
  dplyr::select(symbol,everything()) %>%
  arrange(symbol,desc(logFC)) %>% 
  distinct(symbol,.keep_all = T)

results<- results %>% 
  filter(P.Value<0.05)
dir.create('./Transcriptome/GSE12276')
data.table::fwrite(results,file ='./GSE12276/results_GSE12276.txt',
                   sep='\t',quote = F,row.names = F )

results<- results %>% 
  filter(P.Value<0.05&logFC>=0)

pdf('./GSE12276/result_GSE12276_volcano.pdf',height=6,width=6)
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
##survivalFigure8D=====
rm(list=ls());gc()
eset<-qs::qread('./eset204.qs')
pdata<-qs::qread('./pdata204.qs')
colnames(pdata)
library(hgu133plus2.db)
probeid<-as.data.frame(toTable(hgu133plus2SYMBOL))
exprSet<-exprs(eset)
colnames(exprSet)<-eset$geo_accession
edata1<-exprSet %>%as.data.frame() %>% 
  rownames_to_column('probe_id') %>% 
  inner_join(probeid) %>% 
  dplyr::select(-probe_id) %>%
  dplyr::select(c('symbol',everything()))
exp1 <- edata1 %>% 
  filter(symbol%in%c("SPP1","VSIG4")) %>% 
  group_by(symbol) %>% 
  summarise(across(starts_with("GSM"), mean)) %>% 
  ungroup() %>% 
  column_to_rownames('symbol') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('geo_accession') %>% 
  inner_join(pdata[,c('geo_accession',"Survival_time","Relapse_site")])%>% 
 #filter(Relapse_site=='brain') %>% 
  mutate(spp1group=ifelse(SPP1>mean(SPP1),'high','low'),
         vsig4group=ifelse(VSIG4>mean(VSIG4),'high','low'))
exp1$Status<-1
exp1$Survival_time<-as.numeric(exp1$Survival_time)

qs::qsave(exp1,'./exp_DAMgene_204samples.qs')
exp1<-qs::qread('../GSE12276/exp_DAMgene_204samples.qs')
table(exp1$Site)
table(exp$Site)
colnames(exp1)[4:5]<-c('Time','Site')
colnames(exp2)

sfit <- survfit(Surv(Time, Status)~spp1group, 
                data=exp1[exp1$Site=='brain',])
table(exp$Site)
summary(sfit)
pdf('./SPP1.survival.GSE12276.pdf',height=5,width=5)
ggsurvplot(sfit,pval = T,data=exp1[exp1$Site=='brain',],
           palette = c('#CC3333','#0099CC'),
           ggtheme = theme_bw(),
           risk.table = F,
           legend.title = "SPP1 group",
           legend.labs = c("High", "Low"),
           aspect.ratio = 1,
           legend.position = 'right',
           legend.text = element_text(size=12),
           #font.main = c(16, "bold", "black"),
           font.x = c(14, "bold.italic", "black"),
           font.y = c(14, "bold.italic", "black"),
           font.tickslab = c(12, "plain", "black"))

dev.off()

sfit <- survfit(Surv(Time, Status)~spp1group, 
                data=exp1)
summary(sfit)
pdf('./SPP1.survival.GSE1227.all.pdf',height=5,width=5)
ggsurvplot(sfit,pval = T,data=exp1,
           palette = c('#CC3333','#0099CC'),
           ggtheme = theme_bw(),
           risk.table = F,
           legend.title = "SPP1 group",
           legend.labs = c("High", "Low"),
           aspect.ratio = 1,
           legend.position = 'right',
           legend.text = element_text(size=12),
           #font.main = c(16, "bold", "black"),
           font.x = c(14, "bold.italic", "black"),
           font.y = c(14, "bold.italic", "black"),
           font.tickslab = c(12, "plain", "black"))

dev.off()