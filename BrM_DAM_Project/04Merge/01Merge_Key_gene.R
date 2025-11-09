library(GEOquery);library(tidyverse);library(arrayQualityMetrics);library(limma)
library(Seurat);library(harmony);library(scCustomize);library(scRNAtoolVis);library(clusterProfiler);library(UCell);library(biomaRt)
DAM<-data.table::fread('./Geneset/DAM_signature_hg.txt')
MGnD<-data.table::fread('./Geneset/MGnD_signature_hg.txt')
MG_dataset<-data.table::fread('./Geneset/Gan_et_al_MGdataet1_hg.txt')
DAM_gene<-MG_dataset$Gene[MG_dataset$Category%in%c('DAM_core_program','Stage2_DAM','Stage1_DAM_upregulated')] %>% unique()
DAM_gene<-c(MG_dataset$Gene[MG_dataset$Category%in%c('DAM_core_program','Stage2_DAM','Stage1_DAM_upregulated')] %>% unique(),
           DAM$Gene,
           MGnD$Gene) %>% unique()
##Dataset C6====
c6markers<-data.table::fread('../dataset3/V1/03ConvsMet/03markers_TAM_cluster_type.txt') %>% 
  filter(group=='Met'&cluster=='c6')
library(biomaRt)
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
gene<-c6markers$gene[c6markers$avg_log2FC>=1]
gene<- getLDS(attributes = c("external_gene_name"),
              filters = "external_gene_name",
              values = unique(gene),
              mart = mouse,
              attributesL = c("hgnc_symbol"),
              martL = human,uniqueRows = T
)
c6markers<-gene
saveRDS(c6markers,'./dataset3_DAMc6markers.rds')
c6markers<-readRDS('./DAM/dataset3_DAMc6markers.rds')


##DAM markers TAM celltype-----
TAM_markers.d1<-data.table::fread('./dataset1/Lung/04clusterMDM/TAMcelltype_cluster_markers_all.txt')
TAM_markers.d2<-data.table::fread('./dataset2/Lung/04clusterMDM/TAMcelltype_cluster_markers_all.txt')
TAM_markers.d6<-data.table::fread('./dataset6/Lung/04clusterMDM/TAMcelltype_cluster_markers_all.txt')
TAM_markers.l1<-data.table::fread('./Lung_to_brain/Lung1_GSE123902/03TAM_Brain/Brain_TAMcelltype_cluster_markers_all.txt')
TAM_markers.l2<-data.table::fread('./Lung_to_brain/Lung2_GSE202371/02TAM/TAMcelltype_cluster_markers_all.txt')
TAM_markers.l3<-data.table::fread('./Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/02TAM/TAMcelltype_cluster_markers_all.txt')
table(TAM_markers.d2$cluster)
table(TAM_markers.l3$cluster)
Key_Gene_DAM<-Reduce(intersect,list(TAM_markers.d6$gene[TAM_markers.d6$avg_log2FC>=1&TAM_markers.d6$cluster=='DAM.Microglia'],
                                    TAM_markers.d2$gene[TAM_markers.d2$avg_log2FC>=1&TAM_markers.d2$cluster=='DAM.Microglia'],
                                    TAM_markers.d1$gene[TAM_markers.d1$avg_log2FC>=1&TAM_markers.d1$cluster=='DAM.Microglia'],
                                   TAM_markers.l1$gene[TAM_markers.l1$avg_log2FC>=1&TAM_markers.l1$cluster=='DAM.Microglia'],
                                    TAM_markers.l2$gene[TAM_markers.l2$avg_log2FC>=1&TAM_markers.l2$cluster=='DAM.Microglia'],
                                    TAM_markers.l3$gene[TAM_markers.l3$avg_log2FC>=1&TAM_markers.l3$cluster=='DAM.Microglia'],
                                    DAM_gene,c6markers$HGNC.symbol));

saveRDS(Key_Gene_DAM,'./KEY_Gene_DAM.rds')
#"SPP1" "CCL3" "CCL4" "APOE" "CD74" "CD83"
Key_Gene_DAM<-Reduce(intersect,list(TAM_markers.d6$gene[TAM_markers.d6$avg_log2FC>=1&TAM_markers.d6$cluster=='DAM.Microglia'],
                                    TAM_markers.d2$gene[TAM_markers.d2$avg_log2FC>=1&TAM_markers.d2$cluster=='DAM.Microglia'],
                                    TAM_markers.d1$gene[TAM_markers.d1$avg_log2FC>=1&TAM_markers.d1$cluster=='DAM.Microglia'],
                                    TAM_markers.l1$gene[TAM_markers.l1$avg_log2FC>=1&TAM_markers.l1$cluster=='DAM.Microglia'],
                                    TAM_markers.l2$gene[TAM_markers.l2$avg_log2FC>=1&TAM_markers.l2$cluster=='DAM.Microglia'],
                                    TAM_markers.l3$gene[TAM_markers.l3$avg_log2FC>=1&TAM_markers.l3$cluster=='DAM.Microglia']));
Key_Gene_DAM
saveRDS(Key_Gene_DAM,'./KEY_Gene_DAM_ALL.rds')
# [1] "SPP1"     "CCL3"     "CCL4"     "C1QC"
# [5] "C1QB"     "C1QA"     "VSIG4"    "APOE"
# [9] "CD14"     "TREM2"    "FCGR3A"   "TMIGD3"
# [13] "RNASET2"  "CD74"     "FOLR2"    "GPR34"
# [17] "FCGR1A"   "HLA-DPA1" "LILRB4"   "MS4A7"
# [21] "LGMN"     "MS4A6A"   "OLR1"     "HLA-DRA"
# [25] "C3AR1"    "HLA-DRB1" "CD83"     "MAFB"
# [29] "CSF1R"    "HLA-DMA"  "SLCO2B1"  "SGK1"

##MG_markers all dataframe-----
colnames(TAM_markers.d2)[c(2,5)]<-paste0('d2.',colnames(TAM_markers.d2)[c(2,5)])
colnames(TAM_markers.d6)[c(2,5)]<-paste0('d6.',colnames(TAM_markers.d6)[c(2,5)])
colnames(TAM_markers.l1)[c(2,5)]<-paste0('lung1.',colnames(TAM_markers.l1)[c(2,5)])
colnames(TAM_markers.d1)[c(2,5)]<-paste0('d1.',colnames(TAM_markers.d1)[c(2,5)])
colnames(TAM_markers.l2)[c(2,5)]<-paste0('lung2.',colnames(TAM_markers.l2)[c(2,5)])
colnames(TAM_markers.l3)[c(2,5)]<-paste0('lung3.',colnames(TAM_markers.l3)[c(2,5)])

DAM_df<-TAM_markers.d1[TAM_markers.d1$cluster=='DAM.Microglia',c('gene','d1.avg_log2FC','d1.p_val_adj')] %>% 
  filter(gene%in%Key_Gene_DAM&d1.avg_log2FC>=1) %>% 
  left_join(TAM_markers.d2[TAM_markers.d2$d2.avg_log2FC>=1&TAM_markers.d2$cluster=='DAM.Microglia',c('d2.avg_log2FC','d2.p_val_adj','gene')]) %>% 
  left_join(TAM_markers.d6[TAM_markers.d6$d6.avg_log2FC>=1&TAM_markers.d6$cluster=='DAM.Microglia',c('d6.avg_log2FC','d6.p_val_adj','gene')]) %>% 
  left_join(TAM_markers.l1[TAM_markers.l1$lung1.avg_log2FC>=1&TAM_markers.l1$cluster=='DAM.Microglia',c('lung1.avg_log2FC','lung1.p_val_adj','gene')]) %>% 
  left_join(TAM_markers.l2[TAM_markers.l2$lung2.avg_log2FC>=1&TAM_markers.l2$cluster=='DAM.Microglia',c('lung2.avg_log2FC','lung2.p_val_adj','gene')]) %>% 
  left_join(TAM_markers.l3[TAM_markers.l3$lung3.avg_log2FC>=1&TAM_markers.l3$cluster=='DAM.Microglia',c('lung3.avg_log2FC','lung3.p_val_adj','gene')]) 
dir.create('./output_df/')
data.table::fwrite(DAM_df,'./output_df/')
data.table::fwrite(DAM_df,'./output_df/DAM_markers_32Genes.tsv',sep='\t',
                   quote = F,row.names = F)

##DAM-like.Mac Gene--------
markers_d1<-data.table::fread('./dataset1/Lung/06trajectory.new/TAMcelltypetra_markers_all.txt') 
markers_d2<-data.table::fread('./dataset2/Lung/06trajectory.new/TAMcelltypetra_markers_all.txt') 
markers_d6<-data.table::fread('./dataset6/Lung/06trajectory.new/TAMcelltypetra_markers_all.txt') 
markers_l1<-data.table::fread('./Lung_to_brain/Lung1_GSE123902/06trajectory.new/TAMcelltypetra_markers_all.txt') 
markers_l2<-data.table::fread('./Lung_to_brain/Lung2_GSE202371/06trajectory.new/TAMcelltypetra_markers_all.txt') 
markers_l3<-data.table::fread('./Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/06trajectory.new/TAMcelltypetra_markers_all.txt') 

Key_Gene_MDM.mac<-Reduce(intersect,list(markers_d1$gene[markers_d1$avg_log2FC>=1&markers_d1$cluster=='DAM-like.Mac.Terminal'],
                                        markers_d2$gene[markers_d2$avg_log2FC>=1&markers_d2$cluster=='DAM-like.Mac.Terminal'],
                                        markers_d6$gene[markers_d6$avg_log2FC>=1&markers_d6$cluster=='DAM-like.Mac.Terminal'],
                                        markers_l1$gene[markers_l1$avg_log2FC>=1&markers_l1$cluster=='DAM-like.Mac.Terminal'],
                                        markers_l2$gene[markers_l2$avg_log2FC>=1&markers_l2$cluster=='DAM-like.Mac.Terminal'],
                                        markers_l3$gene[markers_l3$avg_log2FC>=1&markers_l3$cluster=='DAM-like.Mac.Terminal']));Key_Gene_MDM.mac

# [1] "APOC1"    "APOE"     "GPNMB"    "CTSB"    
# [5] "FTL"      "C1QB"     "C1QA"     "C1QC"    
# [9] "TYROBP"   "PSAP"     "FCER1G"   "CD74"    
# [13] "HLA-DRB1" "HLA-DRA"  "HLA-DPA1" "CD163"   
# [17] "CTSS"     "TREM2"    "HLA-DPB1" "CAPG"    
# [21] "HLA-DQB1" "HLA-DMA"  "MS4A7"    "FCGR2A"  
# [25] "CD14"     "VSIG4"    "CST3"     "CCL4L2"  
# [29] "MS4A4A"   "AIF1"     "MS4A6A"   "DAB2"  
saveRDS(Key_Gene_MDM.mac,'./Key_Gene_DAM_MDM_Tra.rds')

##DAM.mac MDM data.frame------
colnames(markers_d1)[c(2,5)]<-paste0('d1.',colnames(markers_d1)[c(2,5)])
colnames(markers_d2)[c(2,5)]<-paste0('d2.',colnames(markers_d2)[c(2,5)])
colnames(markers_d6)[c(2,5)]<-paste0('d6.',colnames(markers_d6)[c(2,5)])
colnames(markers_l1)[c(2,5)]<-paste0('lung1.',colnames(markers_l1)[c(2,5)])
colnames(markers_l2)[c(2,5)]<-paste0('lung2.',colnames(markers_l2)[c(2,5)])
colnames(markers_l3)[c(2,5)]<-paste0('lung3.',colnames(markers_l3)[c(2,5)])
table(markers_d1$cluster)
colnames(markers_d1)
DAM_df<-markers_d1[markers_d1$cluster=='DAM-like.Mac.Terminal',c('gene','d1.avg_log2FC','d1.p_val_adj')] %>% 
  filter(gene%in%Key_Gene_MDM.mac&d1.avg_log2FC>=1) %>% 
  left_join(markers_d2[markers_d2$d2.avg_log2FC>=1&markers_d2$cluster=='DAM-like.Mac.Terminal',c('d2.avg_log2FC','d2.p_val_adj','gene')]) %>% 
  left_join(markers_d6[markers_d6$d6.avg_log2FC>=1&markers_d6$cluster=='DAM-like.Mac.Terminal',c('d6.avg_log2FC','d6.p_val_adj','gene')]) %>% 
  left_join(markers_l1[markers_l1$lung1.avg_log2FC>=1&markers_l1$cluster=='DAM-like.Mac.Terminal',c('lung1.avg_log2FC','lung1.p_val_adj','gene')]) %>% 
  left_join(markers_l2[markers_l2$lung2.avg_log2FC>=1&markers_l2$cluster=='DAM-like.Mac.Terminal',c('lung2.avg_log2FC','lung2.p_val_adj','gene')]) %>% 
  left_join(markers_l3[markers_l3$lung3.avg_log2FC>=1&markers_l3$cluster=='DAM-like.Mac.Terminal',c('lung3.avg_log2FC','lung3.p_val_adj','gene')]) 
data.table::fwrite(DAM_df,'./output_df/DAM.Mac.Terminal_markers_32Genes.tsv',sep='\t',
                   quote = F,row.names = F)



##DAM enrichment-----
dam_d1<-readRDS('./dataset1/Lung/04clusterMDM/DAM_enrichment_GSVA_d1.rds')
dam_d2<-readRDS('./dataset2/Lung/04clusterMDM/DAM_enrichment_GSVA_d2.rds')
dam_d6<-readRDS('./dataset6/Lung/04clusterMDM/DAM_enrichment_GSVA_d6.rds')
dam_l1<-readRDS('./Lung_to_brain/Lung1_GSE123902/03TAM_Brain/DAM_enrichment_GSVA_l1.rds')
dam_l2<-readRDS('./Lung_to_brain/Lung2_GSE202371/02TAM/DAM_enrichment_GSVA_l2.rds')
dam_l3<-readRDS('./Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/02TAM/DAM_enrichment_GSVA_l3.rds')
dam_d3<-readRDS('./dataset3/V2/DAM_enrichment_GSVA_d3.rds')

Reduce(intersect,list(dam_d1$ID[dam_d1$NES>0],
                      dam_d2$ID[dam_d2$NES>0],
                      dam_d6$ID[dam_d6$NES>0],
                      dam_l1$ID[dam_l1$NES>0],
                      dam_l2$ID[dam_l2$NES>0],
                      dam_l3$ID[dam_l3$NES>0],
                      dam_d3$ID[dam_d3$NES>0]))


#"HALLMARK_TNFA_SIGNALING_VIA_NFKB"
##DAM.Mac Enrichment-------
res_d1<-readRDS('./dataset1/Lung/06trajectory.new/DAM.Mac_enrichment_GSVA.rds') %>% arrange(desc(NES)) %>% filter(NES>=0)
res_d2<-readRDS('./dataset2/Lung/06trajectory.new/DAM.Mac_enrichment_GSVA.rds')%>% arrange(desc(NES)) %>% filter(NES>=0)
res_d6<-readRDS('./dataset6/Lung/06trajectory.new/DAM.Mac_enrichment_GSVA.rds')%>% arrange(desc(NES)) %>% filter(NES>=0)
res_l1<-readRDS('./Lung_to_brain/Lung1_GSE123902/06trajectory.new/DAM.Mac_enrichment_GSVA.rds')%>% arrange(desc(NES)) %>% filter(NES>=0)
res_l3<-readRDS('./Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907//06trajectory.new/DAM.Mac_enrichment_GSVA.rds')%>% arrange(desc(NES)) %>% filter(NES>=0)
res_l2<-readRDS('./Lung_to_brain/Lung2_GSE202371/06trajectory.new/DAM.Mac_enrichment_GSVA.rds')%>% arrange(desc(NES)) %>% filter(NES>=0)
DAM.mac.enrichment<-Reduce(intersect,list(res_d1$ID,
                                          # res_d2$ID,
                                          res_d6$ID,
                                          res_l1$ID,
                                          res_l2$ID,
                                          res_l3$ID))

# [1] "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
# [2] "HALLMARK_COMPLEMENT"             
# [3] "HALLMARK_KRAS_SIGNALING_UP"
saveRDS(DAM.mac.enrichment,'./DAM.mac.tra.enrichent.rds')
##DAM related regulon -----
rm(list=ls());gc()
library(tidyverse)
data.regulon.d1<-data.table::fread('./dataset1/dataset1Scenic/TAM.celltype_regulon.txt')
data.regulon.d2<-data.table::fread('./dataset2/dataset2Scenic/TAM.celltype_regulon.txt')
data.regulon.d6<-data.table::fread('./dataset6/dataset6Scenic/TAM.celltype_regulon.txt')
data.regulon.l1<-data.table::fread('./Lung_to_brain/Lung1_GSE123902/05Scenic/TAM.celltype_regulon.txt')
data.regulon.l2<-data.table::fread('./Lung_to_brain/Lung2_GSE202371/04Scenic/TAM.celltype_regulon.txt')
data.regulon.l3<-data.table::fread('./Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907//04Scenic/TAM.celltype_regulon.txt')
data.regulon.d3<-data.table::fread('./dataset3/V2/Scenic/dataset3_seu.TAM_TAMcluster_regulon.txt')
library(biomaRt)
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
gene<-sapply(strsplit(data.regulon.d3$gene,'\\('),'[',1)
gene<- getLDS(attributes = c("external_gene_name"),
              filters = "external_gene_name",
              values = gene,
              mart = mouse,
              attributesL = c("hgnc_symbol"),
              martL = human,
              uniqueRows = T
)
gene$TF<-paste0(gene$HGNC.symbol,'(+)')
colnames(gene)[1]<-'gene'
gene$gene<-paste0(gene$gene,'(+)')
data.regulon.d3<-data.regulon.d3 %>% 
  left_join(gene)
saveRDS(data.regulon.d3,'./DAM/data.regulon.d3.hgnc.rds')

data.regulon.d3<-readRDS('./DAM/data.regulon.d3.hgnc.rds')



DAM_regulon<-Reduce(intersect,list(head(data.regulon.d1$gene[order(desc(data.regulon.d1$DAM))&data.regulon.d1$celltype_pl>0.2],200),
                                   head(data.regulon.d2$gene[order(desc(data.regulon.d2$DAM))],200),
                                   head(data.regulon.d6$gene[order(desc(data.regulon.d6$DAM))],200),
                                  # head(data.regulon.l1$gene[order(desc(data.regulon.l1$DAM))&data.regulon.l1$celltype_pl>0.2],200),
                                   head(data.regulon.l2$gene[order(desc(data.regulon.l2$DAM))],200),
                                   head(data.regulon.l3$gene[order(desc(data.regulon.l3$DAM))],200),
                                   head(data.regulon.d3$TF[order(desc(data.regulon.d3$Metc6))],200)
));DAM_regulon

DAM_up_regulon<-Reduce(intersect,list(DAM_regulon,
                                      data.regulon.d3$TF[str_detect(data.regulon.d3$c6,'\\+')]))
#"ETS2(+)"  "FOS(+)"   "IRF1(+)"  "PRDM1(+)" "RUNX3(+)"
DAM_down_regulon<-Reduce(intersect,list(DAM_regulon,data.regulon.d3$TF[str_detect(data.regulon.d3$c6,'-')]))

# [1] "BCLAF1(+)" "CEBPA(+)"  "CEBPD(+)"  "ERF(+)"   
# [5] "FLI1(+)"   "IKZF1(+)"  "JUND(+)"   "KLF9(+)"  
# [9] "NR2F1(+)"  "NR3C1(+)"  "PURA(+)"   "SOX4(+)"  
# [13] "USF2(+)"   "XBP1(+)"   "ZBTB7A(+)"

saveRDS(DAM_up_regulon,'DAM_up_regulon.rds')
saveRDS(DAM_down_regulon,'DAM_down_regulon.rds')
saveRDS(DAM_regulon,'./DAM_regulon_byrssMat.rds')


DAM_regulon<-data.frame(TF=DAM_regulon,
                        gene=sapply(strsplit(DAM_regulon,'\\('),'[',1))

gene<- getLDS(attributes = c("external_gene_name"),
              filters = "external_gene_name",
              values = DAM_regulon$gene,
              mart = mouse,
              attributesL = c("hgnc_symbol"),
              martL = human,
              uniqueRows = T
)
eregulon.d3$ID<-sapply(strsplit(eregulon.d3$ID,'\\('),'[',1)
colnames(gene)[2]<-'gene'
DAM_regulon<-DAM_regulon %>% 
  left_join(gene)
DAM_regulon$Mgi_TF<-paste0(DAM_regulon$Gene.name,'(+)')
saveRDS(DAM_regulon,'./DAM_regulon_byrssMat_formgi.rds')

##eregulon DAM------
eregulon.d1<-readRDS('./dataset1/dataset1Scenic/KEY_Gene_DAM_regulon_dataset1.rds')
eregulon.d2<-readRDS('./dataset2/dataset2Scenic/KEY_Gene_DAM_regulon_dataset2.rds')
eregulon.d6<-readRDS('./dataset6/dataset6Scenic/KEY_Gene_DAM_regulon_dataset6.rds')

eregulon.l1<-readRDS('./Lung_to_brain/Lung1_GSE123902/05Scenic/KEY_Gene_DAM_regulon_Lung1.rds')
eregulon.l2<-readRDS('./Lung_to_brain/Lung2_GSE202371/04Scenic/KEY_Gene_DAM_regulon_Lung2.rds')
eregulon.l3<-readRDS('./Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/04Scenic/KEY_Gene_DAM_regulon_Lung3.rds')
eregulon.d3<-readRDS('./dataset3/V2/KEY_Gene_DAM_regulon_dataset3.rds')

gene<-sapply(strsplit(eregulon.d3$ID,'\\('),'[',1)
gene<- getLDS(attributes = c("external_gene_name"),
              filters = "external_gene_name",
              values = gene,
              mart = mouse,
              attributesL = c("hgnc_symbol"),
              martL = human,
              uniqueRows = T
)

eregulon.d3$ID<-sapply(strsplit(eregulon.d3$ID,'\\('),'[',1)
colnames(gene)[1]<-'ID'
eregulon.d3<-eregulon.d3 %>% 
  left_join(gene)

eregulon.d1$ID<-sapply(strsplit(eregulon.d1$ID,'\\('),'[',1)
eregulon.d2$ID<-sapply(strsplit(eregulon.d2$ID,'\\('),'[',1)
eregulon.d6$ID<-sapply(strsplit(eregulon.d6$ID,'\\('),'[',1)
eregulon.l1$ID<-sapply(strsplit(eregulon.l1$ID,'\\('),'[',1)
eregulon.l2$ID<-sapply(strsplit(eregulon.l2$ID,'\\('),'[',1)
eregulon.l3$ID<-sapply(strsplit(eregulon.l3$ID,'\\('),'[',1)
DAM_eregulon<-Reduce(intersect,list(eregulon.d1$ID,
                                    eregulon.d2$ID,
                                    eregulon.d6$ID,
                                    eregulon.l1$ID,
                                    eregulon.l2$ID,
                                    eregulon.l3$ID,
                                    eregulon.d3$HGNC.symbol
))
#[1] "SPI1"  "IRF8"  "IRF5"  "REL"   "RUNX3" "FLI1"  "CEBPB" "STAT1"
DAM_eregulon_regulon<-intersect(paste0(DAM_eregulon,'(+)'),DAM_regulon$TF)
# "SPI1(+)"  "IRF8(+)"  "IRF5(+)"  "REL(+)"  "RUNX3(+)" "FLI1(+)" 
saveRDS(DAM_eregulon_regulon,'./DAM_regulon_byeregulon_regulon.rds')




##DAM.mac regulon=======
DAM.mac.d1<-data.table::fread('./dataset1/dataset1Scenic/trajectory_DAM-like.txt') %>% arrange(desc(`DAM-like.Mac`))
DAM.mac.d2<-data.table::fread('./dataset2/dataset2Scenic/trajectory_DAM-like.txt')%>% arrange(desc(`DAM-like.Mac`))
DAM.mac.d6<-data.table::fread('./dataset6/dataset6Scenic/trajectory_DAM-like.txt') %>% arrange(desc(`DAM-like.Mac`))
DAM.mac.l1<-data.table::fread('./Lung_to_brain/Lung1_GSE123902/05Scenic/trajectory_DAM-like.txt')%>% arrange(desc(`DAM-like.Mac`))
DAM.mac.l2<-data.table::fread('./Lung_to_brain/Lung2_GSE202371/04Scenic/trajectory_DAM-like.txt')%>% arrange(desc(`DAM-like.Mac`))
DAM.mac.l3<-data.table::fread('./Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907//04Scenic/trajectory_DAM-like.txt')%>% arrange(desc(`DAM-like.Mac`))
Trajectory.Regulon<-Reduce(intersect,list(head(DAM.mac.d1$gene,50),
                                          #head(DAM.mac.d2$gene,50),
                                          head(DAM.mac.d6$gene,50),
                                          #head(DAM.mac.l1$gene,50),
                                          head(DAM.mac.l2$gene,50),
                                          head(DAM.mac.l3$gene,50)));Trajectory.Regulon
#[1] "SPI1(+)"  "IRF5(+)"  "MAFB(+)"  "CEBPB(+)" "FLI1(+)"  "IRF8(+)"  "NFIL3(+)" "KLF4(+)"  "CREM(+)"  "MAF(+)"  
Trajectory.Regulon<-sub('\\(\\+\\)','',Trajectory.Regulon)
Key_Gene_MDM.mac<-Reduce(intersect,list(markers_d1$gene[markers_d1$avg_log2FC>0&markers_d1$cluster=='DAM-like.Mac.Terminal'],
                                        #markers_d2$gene[markers_d2$avg_log2FC>0&markers_d2$cluster=='DAM-like.Mac.Terminal'],
                                        markers_d6$gene[markers_d6$avg_log2FC>0&markers_d6$cluster=='DAM-like.Mac.Terminal'],
                                        # markers_l1$gene[markers_l1$avg_log2FC>0&markers_l1$cluster=='DAM-like.Mac.Terminal'],
                                        markers_l2$gene[markers_l2$avg_log2FC>0&markers_l2$cluster=='DAM-like.Mac.Terminal'],
                                        markers_l3$gene[markers_l3$avg_log2FC>0&markers_l3$cluster=='DAM-like.Mac.Terminal']));
Key.Trajectory.Regulon<-intersect(Key_Gene_MDM.mac,Trajectory.Regulon)
#"SPI1" "KLF4" "MAFB" "IRF8"
saveRDS(Key.Trajectory.Regulon,'./DAM.MDM.Terminal.regulon.rds')

##ereguon DAM.mac-----
eregulon.d1<-readRDS('../dataset1/dataset1Scenic/KEY_Gene_DAM-like.Mac_regulon_dataset1.rds')
eregulon.d2<-readRDS('../dataset2/dataset2Scenic/KEY_Gene_DAM-like.Mac_regulon_dataset2.rds')
eregulon.d6<-readRDS('../dataset6/dataset6Scenic/KEY_Gene_DAM-like.Mac_regulon_dataset6.rds')

eregulon.l1<-readRDS('../Lung_to_brain/Lung1_GSE123902/05Scenic/KEY_Gene_DAM-like.Mac_regulon_Lung1.rds')
eregulon.l2<-readRDS('../Lung_to_brain/Lung2_GSE202371/04Scenic/KEY_Gene_DAM-like.Mac_regulon_Lung2.rds')
eregulon.l3<-readRDS('../Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/04Scenic/KEY_Gene_DAM-like.Mac_regulon_Lung3.rds')
eregulon.d1$ID<-sapply(strsplit(eregulon.d1$ID,'\\('),'[',1)
eregulon.d2$ID<-sapply(strsplit(eregulon.d2$ID,'\\('),'[',1)
eregulon.d6$ID<-sapply(strsplit(eregulon.d6$ID,'\\('),'[',1)
eregulon.l1$ID<-sapply(strsplit(eregulon.l1$ID,'\\('),'[',1)
eregulon.l2$ID<-sapply(strsplit(eregulon.l2$ID,'\\('),'[',1)
eregulon.l3$ID<-sapply(strsplit(eregulon.l3$ID,'\\('),'[',1)
DAM.Mac_eregulon<-Reduce(intersect,list(eregulon.d1$ID,
                                        eregulon.d2$ID,
                                        eregulon.d6$ID,
                                        eregulon.l1$ID,
                                        eregulon.l2$ID,
                                        eregulon.l3$ID
))
# [1] "SPI1"  "IRF8"  "IRF5"  "IKZF1" "USF2" 
# [6] "MAFB"  "ELF1"  "FLI1"  "RUNX3" "STAT1"
# [11] "MAF"   "IRF2" 
##Regulon of TNFa via NFkB=====
TNFregulon_d1<-readRDS('./dataset1/dataset1Scenic/TNF_TFs_d1.rds')
TNFregulon_d2<-readRDS('./dataset2/dataset2Scenic/TNF_TFs_d2.rds')
TNFregulon_d6<-readRDS('./dataset6/dataset6Scenic/TNF_TFs_d6.rds')
TNFregulon_l1<-readRDS('./Lung_to_brain/Lung1_GSE123902/05Scenic/TNF_TFs_l1.rds')
TNFregulon_l2<-readRDS('./Lung_to_brain/Lung2_GSE202371/04Scenic/TNF_TFs_l2.rds')
TNFregulon_l3<-readRDS('./Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/./04Scenic/TNF_TFs_l3.rds')
TNFregulon_d1$ID<-sapply(strsplit(TNFregulon_d1$ID,'\\('),'[',1)
TNFregulon_d2$ID<-sapply(strsplit(TNFregulon_d2$ID,'\\('),'[',1)
TNFregulon_d6$ID<-sapply(strsplit(TNFregulon_d6$ID,'\\('),'[',1)
TNFregulon_l1$ID<-sapply(strsplit(TNFregulon_l1$ID,'\\('),'[',1)
TNFregulon_l2$ID<-sapply(strsplit(TNFregulon_l2$ID,'\\('),'[',1)
TNFregulon_l3$ID<-sapply(strsplit(TNFregulon_l3$ID,'\\('),'[',1)
TNFregulon_d3<-readRDS('./dataset3/V2/TNF_TFs_dataset3.rds')
gene<-sapply(strsplit(TNFregulon_d3$ID,'\\('),'[',1)

gene<- getLDS(attributes = c("external_gene_name"),
              filters = "external_gene_name",
              values = gene,
              mart = mouse,
              attributesL = c("hgnc_symbol"),
              martL = human,
              uniqueRows = T
)

TNFregulon_d3$ID<-sapply(strsplit(TNFregulon_d3$ID,'\\('),'[',1)
colnames(gene)[1]<-'ID'
TNFregulon_d3<-TNFregulon_d3 %>% 
  left_join(gene)

TNF_eregulon<-Reduce(intersect,list(TNFregulon_d1$ID,
                                    TNFregulon_d2$ID,
                                    TNFregulon_d6$ID,
                                    TNFregulon_l1$ID,
                                    TNFregulon_l2$ID,
                                    TNFregulon_l3$ID
                                    ,TNFregulon_d3$HGNC.symbol
))

DAM_TNF_key_regulon<-intersect(paste0(intersect(TNF_eregulon,DAM_eregulon),'(+)'),
                               DAM_regulon$TF)
# "REL(+)"   "IRF8(+)"  "SPI1(+)" 
# "RUNX3(+)" "FLI1(+)"  "IRF5(+)" 
DAM_TNF_regulon<-intersect(paste0(TNF_eregulon),
                           DAM_eregulon)
saveRDS(DAM_TNF_regulon,'./DAM_TNF_regulon.rds')
# "REL"   "CEBPB" "IRF8"  "SPI1"  "RUNX3" "FLI1" "STAT1" "IRF5" 
DAM_TNF_key_regulon<-intersect(paste0(intersect(TNF_eregulon,DAM_eregulon),'(+)'),
                               DAM_regulon$TF)
# "REL(+)"   "IRF8(+)"  "SPI1(+)"  "RUNX3(+)" "FLI1(+)"  "IRF5(+)" 

##Cellchat Tcell and TAMs----
cc.d1<-data.table::fread('./dataset1/Lung/05cellchat/new/df_net_DAM_T.new.txt') %>% arrange(desc(prob))
cc.d2<-data.table::fread('./dataset2/Lung/05cellchat/new/df_net_DAM_T.new.txt') %>% arrange(desc(prob))
cc.d6<-data.table::fread('./dataset6/Lung/05cellchat/new/df_net_DAM_T.new.txt') %>% arrange(desc(prob))

cc.l1<-data.table::fread('./Lung_to_brain/Lung1_GSE123902/06cellchat/new/df.net_brain_solo_Tcell.txt') %>% arrange(desc(prob))
cc.l2<-data.table::fread('./Lung_to_brain/Lung2_GSE202371/03cellchat/new/df_net_DAM_T.new.txt') %>% arrange(desc(prob))
cc.l3<-data.table::fread('./Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/03cellchat/new/df_net_DAM_T.new.txt') %>% arrange(desc(prob))

###DAM and Tcell(FigureS7C data)=====
cc.d1.f2<-cc.d1 %>% filter(source=='DAM'&target%in%c('CD4.T.cm','CD4.T','CD4.T.cm','CD8.T','CD8.T.em','Treg','DC'))
cc.d2.f2<-cc.d2 %>% filter(source=='DAM'&target%in%c('CD4.T.cm','CD8.T','CD8.T.ex','Treg','DC'))
cc.d6.f2<-cc.d6 %>% filter(source=='DAM'& target%in%c('CD8.T','Treg','DC'))
cc.l1.f2<-cc.l1 %>% filter(source=='DAM'&target%in%c('CD4.T.cm','CD8.T.em','CD8.T.ex','Treg','DC'))
cc.l2.f2<-cc.l2 %>% filter(source=='DAM'&target%in%c('CD4.T.cm','CD8.T','CD8.T.em','CD8.T.ex','Treg','DC'))
cc.l3.f2<-cc.l3 %>% filter(source=='DAM'&target%in%c('CD4.T','CD4.T.activate','CD8.T','CD8.T.em','CD8.T.ex','Treg','DC'))
pathway<-Reduce(intersect,list(cc.d1.f2$interaction_name_2,
                               cc.d2.f2$interaction_name_2,
                               #cc.d3$interaction_name,
                               cc.d6.f2$interaction_name_2,
                               #cc.l1.f2$interaction_name,
                               cc.l2.f2$interaction_name_2,
                               cc.l3.f2$interaction_name_2))
# [1] "LGALS9 - CD45"        "PGE2-PTGES3 - PTGER4"
# [3] "TNF - TNFRSF1B"       "LGALS9 - HAVCR2"     
# [5] "GAS6 - AXL"           "ADORA3 - ENTPD1"   
cc.d1.f3<-cc.d1 %>% filter(source=='DAM-like.Mac'&target%in%c('CD4.T.cm','CD4.T','CD4.T.cm','CD8.T','CD8.T.em','Treg','DC'))
cc.d2.f3<-cc.d2 %>% filter(source=='DAM-like.Mac'&target%in%c('CD4.T.cm','CD8.T','CD8.T.ex','Treg','DC'))
cc.d6.f3<-cc.d6 %>% filter(source=='DAM-like.Mac'& target%in%c('CD8.T','Treg','DC'))
cc.l1.f3<-cc.l1 %>% filter(source=='DAM-like.Mac'&target%in%c('CD4.T.cm','CD8.T.em','CD8.T.ex','Treg','DC'))
cc.l2.f3<-cc.l2 %>% filter(source=='DAM-like.Mac'&target%in%c('CD4.T.cm','CD8.T','CD8.T.em','CD8.T.ex','Treg','DC'))
cc.l3.f3<-cc.l3 %>% filter(source=='DAM-like.Mac'&target%in%c('CD4.T','CD4.T.activate','CD8.T','CD8.T.em','CD8.T.ex','Treg','DC'))
pathway<-Reduce(intersect,list(cc.d1.f3$interaction_name_2,
                               cc.d2.f3$interaction_name_2,
                               cc.d6.f3$interaction_name_2,
                               cc.l1.f3$interaction_name_2,
                               cc.l2. dataf3$interaction_name_2,
                               cc.l3.f3$interaction_name_2))

# [1] "LGALS9 - CD45"   "LGALS9 - HAVCR2"
# [3] "APP - CD74"      "CD99 - CD99" 

##Cellchat TAM&Tumor cell-----
cc.d1<-data.table::fread('../dataset1/Lung/05cellchat/new/df_net_DAM_MTC.new.txt') %>% arrange(desc(prob))
cc.d2<-data.table::fread('../dataset2/Lung/05cellchat/new/df_net_DAM_MTC.new.txt') %>% arrange(desc(prob))
cc.d6<-data.table::fread('../dataset6/Lung/05cellchat/new/df_net_DAM_MTC.new.txt') %>% arrange(desc(prob))

cc.l1<-data.table::fread('../Lung_to_brain/Lung1_GSE123902/06cellchat/new/df_net_DAM_MTC.new.txt') %>% arrange(desc(prob))
cc.l2<-data.table::fread('../Lung_to_brain/Lung2_GSE202371/03cellchat/new/df_net_DAM_MTC.new.txt') %>% arrange(desc(prob))
cc.l3<-data.table::fread('../Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/03cellchat/new/df_net_DAM_MTC.new.txt') %>% arrange(desc(prob))

table(cc.l1$source)
###DAM and tumor (Figure7A,B)-----
table(cc)
cc.d1.f1<-cc.d1 %>% filter(target=='DAM'&source%in%c('EMT_high','EMT_low','REP_high'))
cc.d2.f1<-cc.d2 %>% filter(target=='DAM'&str_detect(source,'EMT'))
cc.d6.f1<-cc.d6 %>% filter(target=='DAM'&source%in%c('EMT_high','EMT_high_hypoxia','EMT_low','Rep_high'))
cc.l1.f1<-cc.l1 %>% filter(target=='DAM'&str_detect(source,'EMT'))
cc.l2.f1<-cc.l2 %>% filter(target=='DAM'&str_detect(source,'EMT'))
cc.l3.f1<-cc.l3 %>% filter(target=='DAM'&str_detect(source,'EMT'))

cc.d1.f1<-cc.d1 %>% filter(target=='DAM-like.Mac'&source%in%c('EMT_high','EMT_low','REP_high'))
cc.d2.f1<-cc.d2 %>% filter(target=='DAM-like.Mac'&str_detect(source,'EMT'))
cc.d6.f1<-cc.d6 %>% filter(target=='DAM-like.Mac'&source%in%c('EMT_high','EMT_high_hypoxia','EMT_low','Rep_high'))
cc.l1.f1<-cc.l1 %>% filter(target=='DAM-like.Mac'&str_detect(source,'EMT'))
cc.l2.f1<-cc.l2 %>% filter(target=='DAM-like.Mac'&str_detect(source,'EMT'))
cc.l3.f1<-cc.l3 %>% filter(target=='DAM-like.Mac'&str_detect(source,'EMT'))

pathway<-Reduce(intersect,list(cc.d1.f1$interaction_name_2,
                               cc.d2.f1$interaction_name_2,
                               cc.d6.f1$interaction_name_2,
                               # cc.l1.f1$interaction_name_2,
                               cc.l2.f1$interaction_name_2,
                               cc.l3.f1$interaction_name_2));pathway

cc.d1.f1<-cc.d1.f1 %>% filter(interaction_name_2%in%pathway)
cc.d2.f1<-cc.d2.f1 %>% filter(interaction_name_2%in%pathway)
cc.d6.f1<-cc.d6.f1 %>% filter(interaction_name_2%in%pathway)
cc.l1.f1<-cc.l1.f1 %>% filter(interaction_name_2%in%pathway)
cc.l2.f1<-cc.l2.f1 %>% filter(interaction_name_2%in%pathway)
cc.l3.f1<-cc.l3.f1 %>% filter(interaction_name_2%in%pathway)
# "APP - CD74"          
# "APP - (TREM2+TYROBP)"
# "MDK - LRP1"
# "MDK - NCL" 
# "CD99 - CD99"
# "CD99 - PILRA"
# "ICAM1 - (ITGAM+ITGB2)"
# "LAMA5 - CD44"
# "COL1A1 - CD44"
cc.d1.f1<-cc.d1 %>% filter(source=='DAM'&target%in%c('EMT_high','EMT_low','REP_high'))
cc.d2.f1<-cc.d2 %>% filter(source=='DAM'&str_detect(target,'EMT'))
cc.d6.f1<-cc.d6 %>% filter(source=='DAM'&target%in%c('EMT_high','EMT_high_hypoxia','EMT_low','Rep_high'))
cc.l1.f1<-cc.l1 %>% filter(source=='DAM'&str_detect(target,'EMT'))
cc.l2.f1<-cc.l2 %>% filter(source=='DAM'&str_detect(target,'EMT'))
cc.l3.f1<-cc.l3 %>% filter(source=='DAM'&str_detect(target,'EMT'))

cc.d1.f1<-cc.d1 %>% filter(source=='DAM-like.Mac'&target%in%c('EMT_high','EMT_low','REP_high'))
cc.d2.f1<-cc.d2 %>% filter(source=='DAM-like.Mac'&str_detect(target,'EMT'))
cc.d6.f1<-cc.d6 %>% filter(source=='DAM-like.Mac'&target%in%c('EMT_high','EMT_high_hypoxia','EMT_low','Rep_high'))
cc.l1.f1<-cc.l1 %>% filter(source=='DAM-like.Mac'&str_detect(target,'EMT'))
cc.l2.f1<-cc.l2 %>% filter(source=='DAM-like.Mac'&str_detect(target,'EMT'))
cc.l3.f1<-cc.l3 %>% filter(source=='DAM-like.Mac'&str_detect(target,'EMT'))


pathway<-Reduce(intersect,list(cc.d1.f1$interaction_name_2,
                               cc.d2.f1$interaction_name_2,
                               cc.d6.f1$interaction_name_2,
                               # cc.l1.f1$interaction_name_2,
                               cc.l2.f1$interaction_name_2,
                               cc.l3.f1$interaction_name_2));pathway


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
saveRDS(DAM.MTC.pathway,'./DAM.MTC.pathway.add.rds')
###DAM and TAMs(Figure7C)-----
table(cc.d1$source)
cc.d1.f1<-cc.d1 %>% filter(target=='DAM'&source%in%c('DAM','DAM-like.Mac','APOE.Mac'))
cc.d2.f1<-cc.d2 %>% filter(target=='DAM'&source%in%c('DAM','DAM-like.Mac','APOE.Mac'))
cc.d6.f1<-cc.d6 %>% filter(target=='DAM'&source%in%c('DAM','DAM-like.Mac','APOE.Mac'))
cc.l1.f1<-cc.l1 %>% filter(target=='DAM'&source%in%c('DAM','DAM-like.Mac','APOE.Mac'))
cc.l2.f1<-cc.l2 %>% filter(target=='DAM'&source%in%c('DAM','DAM-like.Mac','APOE.Mac'))
cc.l3.f1<-cc.l3 %>% filter(target=='DAM'&source%in%c('DAM','DAM-like.Mac','APOE.Mac'))

cc.d1.f1<-cc.d1 %>% filter(target=='DAM-like.Mac'&source%in%c('DAM','DAM-like.Mac','APOE.Mac'))
cc.d2.f1<-cc.d2 %>% filter(target=='DAM-like.Mac'&source%in%c('DAM','DAM-like.Mac','APOE.Mac'))
cc.d6.f1<-cc.d6 %>% filter(target=='DAM-like.Mac'&source%in%c('DAM','DAM-like.Mac','APOE.Mac'))
cc.l1.f1<-cc.l1 %>% filter(target=='DAM-like.Mac'&source%in%c('DAM','DAM-like.Mac','APOE.Mac'))
cc.l2.f1<-cc.l2 %>% filter(target=='DAM-like.Mac'&source%in%c('DAM','DAM-like.Mac','APOE.Mac'))
cc.l3.f1<-cc.l3 %>% filter(target=='DAM-like.Mac'&source%in%c('DAM','DAM-like.Mac','APOE.Mac'))



pathway<-Reduce(intersect,list(cc.d1.f1$interaction_name_2,
                               cc.d2.f1$interaction_name_2,
                               cc.d6.f1$interaction_name_2,
                               cc.l1.f1$interaction_name_2,
                               cc.l2.f1$interaction_name_2,
                               cc.l3.f1$interaction_name_2))

TAM.pathway<-c("APOE - TREM2","SPP1 - CD44",
               "LGALS9 - CD44",
               "LGALS9 - P4HB",
               "LGALS9 - HAVCR2","LGALS9 - CD45",
               "CD99 - CD99","CD99 - PILRA",
               "ANXA1 - FPR1","LAIR1 - LILRB4")
saveRDS(DAM.MTC.pathway,'./DAM.MTC.pathway.rds')
saveRDS(TAM.pathway,'./TAM.pathway.rds')

##Combine Transcriptome====
##BRMt transcriptome-----
diff_brmt_l1<-data.table::fread('./Lung_to_brain/Lung1_GSE123902/01merge/Diff_brainvsprimary_allcelltype.txt')
table(diff_brmt_l1$group)
diff_brmt_l1<-diff_brmt_l1 %>% 
  filter(group=='TAM'&avg_log2FC>=0)
diff_brmt_l3<-data.table::fread('./Lung_to_brain/Lung3_GSE131907/Lung3_GSE131907/05addPrimary/Diff_brainvsprimary_allcelltype.txt')
diff_brmt_l3<-diff_brmt_l3 %>% 
  filter(group=='TAM_DC'&avg_log2FC>=0)
results_GSE100534<-data.table::fread('./Transcriptome/GSE100534/results_GSE100534.txt')
results_GSE125989<-data.table::fread('./Transcriptome/GSE125989/results_GSE125989.txt')
results_GSE248830<-data.table::fread('./Transcriptome/GSE248830/results_GSE248830.txt')
results_GSE14690<-data.table::fread(file ='./Transcriptome/GSE14690/result_GSE14690.txt')

Brmt<-Reduce(intersect,list(results_GSE100534$symbol[results_GSE100534$logFC>0],
                            results_GSE125989$symbol[results_GSE125989$logFC>0],
                            results_GSE248830$symbol[results_GSE248830$logFC>0],
                            results_GSE14690$SYMBOL[results_GSE14690$logFC>0]));Brmt
#SPP1
intersect(Key_Gene_DAM_MDM,Brmt)
Reduce(intersect,list(diff_brmt_l1$gene,
                      diff_brmt_l3$gene,
                      Brmt))
#SPP1
saveRDS(Brmt,'./Brmt.Trans.rds')


