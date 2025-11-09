rm(list=ls())
gc()
library(Seurat);library(tidyverse);library(Nebulosa);
library(clustree);library(GEOquery);library(scCustomize);library(scRNAtoolVis)
library(UCell);library(clusterProfiler);library(harmony)
seu<-qs::qread('./01merge/seu.anno.qs')
rownames(seu@meta.data ) %>% head()
### makeMetaCells----

makeMetaCells <- function(seu, min.cells=10, reduction="umap", dims=1:2, k.param=10, cores=20) {
  seu <- seu %>%
    FindNeighbors(reduction = reduction, dims=dims, k.param=k.param) %>%
    FindClusters(res=50)
  metadata <- seu@meta.data
  metadata$METACELL_ID <- factor(metadata$seurat_clusters)
  dge_mat <- seu[["RNA"]]@counts
  
  dge_mat_mc <- parallel::mclapply(levels(metadata$METACELL_ID), function(xx) {
    cells <- rownames(subset(metadata, METACELL_ID == xx))
    Matrix::rowSums(dge_mat[, cells])
  }, mc.cores = cores)
  dge_mat_mc <- do.call(cbind, dge_mat_mc)
  
  metacell_metadata <-
    metadata[["METACELL_ID"]] %>%
    table() %>%
    as.data.frame()
  colnames(metacell_metadata) <-   c("METACELL_ID", "CELL_COUNT")
  rownames(metacell_metadata) <- metacell_metadata[["METACELL_ID"]]
  
  kept.cells <- subset(metacell_metadata, CELL_COUNT>=min.cells)[["METACELL_ID"]]
  metacells <- list(
    mat = dge_mat_mc[, kept.cells],
    metadata = metacell_metadata[kept.cells, ]
  )
  colnames(metacells$mat) <- paste0(seu@project.name, ".METACELL_", kept.cells)
  rownames(metacells$metadata) <- colnames(metacells$mat)
  metacells
}

DefaultAssay(seu)<-"RNA"
DimPlot(seu, group.by = "seurat_clusters")
DimPlot(seu, group.by = "celltype",label=T)
DimPlot(seu, group.by = "celltype",split.by = 'location')+NoLegend()
metacells <-makeMetaCells(
    seu       = seu,
    min.cells = 10,
    reduction = "umap",
    dims      = 1:2,
    k.param   = 10,
    cores     = 20)

mc.mat <- metacells$mat
summary(metacells$metadata$CELL_COUNT)
dir.create('./04Scenic/output')
saveRDS(mc.mat, "./04Scenic/input/Lung2.mat.rds")
## (1) TF list文件
motif2tfs <- data.table::fread("/Users/huzixin/Brain/brian/dataset1/dataset1Scenic/motifs_v10nr_clust_nr.hgnc_m0.001_o0.0.tbl")
TFs <- sort(unique(motif2tfs$gene_name))
writeLines(TFs, "./04Scenic/cisTarget_db_hsa_hgnc_tfs.txt")

## (2) meta cell matrix (for step1): *.csv or *.loom
mc.mat <- readRDS("./04Scenic/input/Lung2.mat.rds")
expr.in.cells <- rowSums(mc.mat > 0)
mc.mat <- mc.mat[expr.in.cells >= 5, ]
cisdb <- arrow::read_feather("/Users/huzixin/Brain/brian/dataset3/dataset3Scenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
dim(cisdb)
genes.use <- intersect(colnames(cisdb), rownames(mc.mat))
length(genes.use)
dim(mc.mat)
mc.mat <- mc.mat[genes.use, ]
dim(mc.mat)
loom <- SCopeLoomR::build_loom(
  file.name         = "./04Scenic/input/Lung2_mat_for_step1.loom",
  dgem              = mc.mat,
  default.embedding = NULL
)
loom$close()
rm(loom)
#sh------
cd ./04Scenic
conda activate pyscenic_env
## inputs
f_loom_grn=./input/Lung2_mat_for_step1.loom
## outputs
grn_output=./output/s1_adj.tsv
ctx_output=./output/s2_reg.tsv
dir.create('./05Scenic/output')


## reference
f_tfs=./cisTarget_db_hsa_hgnc_tfs.txt
f_motif_path=./motifs_v10nr_clust_nr.hgnc_m0.001_o0.0.tbl
f_db_names=`find . -name "hg38*.feather"`


#### 1. build GRN
arboreto_with_multiprocessing.py \
$f_loom_grn \
$f_tfs \
--method grnboost2 \
--output $grn_output \
--num_workers 10 \
--seed 777


#### 2. cisTarget
pyscenic ctx \
$grn_output \
$f_db_names \
--annotations_fname $f_motif_path \
--expression_mtx_fname $f_loom_grn \
--output $ctx_output \
--num_workers 10




##Define ComputeModuleScore--------
ComputeModuleScore <- function(x, ...) UseMethod('ComputeModuleScore')

ComputeModuleScore.default <- function(x, gene.sets, min.size=20, batch.size=500, cores=1, ...) {
  if (!is.list(gene.sets)) {
    stop("'gene.sets' should be a list or data.frame!")
  }
  gene.sets <- gene.sets[sapply(gene.sets, length) >= min.size]
  n.cells <- ncol(x)
  batches <- floor((1:n.cells-1) / batch.size)
  batch.levels <- unique(batches)
  aucell <- function(i) {
    dge.tmp <- x[, batches == i]
    cr <- AUCell::AUCell_buildRankings(dge.tmp, nCores=1, plotStats=F, verbose = F)
    auc <- AUCell::AUCell_calcAUC(gene.sets, cr, nCores=1, verbose = F)
    AUCell::getAUC(auc)
  }
  auc_scores <- parallel::mclapply(batch.levels, aucell, mc.cores = cores)
  do.call(cbind, auc_scores)
}

ComputeModuleScore.Seurat <- function(x, gene.sets, 
                                      min.size=20, batch.size=500, cores=1, 
                                      assay = Seurat::DefaultAssay(x), ...) {
  dge <- x[[assay]]@counts
  ras_mat <- ComputeModuleScore.default(x = dge, gene.sets, min.size, batch.size, cores)
  x[["AUCell"]] <- Seurat::CreateAssayObject(data = ras_mat)
  return(x)
}
##AUCell-----

regulons <- clusterProfiler::read.gmt("./04Scenic/output/Lung2BRME.regulons.gmt")

rg.names <- unique(regulons$term)
regulon.list <- lapply(rg.names, function(rg) {
  subset(regulons, term == rg)$gene
})
names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
summary(sapply(regulon.list, length))
saveRDS(regulon.list, "04Scenic/output/regulons.rds")

seu.TAM<-readRDS('./Figuresubmit/seu.TAM.celltype_pl.rds')
DefaultAssay(seu.TAM) <- "RNA"
seu.TAM <- ComputeModuleScore(seu.TAM, 
                              gene.sets = regulon.list, 
                              min.size = 10, 
                              cores = 10)
DefaultAssay(seu.TAM)<-'AUCell'
head(rownames(seu.TAM))

seu.TAM <- RunUMAP(object = seu.TAM,
                   features = rownames(seu.TAM),
                   metric = "correlation", # 注意这里用correlation效果最好
                   reduction.name = "umapRAS",
                   reduction.key = "umapRAS_")

seu.TAM <- ScaleData(seu.TAM)
seu.TAM <- RunPCA(object = seu.TAM,
                  features = rownames(seu.TAM),
                  reduction.name = "pcaRAS",
                  reduction.key = "pcaRAS_")

table(seu.TAM$celltype_pl %>% is.na())
(DimPlot(seu.TAM, reduction = "umap", group.by = "celltype_pl",label = T) + ggsci::scale_color_d3("category20")) +
  (DimPlot(seu.TAM, reduction = "umapRAS", group.by = "celltype_pl",label = T) + ggsci::scale_color_d3("category20")) 
DimPlot(seu.TAM, reduction = "pcaRAS", group.by = "celltype_pl") + ggsci::scale_color_d3("category20")
seu.TAM[["pcaRAS"]]@stdev ^ 2 / sum(seu.TAM[["pcaRAS"]]@stdev ^ 2)
ElbowPlot(seu.TAM, reduction = "pcaRAS")
Loadings(seu.TAM, reduction = "pcaRAS")
VizDimLoadings(seu.TAM, reduction = "pcaRAS", dims = 2, balanced = T, projected = F)

qs::qsave(seu.TAM, "./04Scenic/seu.TAM_addScenic.qs")

##寻找与MET相关的TF方差分解-------
VarDecompose <- function(data, meta.data, vd.vars, genes = "all", cores = -1) {
  
  if (length(genes) == 1) {
    if (genes == "all") {
      genes <- colnames(data)
    } else {
      stop("'genes' should be 'all' or a vector.")
    }
  } else {
    genes.404 <- setdiff(genes, colnames(data))
    if (length(genes.404) > 0) {
      warning(paste(length(genes.404), "gene(s) are not found in 'data'."))
      genes <- setdiff(genes, genes.404)
    }
  }
  cores <- ifelse(cores > 0, cores, parallel::detectCores())
  ## prepare data for VD
  vd.vars.str <- sapply(vd.vars, function(xx) sprintf("(1|%s)", xx))
  modelFormulaStr <- paste("expression ~", paste(vd.vars.str, collapse = "+"))
  data.use <- cbind(data[, genes], meta.data)
  ## exe VD
  vd.res <- do.call(rbind, parallel::mclapply(genes, function(genename) {
    data.model <- data.use[, c(vd.vars, genename)]
    colnames(data.model) <- c(vd.vars, "expression")
    tryCatch({
      model <- suppressWarnings(lme4::lmer(stats::as.formula(modelFormulaStr), data = data.model, REML = TRUE, verbose = FALSE))
      results <- as.data.frame(lme4::VarCorr(model))
      rownames(results) <- results$grp
      results <- results[c(vd.vars, "Residual"), ]
      frac.var <- results$vcov / sum(results$vcov)
      
      res.tmp <- c("OK", frac.var)
    },
    error = function(e) {
      print(e)
      res.tmp <- c("FAIL", rep(-1, length(vd.vars)+1))
    })
    names(res.tmp) <- c("status", vd.vars, "residual")
    as.data.frame(as.list(res.tmp)) # return
  }, mc.cores = cores)
  )
  rownames(vd.res) <- genes
  vd.res<- vd.res  %>% as.data.frame()
  vd.res <- vd.res %>% dplyr::mutate(gene = rownames(vd.res), .before=1)
  # vd.res <- vd.res %>% as.data.frame() %>% dplyr::mutate(gene = rownames(.), .before=1)
  for (i in 3:ncol(vd.res)) {
    vd.res[[i]] <- vd.res[[i]]  %>% as.numeric()
  }
  return(vd.res)
}
DefaultAssay(seu.TAM) <- "AUCell"

vd.vars <- c("celltype_pl")
meta.data <- data.frame(celltype_pl=seu.TAM@meta.data[, vd.vars])
ras.data <- FetchData(seu.TAM, vars = rownames(seu.TAM))
vd.res <- VarDecompose(data = ras.data, meta.data = meta.data, vd.vars = vd.vars, cores = 5)
data.table::fwrite(vd.res,
                   './04Scenic/vd.res.txt',sep='\t',
                   quote = F,row.names = F)
## PCA和VD结果的一致性
vd.res$PC1.loadings <- Loadings(seu.TAM, reduction = "pcaRAS")[vd.res$gene, 1]
vd.res$PC2.loadings <- Loadings(seu.TAM, reduction = "pcaRAS")[vd.res$gene, 2]



##DAM directed TFs-------


calIndMat <- function(celltype.vector) {
  cell.types <- as.character(unique(celltype.vector))
  ctMat <- lapply(cell.types, function(ii) {
    as.numeric(celltype.vector == ii)
  }) %>% do.call(cbind, .)
  colnames(ctMat) <- cell.types
  rownames(ctMat) <- names(celltype.vector)
  return(ctMat)
}
calRSSMat <- function(rasMat, ctMat){
  rssMat <- pbapply::pblapply(colnames(rasMat), function(i) {
    sapply(colnames(ctMat), function(j) {
      suppressMessages(
        1 - philentropy::JSD(rbind(rasMat[, i], ctMat[, j]), unit = 'log2', est.prob = "empirical")
      )
    })
  }) %>% do.call(rbind, .)
  rownames(rssMat) <- colnames(rasMat)
  colnames(rssMat) <- colnames(ctMat)
  return(rssMat)
}

rasMat <- t(seu.TAM[["AUCell"]]@data)
dim(rasMat)
#计算cell type indicate matrix
ctMat <- calIndMat(seu.TAM$celltype_pl)
dim(ctMat)
#计算RSS matrix
rssMat <- calRSSMat(rasMat, ctMat)
dim(rssMat)
rssMat<-rssMat %>% as.data.frame()
rssMat$gene<-rownames(rssMat)
data.use<-rssMat %>% inner_join(vd.res)
to.label <- subset(data.use, DAM>0.4)

ggplot(data.use, aes(DAM, `DAM-like.Mac`, 
                     color = abs(celltype_pl))) +
  geom_point(size = 2) +
  #geom_hline(yintercept = 0.75, color = "blue", linetype="dashed") +
  geom_vline(xintercept = 0.9, color = "blue", linetype="dashed") +
  ggrepel::geom_text_repel(inherit.aes = F, data = to.label,
                           max.overlaps = Inf, 
                           aes(DAM, `DAM-like.Mac`,
                               label=gene)) +
  xlab("Fraction of variance by DAM") +
  ylab("Fraction of variance by DAM-like.Mac") +
  scale_color_viridis_c() +
  theme_bw(base_size = 15)


data.table::fwrite(rssMat,
                   './04Scenic/seu.TAM_rssMat.txt',sep='\t',
                   quote = F,row.names = F)
data.table::fwrite(data.use,
                   './04Scenic/TAM.celltype_regulon.txt',sep='\t',
                   quote = F,row.names = F)

## Key gene regulon Figure4G----

genes<-readRDS('~/Brain/brian/DAM/KEY_Gene_DAM_ALL.rds')
regulons <- read.gmt("./04Scenic/output/Lung2BRME.regulons.gmt")
e.regulon <- enricher(gene = genes, TERM2GENE = regulons, minGSSize = 10, maxGSSize = 5000)

pdf('./04Scenic/KEY_Gene_DAM_regulon_Lung2.pdf',height=4,width=6)
dotplot(e.regulon,showCategory = 50)
dev.off()
saveRDS(e.regulon@result,'./04Scenic/KEY_Gene_DAM_regulon_Lung2.rds')

genes<-readRDS('~/Brain/brian/DAM/Key_Gene_DAM_MDM_ALL.rds')
regulons <- read.gmt("./04Scenic/output/Lung2BRME.regulons.gmt")
e.regulon <- enricher(gene = genes, TERM2GENE = regulons, minGSSize = 10, maxGSSize = 5000)

pdf('./04Scenic/KEY_Gene_DAM-like.Mac_regulon_Lung2.pdf',height=4,width=6)
dotplot(e.regulon,showCategory = 50)
dev.off()
saveRDS(e.regulon@result,'./04Scenic/KEY_Gene_DAM-like.Mac_regulon_Lung2.rds')
DAM_eregulon<-readRDS('~/Brain/brian/DAM/DAM_regulon_byeregulon.rds')

enrichment_data<-e.regulon@result %>% 
  mutate(gene=sub("\\(.*", "", ID),
         count=as.numeric(sub("/.*", "", BgRatio))) %>% 
  filter(gene%in%DAM_eregulon) %>% 
  arrange(desc(Count))
enrichment_data$gene<-factor(enrichment_data$gene,levels=rev(enrichment_data$gene))
pdf('./Figuresubmit/DAM_eregulon.pdf',height=4,width=7)
ggplot(enrichment_data, 
       aes(x = gene, 
           y = Count)) +
  geom_point(aes(size = Count, 
                 color =-log10(p.adjust)) )+  # 用 size 映射 p 值
  scale_size(range = c(3, 10), name = "-log10(p.adjust)") +  
  scale_color_gradient(low = '#0099CC',high = '#CC3333')+
  labs(x = "", y = "Gene Count", title = "") +
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
## HALLMARK_TNFA_SIGNALING_VIA_NFKB" Figure4F-------
regulons <- clusterProfiler::read.gmt("./04Scenic/output/Lung2BRME.regulons.gmt")
hallmarks <- read.gmt("~/Brain/brian/Geneset/h.all.v2022.1.Hs.symbols.gmt")
genes <- subset(hallmarks, term == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")$gene
e.regulon <- enricher(gene = genes, TERM2GENE = regulons, minGSSize = 10, maxGSSize = 5000)
saveRDS(e.regulon@result,file='./04Scenic/TNF_TFs_l2.rds')
dotplot(e.regulon,showCategory = 30)

DAM_TNF<-readRDS('~/Brain/brian/DAM/DAM_TNF_regulon.rds')
DAM_TNF<-sub("\\(.*", "",DAM_TNF)
enrichment_data<-e.regulon@result %>% 
  mutate(gene=sub("\\(.*", "", ID),
         count=as.numeric(sub("/.*", "", BgRatio))) %>% 
  filter(gene%in%DAM_TNF) %>% 
  arrange(desc(Count))
enrichment_data$gene<-factor(enrichment_data$gene,levels=rev(enrichment_data$gene))
pdf('./Figuresubmit/DAM_TNF_regulon.pdf',height=4,width=7)
ggplot(enrichment_data, 
       aes(x = gene, 
           y = Count)) +
  geom_point(aes(size = Count, 
                 color =-log10(p.adjust)) )+  # 用 size 映射 p 值
  scale_size(range = c(3, 10), name = "-log10(p.adjust)") +  
  scale_color_gradient(low = '#0099CC',high = '#CC3333')+
  labs(x = "", y = "Gene Count", title = "") +
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
##DAM regulonFigure4I----
data<-data.table::fread('./04Scenic/TAM.celltype_regulon.txt')
data<-data%>% 
  arrange(desc(DAM)) %>% 
  mutate(Regulons=1:nrow(data),
         label=sapply(str_split(gene,'\\('),'[',1),
         RSS=DAM) %>% 
  dplyr::select(c('Regulons','RSS','label'))

DAM_regulon<-readRDS('~/Brain/brian/DAM/DAM_regulon_byrssMat_formgi.rds')
DAM_keyregulon<-readRDS('~/Brain/brian/DAM/DAM_regulon_byeregulon_regulon.rds')
DAM_keyregulon<-sapply(str_split(DAM_keyregulon,'\\('),'[',1)
data <- head(data, n=250)
data$regulate[data$label%in%DAM_regulon$gene]<-'DAM related regulon'
#data$regulate[data$label%in%DAM_keyregulon]<-'DAM key regulon'
data$regulate[!(data$label%in%DAM_regulon$gene)]<-'NoSig'

data.label<-subset(data,label%in%DAM_keyregulon)
table(data$regulate)


pdf('./Figuresubmit/DAM_regulon.pdf',height=5,width=6)
ggplot(data, aes(Regulons, RSS)) +
  geom_point(size=3,aes(color=regulate)) +
  ggrepel::geom_text_repel(inherit.aes = FALSE, max.overlaps = Inf,
                           data = data.label, 
                           aes(Regulons,RSS,label=label), 
                           color='red',
                           size=6, show.legend = FALSE) +
  scale_color_manual(values=c(`DAM related regulon`="#B294C7",
                              NoSig="#BECEE3"))+
  #geom_point(data=data.label,size=5,aes(color=regulate))+
  #ggtitle("Metc6") + 
  ylab("Specificity score") +
  theme_bw(base_size = 12) +
  theme(panel.grid  = element_blank(),
        aspect.ratio = 1,
        legend.text = element_text(size=12),  
        plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
        axis.text = element_text(color = "black",size=14),
        axis.title = element_text(color = "black",size=14),
        legend.position = "right",
        legend.title = element_blank())+
  labs(color = NULL) 
dev.off()

