rm(list=ls())
gc()
library(Seurat);library(tidyverse);library(Nebulosa);
library(clustree);library(GEOquery);library(scCustomize);library(scRNAtoolVis)
library(UCell);library(clusterProfiler)
seu<-qs::qread('./Lung/03addTNBC/seu.all.lung.breast.qs')
rownames(seu@meta.data)<-seu$sample

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


table(seu$type)
seu.list <- SplitObject(seu, split.by = "type")
seu.list <- lapply(seu.list, function(object) {
  object@project.name <- unique(object$type)
  return(object)
})
metacells.list <- lapply(seq_along(seu.list), function(ii) {
  makeMetaCells(
    seu       = seu.list[[ii]],
    min.cells = 10,
    reduction = "umap",
    dims      = 1:2,
    k.param   = 10,
    cores     = 20)
})
mc.mat <- lapply(metacells.list, function(mc) mc$mat) %>% Reduce(cbind, .)
mc.cellmeta <- lapply(metacells.list, function(mc) mc$metadata) %>% Reduce(rbind, .)
saveRDS(mc.mat, "./dataset6Scenic/input/dataset6.luad.tnbc.mat.rds")
## (1) TF list文件
motif2tfs <- data.table::fread("/Users/huzixin/Brain/brian/dataset2/dataset2Scenic/motifs_v10nr_clust_nr.hgnc_m0.001_o0.0.tbl")
TFs <- sort(unique(motif2tfs$gene_name))
writeLines(TFs, "./dataset6Scenic/cisTarget_db_hsa_hgnc_tfs.txt")

## (2) meta cell matrix (for step1): *.csv or *.loom
mc.mat <- readRDS("./dataset6Scenic/input/dataset6.luad.tnbc.mat.rds")
expr.in.cells <- rowSums(mc.mat > 0)
mc.mat <- mc.mat[expr.in.cells >= 5, ]
cisdb <- arrow::read_feather("../dataset3/dataset3Scenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
dim(cisdb)
genes.use <- intersect(colnames(cisdb), rownames(mc.mat))
length(genes.use)
dim(mc.mat)
mc.mat <- mc.mat[genes.use, ]
dim(mc.mat)


dir.create("./dataset6Scenic/input")
loom <- SCopeLoomR::build_loom(
  file.name         = "./dataset6Scenic/input/dataset6_mat_for_step1.loom",
  dgem              = mc.mat,
  default.embedding = NULL
)
loom$close()
rm(loom)


#sh------
cd ./dataset6Scenic
conda activate pyscenic_env
## inputs
f_loom_grn=./input/dataset6_mat_for_step1.loom
## outputs
grn_output=./output/s1_adj.tsv
ctx_output=./output/s2_reg.tsv



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
regulons <- clusterProfiler::read.gmt("./dataset6Scenic/output/Dataset6BRME.regulons.gmt")

rg.names <- unique(regulons$term)
regulon.list <- lapply(rg.names, function(rg) {
  subset(regulons, term == rg)$gene
})
names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
summary(sapply(regulon.list, length))
saveRDS(regulon.list, "dataset6Scenic/output/regulons.rds")

seu.TAM<-readRDS('./Figuresubmit/seu.TAM.rds')
DefaultAssay(seu.TAM) <- "RNA"
seu.TAM <- ComputeModuleScore(seu.TAM, 
                              gene.sets = regulon.list, 
                              min.size = 10, 
                              cores = 10)
DefaultAssay(seu.TAM)<-'AUCell'
head(rownames(seu.TAM))


qs::qsave(seu.TAM, "./dataset6Scenic/seu.TAM_addScenic.qs")

##Variance Deposition-------
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
table(seu.TAM$type)
vd.vars <- c("celltype_pl", "type")
meta.data <- seu.TAM@meta.data[, vd.vars]
ras.data <- FetchData(seu.TAM, vars = rownames(seu.TAM))
vd.res <- VarDecompose(data = ras.data, meta.data = meta.data, vd.vars = vd.vars, cores = 5)
data.table::fwrite(vd.res,
                   './dataset6Scenic/vd.res.txt',sep='\t',
                   quote = F,row.names = F)
## PCA和VD结果的一致性
vd.res$PC1.loadings <- Loadings(seu.TAM, reduction = "pcaRAS")[vd.res$gene, 1]
vd.res$PC2.loadings <- Loadings(seu.TAM, reduction = "pcaRAS")[vd.res$gene, 2]

to.label <- subset(vd.res, celltype_pl > 0.25 & type > 0.10)
to.label <- subset(vd.res, celltype_pl > 0.8)

ggplot(vd.res, aes(celltype_pl, type, color = abs(PC2.loadings))) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.10, color = "blue", linetype="dashed") +
  geom_vline(xintercept = 0.25, color = "blue", linetype="dashed") +
  ggrepel::geom_text_repel(inherit.aes = F, 
                           data = to.label,
                           max.overlaps = Inf, 
                           aes(celltype_pl, type,label=gene)) +
  xlab("Fraction of variance by TAM cluster") +
  ylab("Fraction of variance by type") +
  scale_color_viridis_c() +
  theme_bw(base_size = 15)
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


## Key gene regulon Figure4G----
seu.TAM<-qs::qread( "./dataset6Scenic/seu.TAM_addScenic.qs")
genes<-readRDS('../DAM/KEY_Gene_DAM_ALL.rds')
regulons <- read.gmt("./dataset6Scenic/output/Dataset6BRME.regulons.gmt")
e.regulon <- enricher(gene = genes, TERM2GENE = regulons, minGSSize = 10, maxGSSize = 5000)

pdf('./dataset6Scenic/KEY_Gene_DAM_regulon_dataset6.pdf',height=4,width=6)
dotplot(e.regulon,showCategory = 50)
dev.off()
saveRDS(e.regulon@result,'./dataset6Scenic/KEY_Gene_DAM_regulon_dataset1.rds')

genes<-readRDS('../DAM/Key_Gene_DAM_MDM_ALL.rds')
regulons <- read.gmt("./dataset6Scenic/output/Dataset6BRME.regulons.gmt")
e.regulon <- enricher(gene = genes, TERM2GENE = regulons, minGSSize = 10, maxGSSize = 5000)

pdf('./dataset6Scenic/KEY_Gene_DAM-like.Mac_regulon_dataset6.pdf',height=4,width=6)
dotplot(e.regulon,showCategory = 50)
dev.off()
saveRDS(e.regulon@result,'./dataset6Scenic/KEY_Gene_DAM-like.Mac_regulon_dataset6.rds')

DAM_eregulon<-readRDS('./DAM/DAM_regulon_byeregulon.rds')

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
##HALLMARK_TNFA_SIGNALING_VIA_NFKB Figure4F-------
hallmarks <- read.gmt("../Geneset/h.all.v2022.1.Hs.symbols.gmt")
genes <- subset(hallmarks, term == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")$gene
e.regulon <- enricher(gene = genes, TERM2GENE = regulons, minGSSize = 10, maxGSSize = 5000)
saveRDS(e.regulon@result,file='./dataset6Scenic/TNF_TFs_d6.rds')
dotplot(e.regulon,showCategory = 30)
DAM_TNF<-readRDS('../DAM/DAM_TNF_regulon.rds')
DAM_TNF<-sub("\\(.*", "",DAM_TNF)
enrichment_data<-e.regulon@result %>% 
  mutate(gene=sub("\\(.*", "", ID),
         count=as.numeric(sub("/.*", "", BgRatio))) %>% 
  filter(gene%in%DAM_TNF) %>% 
  arrange(desc(Count))
enrichment_data$gene<-factor(enrichment_data$gene,levels=rev(enrichment_data$gene))
pdf('./Figu Figure4Fresubmit/DAM_TNF_regulon.pdf',height=4,width=7)
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
##DAM regulon Figure4I-----
data<-data.table::fread('./dataset6Scenic/TAM.celltype_regulon.txt')
data<-data%>% 
  arrange(desc(DAM)) %>% 
  mutate(Regulons=1:nrow(data),
         label=sapply(str_split(gene,'\\('),'[',1),
         RSS=DAM) %>% 
  dplyr::select(c('Regulons','RSS','label'))

DAM_regulon<-readRDS('../DAM/DAM_regulon_byrssMat_formgi.rds')
DAM_keyregulon<-readRDS('../DAM/DAM_regulon_byeregulon_regulon.rds')
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

