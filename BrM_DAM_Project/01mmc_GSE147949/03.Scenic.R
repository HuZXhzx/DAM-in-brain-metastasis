rm(list=ls())
library(Seurat);library(tidyverse);library(Nebulosa);library(clustree);library(GEOquery)
library(scCustomize);library(clusterProfiler);library(AUCell)
## pre-Scenic-----
seu<-qs::qread('./dataset3/tmp/02-4.seurat_mmc_anno.seurat.qs')
DefaultAssay(seu)<-"integrated"
DimPlot(seu, group.by = "seurat_clusters")
table(seu$type)
seu.list <- SplitObject(seu, split.by = "type")
seu.list <- lapply(seu.list, function(object) {
  object@project.name <- unique(object$type)
  return(object)
})
makeMetaCells <- function(seu, min.cells=10, reduction="umap", dims=1:2, 
                          k.param=10, cores=20) {
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
summary(mc.cellmeta$CELL_COUNT)
seu2 <- CreateSeuratObject(mc.mat)
seu2 <- NormalizeData(seu2)
FeatureScatter(seu, feature1 = "Cd8a", feature2 = "Gzmb") +
  FeatureScatter(seu2, feature1 = "Cd8a", feature2 = "Gzmb")

saveRDS(mc.mat, "./dataset3/tmp/02-5.seurat_mmc_scenic_mc.mat_v2.rds")
## (1) TF list文件
motif2tfs <- data.table::fread("dataset3Scenic/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
TFs <- sort(unique(motif2tfs$gene_name))
writeLines(TFs, "./dataset3Scenic/cisTarget_db_mm10_tfs.txt")

## (2) meta cell matrix (for step1): *.csv or *.loom
mc.mat <- readRDS("./dataset3/tmp/02-5.seurat_mmc_scenic_mc.mat_v2.rds")
dim(mc.mat)
## 过滤低表达基因
expr.in.cells <- rowSums(mc.mat > 0)
mc.mat <- mc.mat[expr.in.cells >= 5, ]
dim(mc.mat)
## 过滤不在cisTargetDB中的基因
cisdb <- arrow::read_feather("./dataset3Scenic/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
dim(cisdb)
#cisdb <- arrow::read_feather("./dataset3Scenic/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")


genes.use <- intersect(colnames(cisdb), rownames(mc.mat))
length(genes.use)
dim(mc.mat)
mc.mat <- mc.mat[genes.use, ]
dim(mc.mat)


dir.create("./dataset3Scenic/input")
loom <- SCopeLoomR::build_loom(
  file.name         = "./dataset3Scenic/input/mc_mat_for_step1.loom",
  dgem              = mc.mat,
  default.embedding = NULL
)
loom$close()
rm(loom)
## AfterScenic --------------
seu<-qs::qread('./dataset3/tmp/02-4.seurat_mmc_anno.seurat.qs')
table(subset(seu,celltype_cluster==23)$combined_cell_type)
table(subset(seu,seurat_clusters==25)$combined_cell_type)

table(seu$type)
table(seu$combined_cell_type)
regulons <- clusterProfiler::read.gmt("./dataset3Scenic/output/02_mmc_BRBMET.regulons.gmt")

rg.names <- unique(regulons$term)
regulon.list <- lapply(rg.names, function(rg) {
  subset(regulons, term == rg)$gene
})
names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
summary(sapply(regulon.list, length))
saveRDS(regulon.list, "dataset3Scenic/output/03-1.mmc.regulons.rds")

seu.TAM<-readRDS('./Figuresubmit/seu.TAM.rds')

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
regulon.list<-qs::qread( "dataset3Scenic/output/03-1.mmc.regulons.rds")
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


##Variance deposition-------
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
vd.vars <- c("TAM_cluster", "type")
meta.data <- seu.TAM@meta.data[, vd.vars]
ras.data <- FetchData(seu.TAM, vars = rownames(seu.TAM))
vd.res <- VarDecompose(data = ras.data, meta.data = meta.data, vd.vars = vd.vars, cores = 5)

##MET related TF-------
DefaultAssay(seu.TAM) <- "AUCell"
seu.TAM$type <- factor(seu.TAM$type, levels = c("Met","Con"))
seu.TAM$new.group <- paste(seu.TAM[['TAM_cluster', drop=T]], 
                           seu.TAM[['type', drop=T]], sep = "_")
Idents(seu.TAM) <- factor(seu.TAM$new.group)

celltypes <- setdiff(unique(seu.TAM[['TAM_cluster', drop=T]]),'c11')
groups <- levels(seu.TAM[['type', drop=T]])

de.list <- lapply(celltypes, function(ct) {
  message(glue::glue("processing {ct} ..."))
  ct1 <- paste(ct, groups[1], sep = "_")
  ct2 <- paste(ct, groups[2], sep = "_")
  de <- FindMarkers(seu.TAM, ident.1 = ct1, ident.2 = ct2, test.use = "wilcox", fc.name = "avg_diff", logfc.threshold = 0)
  de$change <- ifelse(de$avg_diff > 0,
                      paste("up in", groups[1]),
                      paste("up in", groups[2]))
  
  #de$diff_rate <- de$avg_diff / baseline.levels[rownames(de), ct]
  de$group <- ct
  return(de)
})
names(de.list) <- celltypes

de.regulons <- lapply(de.list) %>%
  Reduce(rbind, .) %>%
  pivot_wider(names_from = "group", values_from = "change")


de.regulons <- de.regulons %>% dplyr::select(c("gene", all_of(names(de.list))))
de.regulons[is.na(de.regulons)] <- ""


##celltype TFs----

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
PlotRegulonRank <- function(rssMat, cell.type, topn=5) {
  data <- data.frame(
    Regulons = 1:nrow(rssMat),
    RSS = sort(rssMat[, cell.type], decreasing = T),
    label = sub("(+)", "", names(sort(rssMat[, cell.type], decreasing = T)), fixed = T)
  )
  
  data$pt.col <- ifelse(data$Regulons <= topn, "#007D9B", "#BECEE3")
  data <- head(data, n=200)
  data.label <- head(data, n=topn)
  
  ggplot(data, aes(Regulons, RSS)) +
    geom_point(size=3, color=data$pt.col) +
    ggrepel::geom_text_repel(inherit.aes = FALSE, data = data.label, aes(Regulons, RSS, label=label), size=4) +
    ggtitle(cell.type) + ylab("Specificity score") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = .5)
    )
}
rasMat <- t(seu.TAM[["AUCell"]]@data)
dim(rasMat)
ctMat <- calIndMat(seu.TAM$TAMcluster_type)
dim(ctMat)
rssMat <- calRSSMat(rasMat, ctMat)
dim(rssMat)
data<-rssMat[,c('Metc6','gene')] %>% arrange(desc(Metc6)) %>% 
  mutate(Regulons=1:nrow(rssMat),
         label=sub("(+)", "",gene),
         RSS=Metc6) %>% 
  dplyr::select(c('Regulons','RSS','label'))
data$label<-sapply(str_split(data$label,'\\('),'[',1)

DAM_regulon<-readRDS('./DAM/DAM_regulon_byrssMat_formgi.rds')
data.use<- vd.res%>% 
  inner_join(rssMat [,c('Metc6','gene')]) %>% 
  inner_join(de.regulons[de.regulons$group=='c6',]) 

data.label<-subset(data.use,gene%in%DAM_regulon$Mgi_TF)

pdf('./V2/Scenic/DAM_regulon_v2.pdf',height=8,width=8)
ggplot() +
  geom_point(data = data.use, 
             aes(x = TAM_cluster, y = type), 
             color = "#BECEE3", size = 2) +
  geom_point(data = data.label, 
             aes(x = TAM_cluster, y = type, color = change), 
             size = 3) +
  ggrepel::geom_text_repel(data = data.label, 
                           aes(x = TAM_cluster, y = type, label = gene, color = change), 
                           size = 6, 
                           inherit.aes = FALSE, 
                           max.overlaps = Inf, 
                           show.legend = FALSE) +
  scale_color_manual(values = c(`up in Met` = '#CC3333', 
                                `up in Con` = '#0099CC')) +
  # geom_hline(yintercept = 0.10, color = "blue", linetype = "dashed") +
  # geom_vline(xintercept = 0.25, color = "blue", linetype = "dashed") +
  xlab("Fraction of variance by TAM cluster") +
  ylab("Fraction of variance by Group") +
  theme_bw(base_size = 15)+
  theme(panel.grid  = element_blank(),
        aspect.ratio = 1,
        legend.text = element_text(size=12),  
        plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
        axis.text = element_text(color = "black",size=14),
        axis.title = element_text(color = "black",size=14),
        legend.position = "right",
        legend.title = element_blank())

dev.off()
##Enrichment of DAM ------
DefaultAssay(seu.TAM)<-'RNA'
Idents(seu.TAM)<-'TAMcluster_type'
TAM_markers<-data.table::fread('./V2/02markers_TAM_clusters_type.txt')%>%
  filter(cluster=='c6'&group=='Met')

geneList <- TAM_markers$avg_log2FC
names(geneList) =  TAM_markers$gene

geneList = sort(geneList, decreasing = TRUE)
head(geneList)
library(clusterProfiler)
library(enrichplot)
hallmarks <- read.gmt("~/Brain/brian/Geneset/mh.all.v2023.2.Mm.symbols.gmt.txt")
y <- GSEA(geneList, TERM2GENE = hallmarks)
yd <- as.data.frame(y)


pdf('./Figuresubmit/HALLMARK_TNFA_SIGNALING_VIA_NFKB_gsea.pdf',height=8,width=10)
gseaplot2(y, "HALLMARK_TNFA_SIGNALING_VIA_NFKB",color = "red", pvalue_table = T)
dev.off()


##TNF----
hallmarks <- read.gmt("~/Brain/brian/Geneset/mh.all.v2023.2.Mm.symbols.gmt.txt")

genes <- subset(hallmarks, term == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")$gene
e.regulon <- enricher(gene = genes, TERM2GENE = regulons, minGSSize = 10, maxGSSize = 5000)
saveRDS(e.regulon@result,file='./V2/TNF_TFs_dataset3.rds')
DAM_TNF<-readRDS('../DAM/DAM_TNF_regulon.rds')

library(biomaRt)
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

DAM_TNF<- getLDS(attributes = c("external_gene_name"),
                 filters = "external_gene_name",
                 values = readRDS('../DAM/DAM_TNF_regulon.rds'),
                 mart = human,
                 attributesL = c("mgi_symbol"),
                 martL = mouse,uniqueRows = T
)

DAM_TNF<-sub("\\(.*", "",DAM_TNF$MGI.symbol)
saveRDS(DAM_TNF,'../DAM/DAM_TNF_regulonformgi.rds')
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
    
    legend.text = element_text(size=12),  
    plot.margin = margin(t = 1, r = 1,b = 1, l = 1,unit = "cm"),
    axis.text = element_text(color = "black",size=14),
    axis.title = element_text(color = "black",size=14),
    legend.position = "right",
    legend.title = element_blank())
dev.off()

##KEY GENE regulon----
library(biomaRt)
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
genes<-readRDS('../DAM/KEY_Gene_DAM_ALL.rds')
gene<- getLDS(attributes = c("external_gene_name"),
              filters = "external_gene_name",
              values = genes,
              mart = human,
              attributesL = c("mgi_symbol"),
              martL = mouse,uniqueRows = T
)
saveRDS(gene,'../DAM/KEY_Gene_DAM_ALL_MGI.rds')
gene<-readRDS('../DAM/KEY_Gene_DAM_ALL_MGI.rds')
regulons <- read.gmt("./dataset3Scenic/output/02_mmc_BRBMET.regulons.gmt")
e.regulon <- enricher(gene = gene$MGI.symbol, TERM2GENE = regulons, minGSSize = 10, maxGSSize = 5000)
DAM_eregulon<-readRDS('../DAM/DAM_regulon_byeregulon.rds')
DAM_eregulon<- getLDS(attributes = c("external_gene_name"),
                      filters = "external_gene_name",
                      values = DAM_eregulon,
                      mart = human,
                      attributesL = c("mgi_symbol"),
                      martL = mouse,uniqueRows = T
)


pdf('./V1/Keygene_regulon_dataset3.pdf.pdf',height=4,width=6)
dotplot(e.regulon,showCategory = 50)
dev.off()
saveRDS(e.regulon@result,'./V2/KEY_Gene_DAM_regulon_dataset3.rds')

enrichment_data<-e.regulon@result %>% 
  mutate(gene=sub("\\(.*", "", ID),
         count=as.numeric(sub("/.*", "", BgRatio))) %>% 
  filter(gene%in%DAM_eregulon$MGI.symbol) %>% 
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

