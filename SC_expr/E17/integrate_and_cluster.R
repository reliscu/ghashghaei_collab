setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/E17/cluster")

library(Seurat)
library(data.table)
library(dplyr)

source("/home/rebecca/code/misc/upper_first.R")

expr <- fread("../counts/QC_counts/expression_QC.csv", data.table=F)
sampleinfo <- read.csv("../sampleinfo_SC_E17.csv")

## Prep metadata:

sample_ids <- sapply(strsplit(colnames(expr[-c(1)]), "_"), function(x) paste(x[-length(x)], collapse="_"))
cellinfo <- data.frame(Cell_ID=colnames(expr)[-c(1)], Label=sample_ids)
sampleinfo <- sampleinfo[!duplicated(sampleinfo$Label),]
cellinfo <- merge(cellinfo, sampleinfo, by="Label")
cellinfo <- cellinfo %>% dplyr::select(-c(Sequencing_Date, Depth, Number_of_Cells_Loaded))

genes <- expr[,c(1)]
expr <- as(expr[,-c(1)], "matrix")
rownames(expr) <- genes
colnames(expr) <- cellinfo$Cell_ID
rownames(cellinfo) <- cellinfo$Cell_ID

############################################# Cluster all #############################################

## Ref: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

expr <- CreateSeuratObject(counts=expr, meta.data=cellinfo, row.names=genes)
expr[["percent.mt"]] <- PercentageFeatureSet(expr, pattern="^mt-")
expr <- subset(expr, subset=percent.mt<10&nFeature_RNA<6000)

## Get genes for regressing out cell cycle effects:

s.genes <- sapply(tolower(cc.genes$s.genes), upper_first)
g2m.genes <- sapply(tolower(cc.genes$g2m.genes), upper_first)

# ## First look at data w/o integration:
# 
# expr1 <- expr
# expr1 <- SCTransform(expr1, vars.to.regress=c('percent.mt', 'nFeature_RNA', 'nCount_RNA'), vst.flavor="v2")
# expr1 <- CellCycleScoring(expr1, s.features=s.genes, g2m.features=g2m.genes, assay='SCT', set.ident=T)
# expr1 <- SCTransform(expr1, vars.to.regress=c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'), vst.flavor="v2")
# expr1 <- RunPCA(expr1, npcs=30)
# expr1 <- RunUMAP(expr1, reduction="pca", dims=1:20)
# 
# pdf("sct_cell_cycle_regressed_no_integration_UMAP.pdf", width=12, height=7)
# DimPlot(expr1, reduction="umap", pt.size=.3, group.by="Library_Prep_Date")
# DimPlot(expr1, reduction="umap", pt.size=.3, group.by="Litter")
# FeaturePlot(expr1, features="nFeature_RNA") 
# FeaturePlot(expr1, features="percent.mt")
# dev.off()
# 
# ## Clearly there are still batch effects w/o integration.
# 
# rm(expr1)

## Integrate samples:

## Ref: https://satijalab.org/seurat/articles/integration_introduction.html

split <- SplitObject(expr, split.by="Library_Prep_Date")
split <- lapply(X=split, FUN=function(x){
  x <- SCTransform(x, vars.to.regress=c('percent.mt', 'nFeature_RNA', 'nCount_RNA'), vst.flavor="v2")
  x <- CellCycleScoring(x, s.features=s.genes, g2m.features=g2m.genes, assay='SCT', set.ident=T)
  x <- SCTransform(x, vars.to.regress=c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'), vst.flavor="v2")
})
features <- SelectIntegrationFeatures(object.list=split, nfeatures=3000)
split <- PrepSCTIntegration(object.list=split, anchor.features=features)
anchors <- FindIntegrationAnchors(object.list=split, anchor.features=features, normalization.method="SCT")
expr_int <- IntegrateData(anchorset=anchors, normalization.method="SCT")
expr_int <- RunPCA(expr_int, npcs=30)
expr_int <- RunUMAP(expr_int, reduction="pca", dims=1:30)
expr_int <- FindNeighbors(expr_int, reduction="pca", dims=1:30)
res <- .8
expr_int <- FindClusters(expr_int, resolution=res)

pdf(paste0("sct_cell_cycle_regressed_UMAP_cluster_res_", res, ".pdf"), width=12, height=7)
DimPlot(expr_int, reduction="umap", pt.size=.3, label=T)
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Library_Prep_Date")
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Litter")
FeaturePlot(expr_int, features="nFeature_RNA") 
FeaturePlot(expr_int, features="percent.mt")
dev.off()

## Visualize canonical markers:

pdf("sct_cell_cycle_regressed_canonical_markers.pdf", width=12, height=7)

FeaturePlot(
  object=expr_int, 
  features=c("Tmem119", "Cd11b", "Cd45"),
  label=T,
  repel=F,
  label.size=2,
  slot="data",
  pt.size=.1,
  order=T
)

FeaturePlot(
  object=expr_int, 
  features=c("Olig1", "Olig2", "Ccnd1"),
  label=T,
  repel=F,
  label.size=2,
  slot="data",
  pt.size=.1,
  order=T
)

FeaturePlot(
  object=expr_int, 
  features=c("Mt3", "Ddah1", "Slc1a3"), # "Vim"
  label=T,
  repel=F,
  label.size=2,
  slot="data",
  pt.size=.1,
  order=T
)

FeaturePlot(
  object=expr_int, 
  features=c("Neurod2", "Neurod6", "Zbtb18"),
  label=T,
  repel=F,
  label.size=2,
  pt.size=.1,
  order=T
)

FeaturePlot(
  object=expr_int, 
  features=c("Gad2", "Nrxn3", "Dlx6os1"), #  "Sp9"
  label=T,
  repel=F,
  label.size=2,
  pt.size=.1,
  order=T
)

dev.off()

## Prepare to find cluster markers:

expr_int <- PrepSCTFindMarkers(expr_int)

## Consolidate clusters

expr_int$Cluster <- "NA"
expr_int$Cluster[is.element(expr_int$seurat_clusters, c(9, 20))] <- "ASC"
expr_int$Cluster[is.element(expr_int$seurat_clusters, 18)] <- "OG"
expr_int$Cluster[is.element(expr_int$seurat_clusters, c(0, 1, 2, 3, 4, 5, 13, 15, 16, 17, 19))] <- "EXC"
expr_int$Cluster[is.element(expr_int$seurat_clusters, c(6, 7, 8, 10, 11, 12, 14))] <- "INH"
Idents(expr_int) <- "Cluster"

markers <- FindAllMarkers(expr_int, assay="SCT")
markers %>% dplyr::group_by(cluster) %>% slice_head(n=5) %>% as.data.frame()

#     p_val avg_log2FC pct.1 pct.2 p_val_adj cluster    gene
# 1      0  1.1601795 0.927 0.179         0     EXC Neurod6
# 2      0  1.0112427 0.855 0.103         0     EXC Neurod2
# 3      0  0.9895516 0.950 0.375         0     EXC  Zbtb18
# 4      0  0.8565787 0.934 0.466         0     EXC   Ttc28
# 5      0  0.7999780 0.815 0.242         0     EXC     Id2
# 6      0  1.1634940 0.864 0.042         0     INH   Nrxn3
# 7      0  1.0610238 0.752 0.020         0     INH Dlx6os1
# 8      0  1.0069986 0.736 0.131         0     INH    Meg3
# 9      0  0.8954642 0.675 0.023         0     INH    Dlx1
# 10     0  0.8867920 0.678 0.016         0     INH     Sp9
# 11     0  1.7542715 0.993 0.413         0     ASC   Fabp7
# 12     0  1.5983275 0.976 0.116         0     ASC     Mt3
# 13     0  1.5115550 1.000 0.413         0     ASC     Dbi
# 14     0  1.4149106 0.996 0.259         0     ASC     Vim
# 15     0  1.3204331 0.969 0.086         0     ASC   Ddah1
# 16     0  1.2979099 0.829 0.007         0      OG   Olig1
# 17     0  1.2508378 0.991 0.434         0      OG   Fabp7
# 18     0  1.0865712 0.936 0.287         0      OG     Vim
# 19     0  1.0790730 0.832 0.015         0      OG   Olig2
# 20     0  1.0093051 0.978 0.434         0      OG     Dbi
# 21     0  1.1705530 0.464 0.015         0      NA    Hexb
# 22     0  1.1112164 0.580 0.003         0      NA Selenop
# 23     0  1.0279714 0.562 0.022         0      NA    Cyba
# 24     0  1.0230364 0.598 0.020         0      NA     B2m
# 25     0  0.9173006 0.411 0.002         0      NA    Atf3


saveRDS(expr_int, file="expr_integrated_annotated.RDS")

############################################# DE genes within cell types #############################################

## Ref: https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/pseudobulk-expression.html