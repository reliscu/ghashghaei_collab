setwd("/mnt/bdata/rebecca/collaborations/troy/bulk_expr/E17/analyses/modeling/modeling_genes")

source("model_genes_fxn.R")

projectname <- "featureCounts_E17"

sampleinfo <- read.csv("../../../sampleinfo_bulk_E17.csv")
expr <- fread("../../SN/featureCounts_SampleNetworks/All_03-16-04/featureCounts_All_12_ComBat.csv", data.table=F)
meta_cols <- 1:2
gene_name_col <- 2
expr <- expr[,c(meta_cols, match(make.names(sampleinfo$Label), colnames(expr)))]
all.equal(colnames(expr)[-c(meta_cols)], make.names(sampleinfo$Label))
# [1] TRUE

expr <- expr[!grepl("^mt-", expr$Gene, ignore.case=T),]

genotype_vec <- as.factor(sampleinfo$Genotyp)
levels(genotype_vec)
# [1] "Control"    "NcreEGFRFF"

region_vec <- as.factor(sampleinfo$Region)
levels(region_vec)
# [1] "Caudal"  "Rostral"

scale <- T

model_genes(
  projectname, 
  expr, 
  meta_cols, 
  gene_name_col, 
  genotype_vec, 
  region_vec, 
  scale
)
