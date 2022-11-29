setwd("/mnt/bdata/rebecca/collaborations/troy/bulk_expr/P5/analyses/modeling/modeling_genes")

source("model_genes_fxn.R")

projectname <- "featureCounts_P5"

sampleinfo <- read.csv("../../../sampleinfo_bulk_P5.csv")
sampleinfo$Genotype <- sapply(strsplit(sampleinfo$Genotype, " "), "[", 1)
sampleinfo <- sampleinfo %>% arrange(Genotype, Region)

expr <- fread("../../SN/featureCounts_P5_SampleNetworks/All_06-26-21/featureCounts_P5_All_12_outliers_removed.csv", data.table=F)
meta_cols <- 1:2
gene_name_col <- 2
expr <- expr[,c(meta_cols, match(make.names(sampleinfo$Plot_Label), colnames(expr)))]
all.equal(colnames(expr)[-c(meta_cols)], make.names(sampleinfo$Plot_Label))
# [1] TRUE

expr <- expr[!grepl("^mt-", expr$Gene, ignore.case=T),]

genotype_vec <- as.factor(sampleinfo$Genotype)
levels(genotype_vec)
# [1] "Control" "Nes:F/F"

region_vec <- as.factor(sampleinfo$Region)
levels(region_vec)
# [1] "Caudal"  "Rostral"

scale <- T

##

gene_name <- "Plch1"
x <- as.numeric(expr[is.element(expr$Gene, gene_name), -c(1:2)])
names(x) <- sampleinfo$Plot_Label
barplot(x, breaks=20)

model_genes(
  projectname, 
  expr, 
  meta_cols, 
  gene_name_col, 
  genotype_vec, 
  region_vec, 
  scale
)
