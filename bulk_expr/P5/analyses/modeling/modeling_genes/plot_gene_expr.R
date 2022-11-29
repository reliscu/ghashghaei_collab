setwd("/mnt/bdata/rebecca/collaborations/troy/bulk_expr/P5/analyses/modeling/modeling_genes")

source("plot_gene_expr_fxns.R")

projectname <- "featureCounts_P5"

sampleinfo <- read.csv("../../../sampleinfo_bulk_P5.csv")
sampleinfo$Genotype <- sapply(strsplit(sampleinfo$Genotype, " "), "[", 1)
sampleinfo <- sampleinfo %>% arrange(Genotype, Region)

meta_cols <- 1:2
gene_name_col <- 2
expr <- fread("../../SN/featureCounts_P5_SampleNetworks/All_06-26-21/featureCounts_P5_All_12_outliers_removed.csv", data.table=F)
expr <- expr[,c(meta_cols, match(make.names(sampleinfo$Plot_Label), colnames(expr)))]
all.equal(colnames(expr)[-c(meta_cols)], make.names(sampleinfo$Plot_Label))
# [1] TRUE

genotype_vec <- sampleinfo$Genotype
region_vec <- sampleinfo$Region

genotype_vec <- gsub("Control", "Ctrl", genotype_vec)
region_vec <- gsub("Rostral", "Rost", region_vec)
region_vec <- gsub("Caudal", "Caud", region_vec)

percentile_cut <- 70
model_df <- fread("featureCounts_P5_modeling_42863_genes_scaled.csv", data.table=F)

predictor <- "Genotype"
score_type <- "neg"

for(predictor in c("Genotype", "Region", "Interaction")){
  for(score_type in c("neg", "pos")){
    
    plot_gene_expr_violin(projectname, model_df, predictor, score_type, percentile_cut, gene_name, expr, meta_cols, gene_name_col, region_vec, genotype_vec)
    plot_gene_expr_barplot(projectname, model_df, predictor, score_type, percentile_cut, gene_name, expr, meta_cols, gene_name_col, region_vec, genotype_vec)
    
  }
}