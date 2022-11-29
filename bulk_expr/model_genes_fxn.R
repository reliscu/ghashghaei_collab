library(dplyr)
library(data.table)
library(parallel)

source("/home/rebecca/code/misc/rank_percentile.R")

model_genes <- function(projectname, expr, meta_cols, gene_name_col, genotype_vec, region_vec, scaled=F){
  
  expr1 <- expr[,-meta_cols]
  mean_expr <- rowSums(expr1)
  mean_expr_percentile <- rank_percentile(mean_expr) 
  
  if(scaled){
    expr1 <- t(scale(t(expr1)))
  }
  
  interaction_model <- mclapply(1:nrow(expr), function(x){
    lm(as.numeric(expr1[x,]) ~ genotype_vec + region_vec + genotype_vec:region_vec)
  }, mc.cores=15)
  
  basic_model <- mclapply(1:nrow(expr), function(x){
    lm(as.numeric(expr1[x,]) ~ genotype_vec + region_vec)
  }, mc.cores=15)
  
  genotype_est <- sapply(interaction_model, function(x) summary(x)$coeff[2,1])
  genotype_pval <- sapply(interaction_model, function(x) summary(x)$coeff[2,4])
  region_est <- sapply(interaction_model, function(x) summary(x)$coeff[3,1])
  region_pval <- sapply(interaction_model, function(x) summary(x)$coeff[3,4])
  interaction_est <- sapply(interaction_model, function(x) summary(x)$coeff[4,1])
  interaction_pval <- sapply(interaction_model, function(x) summary(x)$coeff[4,4])
  basic_mdl_r2 <- sapply(basic_model, function(x) summary(x)$adj.r.squared)
  interaction_mdl_r2 <- sapply(interaction_model, function(x) summary(x)$adj.r.squared)
  
  model_df <- data.frame(
    Gene=expr[,gene_name_col],
    Genotype_Est=genotype_est,
    Genotype_Pval=genotype_pval,
    Region_Est=region_est,
    Region_Pval=region_pval,
    Interaction_Est=interaction_est,
    Interaction_Pval=interaction_pval,
    Basic_Model_R2=basic_mdl_r2,
    Interaction_Model_R2=interaction_mdl_r2,
    Mean_Expr=mean_expr,
    Mean_Expr_Percentile=mean_expr_percentile
  )
  
  file_name <- paste(projectname, "modeling", nrow(expr), "genes.csv", sep="_")
  
  if(scaled){
    file_name <- gsub(".csv", "_scaled.csv", file_name)
  }
  
  fwrite(model_df, file=file_name)
  
}

