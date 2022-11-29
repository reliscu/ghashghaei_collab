setwd("/mnt/bdata/rebecca/collaborations/troy/bulk_expr/E17/analyses/modeling/modeling_genes/EA")

library(data.table)

source("EA_fxns.R")

projectname <- "featureCounts_E17"
min_size <- 5
max_size <- 500

model_df <- fread("../featureCounts_E17_modeling_37545_genes_scaled.csv", data.table=F)

set_source <- "MO"

if(set_source=="MO"){
  
  load("/home/rebecca/gene_sets/MO/MO_sets_mapped.RData")
  
  gene_sets <- MO_sets_mapped
  set_legend <- MO_legend
  
} else {
  
  load("/home/rebecca/gene_sets/broad/broad_sets_v7_mapped.RData")
  
  gene_sets <- broad_sets_mapped
  set_legend <- broad_legend
  
}

for(predictor in c("Genotype", "Region", "Interaction")){
  
  for(score_type in c("pos", "neg")){
    
    print(paste(predictor, toupper(score_type)))
    cat("\n")
    
    EA_fxn(
      projectname,
      model_df,
      gene_sets,
      set_legend,
      set_source,
      predictor,
      min_size, max_size,
      score_type
    )
    
  } ## for(score_type in c("pos", "neg")){
  
} ## for(predictor in c("Genotype", "Region", "Interaction")){
