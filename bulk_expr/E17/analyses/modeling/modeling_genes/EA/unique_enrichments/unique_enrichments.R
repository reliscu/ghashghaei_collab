setwd("/mnt/bdata/rebecca/collaborations/troy/bulk_expr/E17/analyses/modeling/modeling_genes/EA/unique_enrichments")

source("unique_enrichments_fxn.R")

projectname <- "featureCounts_E17"
all_pred <- c("genotype", "region", "interaction")
pval_cut <- .05

set_source <- "MO"

# score_type <- "neg"
# predictor <- "genotype"

for(predictor in all_pred){
  
  for(score_type in c("pos", "neg")){
    
    print(paste(predictor, toupper(score_type)))
    cat("\n")
    
    unique_enrichments(projectname, predictor, all_pred, set_source, score_type, pval_cut)
    
  } ## for(score_type in c("pos", "neg")){
  
} ## for(predictor in all_pred){
