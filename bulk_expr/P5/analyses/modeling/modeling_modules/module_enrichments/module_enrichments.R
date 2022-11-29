setwd("/mnt/bdata/rebecca/collaborations/troy/bulk_expr/P5/analyses/modeling/modeling_modules/module_enrichments")

source("module_enrichments_fxn.R")

projectname <- "featureCounts_P5"
model_df <- fread("../featureCounts_P5_modeling_3205_modules.csv", data.table=F)
fm_dir <- "/mnt/bdata/rebecca/collaborations/troy/bulk_expr/P5/analyses/FM/featureCounts_P5_Modules"

module_enrichments(
  projectname,
  model_df,
  fm_dir
)