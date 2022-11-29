setwd("/mnt/bdata/rebecca/collaborations/troy/bulk_expr/E17/analyses/modeling/modeling_modules/module_enrichments")

source("module_enrichments_fxn.R")

projectname <- "featureCounts_E17"
model_df <- fread("../featureCounts_E17_modeling_1694_modules.csv", data.table=F)
fm_dir <- "/mnt/bdata/rebecca/collaborations/troy/bulk_expr/E17/analyses/FM/featureCounts_Modules"

module_enrichments(
  projectname,
  model_df,
  fm_dir
)