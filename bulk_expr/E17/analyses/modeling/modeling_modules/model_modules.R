setwd("/mnt/bdata/rebecca/collaborations/troy/bulk_expr/E17/analyses/modeling/modeling_modules")

source("model_modules_fxn.R")

fm_dir <- "/mnt/bdata/rebecca/collaborations/troy/bulk_expr/E17/analyses/FM/featureCounts_Modules"

projectname <- "featureCounts_E17"

sampleinfo <- read.csv("../../../sampleinfo_bulk_E17.csv")

genotype_vec <- as.factor(sampleinfo$Genotype)
levels(genotype_vec)
# [1] "Control" "Nes:F/F"

region_vec <- as.factor(sampleinfo$Region)
levels(region_vec)
# [1] "Caudal"  "Rostral"

model_modules(
  projectname, 
  fm_dir,
  sampleinfo, 
  genotype_vec, 
  region_vec
)
