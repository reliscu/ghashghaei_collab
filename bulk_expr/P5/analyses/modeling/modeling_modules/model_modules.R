setwd("/mnt/bdata/rebecca/collaborations/troy/bulk_expr/P5/analyses/modeling/modeling_modules")

source("model_modules_fxn.R")

fm_dir <- "/mnt/bdata/rebecca/collaborations/troy/bulk_expr/P5/analyses/FM/featureCounts_P5_Modules"

projectname <- "featureCounts_P5"

sampleinfo <- read.csv("../../../sampleinfo_bulk_P5.csv")

genotype <- sapply(strsplit(sampleinfo$Genotype, " "), "[", 1)
genotype_vec <- as.factor(genotype)
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
