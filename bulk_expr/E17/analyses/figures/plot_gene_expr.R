setwd("/mnt/bdata/rebecca/collaborations/troy/bulk_expr/E17/analyses/figures")

source("plot_gene_expr_fxns.R")

load("/home/rebecca/gene_sets/MO/MO_sets_mapped.RData")

projectname <- "featureCounts_E17"

sampleinfo <- read.csv("../../sampleinfo_bulk_E17.csv")
sampleinfo$Genotype <- sapply(strsplit(sampleinfo$Genotype, " "), "[", 1)
sampleinfo <- sampleinfo %>% arrange(Genotype, Region)

meta_cols <- 1:2
gene_name_col <- 2
expr <- fread("../SN/featureCounts_SampleNetworks/All_03-16-04/featureCounts_All_12_ComBat.csv", data.table=F)
expr <- expr[,c(meta_cols, match(make.names(sampleinfo$Label), colnames(expr)))]
all.equal(colnames(expr)[-c(meta_cols)], make.names(sampleinfo$Label))
# [1] TRUE

genotype_vec <- sampleinfo$Genotype
region_vec <- sampleinfo$Region

genotype_vec <- gsub("Control", "Ctrl", genotype_vec)
genotype_vec <- gsub("NcreEGFRFF", "Nes:F/F", genotype_vec)
region_vec <- gsub("Rostral", "Rost", region_vec)
region_vec <- gsub("Caudal", "Caud", region_vec)

model_df <- fread("../modeling/modeling_genes/featureCounts_E17_modeling_37545_genes_scaled.csv", data.table=F)

enrich_files <- list.files(path="../modeling/modeling_genes/EA/", pattern="MO", full.names=T)

############################################# Genotype NEG #############################################

predictor <- "Genotype"
score_type <- "neg"

enrich_idx <- intersect(grep(predictor, enrich_files, ignore.case=T), grep(score_type, enrich_files, ignore.case=T))
enrich <- fread(enrich_files[enrich_idx], data.table=F)

## Oligodendrocytes:

enrichment <- "oligodendrocyte"

set_id <- MO_legend$SetID[is.element(MO_legend$SetName, c("DOYLE_CORTICAL_OLIGODENDROCYTES_OLIG2", "DOYLE_CEREBELLAR_OLIGODENDROCYTES_OLIG2", "DOYLE_CORTICAL_OLIGODENDROCYTES_CMTM5", "kelley_oligo"))]
working_set <- unique(unname(unlist(MO_sets_mapped[is.element(names(MO_sets_mapped), set_id)])))
model_set_genes <- model_df[is.element(toupper(model_df$Gene), working_set),] %>%
  arrange(Genotype_Est) %>%
  filter(Genotype_Pval<.05)

genes <- model_set_genes$Gene
#genes <- c("S100b", "Mbp", "Cnp", "Enpp6", "Gal3st1", "Plp1", "Mag", "Abcb9")
plot_gene_expr_violin(projectname, expr, meta_cols, gene_name_col, genes, enrichment, predictor, score_type)

## Egfr:

enrichment <- "Egfr"
genes <- "Egfr"
plot_gene_expr_violin(projectname, expr, meta_cols, gene_name_col, genes, enrichment, predictor, score_type)

## Olig1/2:

enrichment <- "Olig1_2"
genes <- c("Olig1", "Olig2")
plot_gene_expr_violin(projectname, expr, meta_cols, gene_name_col, genes, enrichment, predictor, score_type)

############################################# Interaction NEG #############################################

predictor <- "Interaction"
score_type <- "neg"

enrich_idx <- intersect(grep(predictor, enrich_files, ignore.case=T), grep(score_type, enrich_files, ignore.case=T))
enrich <- fread(enrich_files[enrich_idx], data.table=F)

## Oligodendrocytes:

enrichment <- "oligodendrocyte"

set_id <- MO_legend$SetID[is.element(MO_legend$SetName, c("MORABITO_2021_OLIGO12_DE_SUBTYPE_CLUSTERS", "NAGY_2020_OLIG3", "MORABITO_2021_OLIGO9_DE_SUBTYPE_CLUSTERS", "ABSINTA_2021_MS_OLIGODENDROCYTE_DE_OG_CLUSTERS"))]
working_set <- unique(unname(unlist(MO_sets_mapped[is.element(names(MO_sets_mapped), set_id)])))
model_set_genes <- model_df[is.element(toupper(model_df$Gene), working_set),] %>%
  arrange(Genotype_Est) %>%
  filter(Genotype_Pval<.05)

genes <- c("S100b", "Mbp", "Cnp", "Enpp6", "Gal3st1", "Plp1", "Mag", "Abcb9")
plot_gene_expr_violin(projectname, expr, meta_cols, gene_name_col, genes, enrichment, predictor, score_type)

