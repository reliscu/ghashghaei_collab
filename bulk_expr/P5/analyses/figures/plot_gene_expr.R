setwd("/mnt/bdata/rebecca/collaborations/troy/bulk_expr/P5/analyses/figures")

source("plot_gene_expr_fxns.R")

load("/home/rebecca/gene_sets/MO/MO_sets_mapped.RData")

projectname <- "featureCounts_P5"

sampleinfo <- read.csv("../../sampleinfo_bulk_P5.csv")
sampleinfo$Genotype <- sapply(strsplit(sampleinfo$Genotype, " "), "[", 1)
sampleinfo <- sampleinfo %>% arrange(Genotype, Region)

meta_cols <- 1:2
gene_name_col <- 2
expr <- fread("../SN/featureCounts_P5_SampleNetworks/All_06-26-21/featureCounts_P5_All_12_outliers_removed.csv", data.table=F)
expr <- expr[,c(meta_cols, match(make.names(sampleinfo$Plot_Label), colnames(expr)))]
all.equal(colnames(expr)[-c(meta_cols)], make.names(sampleinfo$Plot_Label))
# [1] TRUE

genotype_vec <- sampleinfo$Genotype
region_vec <- sampleinfo$Region

genotype_vec <- gsub("Control", "Ctrl", genotype_vec)
region_vec <- gsub("Rostral", "Rost", region_vec)
region_vec <- gsub("Caudal", "Caud", region_vec)

model_df <- fread("../modeling/modeling_genes/featureCounts_P5_modeling_42863_genes_scaled.csv", data.table=F)

enrich_files <- list.files(path="../modeling/modeling_genes/EA/", pattern="MO", full.names=T)

############################################# Genotype NEG #############################################

predictor <- "Genotype"
score_type <- "neg"

enrich_idx <- intersect(grep(predictor, enrich_files, ignore.case=T), grep(score_type, enrich_files, ignore.case=T))
enrich <- fread(enrich_files[enrich_idx], data.table=F)

## Astrocytes:

enrichment <- "astrocyte"

set_id <- MO_legend$SetID[is.element(MO_legend$SetName, c("BARRES_ASTROCYTES", "ZEISEL_ASTROCYTE", "kelley_astro", "ZHANG_MOUSE_ASTROCYTES_2014"))]
working_set <- unique(unname(unlist(MO_sets_mapped[is.element(names(MO_sets_mapped), set_id)])))
model_set_genes <- model_df[is.element(toupper(model_df$Gene), working_set),] %>%
  arrange(Genotype_Est) %>%
  filter(Genotype_Pval<.05)

# genes <- model_set_genes$Gene
genes <- c("Ephx2", "Lcat", "Fgfr3", "Adora2b", "Ppp1r3d", "Papss2", "Baalc", "Tm7sf2", "Aldh1l1", "Slc7a10", "Fabp7", "Aqp4", "Acsbg1", "Fgfbp3", "Ednrb", "S100a16")
plot_gene_expr_violin(projectname, expr, meta_cols, gene_name_col, genes, enrichment, predictor, score_type)

## Oligodendrocytes:

enrichment <- "oligodendrocyte"

set_id <- MO_legend$SetID[is.element(MO_legend$SetName, c("JAKEL_2019_OLIGODENDROCYTE5_DE_OG_CLUSTERS", "JAKEL_2019_OLIGODENDROCYTE6_DE_OG_CLUSTERS", "JAKEL_2019_OLIGODENDROCYTE5", "ZHANG_MOUSE_MYELINATING_OLIGOS_2014"))]
working_set <- unique(unname(unlist(MO_sets_mapped[is.element(names(MO_sets_mapped), set_id)])))
model_set_genes <- model_df[is.element(toupper(model_df$Gene), working_set),] %>%
  arrange(Genotype_Est) %>%
  filter(Genotype_Pval<.05)

# genes <- model_set_genes$Gene
genes <- c("S100b", "Mbp", "Cnp", "Enpp6", "Gal3st1", "Plp1", "Mag", "Abcb9")
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

## Prolif:

enrichment <- "proliferative"

set_id <- MO_legend$SetID[is.element(MO_legend$SetName, c("PHILLIPS_MOST_DIFF_EXP_IN_PROLIF_SUBTYPE_217_GENES", "MASICA_Upregulated_IDH1_mutation"))]
working_set <- unique(unname(unlist(MO_sets_mapped[is.element(names(MO_sets_mapped), set_id)])))
model_set_genes <- model_df[is.element(toupper(model_df$Gene), working_set),] %>%
  arrange(Genotype_Est) %>%
  filter(Genotype_Pval<.05)

#genes <- model_set_genes$Gene
genes <- c("Dll3", "Abhd3", "Rgcc", "Emp2", "Tmem100", "Ccnb1", "Pbk", "Dhfr", "Smc2")
plot_gene_expr_violin(projectname, expr, meta_cols, gene_name_col, genes, enrichment, predictor, score_type)

############################################# Genotype POS #############################################

predictor <- "Genotype"
score_type <- "pos"

enrich_idx <- intersect(grep(predictor, enrich_files, ignore.case=T), grep(score_type, enrich_files, ignore.case=T))
enrich <- fread(enrich_files[enrich_idx], data.table=F)

## Microglia:

enrichment <- "microglia"

set_id <- MO_legend$SetID[is.element(MO_legend$SetName, c("ZHANG_MOUSE_MICROGLIA_2014", "NAGY_2020_MICRO_MACRO", "AYHAN_2021_MICRO1" , "VELMESHEV_2019_MICROGLIA", "LAKE_2018_CER_MIC"))]
working_set <- unique(unname(unlist(MO_sets_mapped[is.element(names(MO_sets_mapped), set_id)])))
model_set_genes <- model_df[is.element(toupper(model_df$Gene), working_set),] %>%
  arrange(Genotype_Est) %>%
  filter(Genotype_Pval<.05)

#genes <- model_set_genes$Gene
genes <- c("Cd68", "C5ar1", "C3ar1", "Meg3", "Tet2", "Apobec1", "Grk3", "Csf2rb")
plot_gene_expr_violin(projectname, expr, meta_cols, gene_name_col, genes, enrichment, predictor, score_type)
