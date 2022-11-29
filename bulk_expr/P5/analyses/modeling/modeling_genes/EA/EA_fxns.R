library(dplyr)
library(fgsea) 
library(GSEABase)
library(data.table)

EA_fxn <- function(
  projectname,
  model_df,
  gene_sets,
  set_legend,
  set_source,
  predictor,
  min_size=5, max_size=500,
  score_type=c("pos", "neg")
  ){
  
  est_col <- intersect(
    grep(predictor, colnames(model_df)), grep("Est", colnames(model_df))
  )
  
  model_df <- model_df %>% arrange(desc(!!as.name(colnames(model_df)[est_col])))
  
  ranks <- model_df[,est_col]
  names(ranks) <- toupper(model_df$Gene)
  
  gene_sets <- lapply(gene_sets, toupper)
  
  enrich <- fgsea(
    pathways=gene_sets, 
    stats=ranks, 
    minSize=min_size, 
    maxSize=max_size, 
    scoreType=score_type
    )
  
  set_legend <- set_legend %>% dplyr::select(SetID, SetName)
  enrich <- merge(set_legend, enrich, by.x="SetID", by.y="pathway")
  enrich <- enrich %>% dplyr::arrange(padj)
  
  file_name <- paste(
    projectname, 
    toupper(set_source),
    "SETS",
    toupper(predictor), 
    toupper(score_type), 
    nrow(model_df), 
    "genes.csv", sep="_"
    )
  
  fwrite(enrich, file=paste0("fgsea_", file_name))
  
}