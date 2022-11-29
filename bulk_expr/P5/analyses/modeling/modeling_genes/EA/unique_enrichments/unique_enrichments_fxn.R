library(dplyr)
library(data.table)

unique_enrichments <- function(
  projectname, predictor, all_pred, set_source, score_type, pval_cut){
 
  idx_all_pred <- intersect(
    grep(toupper(set_source), list.files(path="..")),
    grep(toupper(score_type), list.files(path=".."))
  )
  
  rest_pred <- all_pred[!all_pred%in%predictor]
  
  idx_pred <- intersect(idx_all_pred, grep(toupper(predictor), list.files(path="..")))
  idx_all_pred1 <- intersect(idx_all_pred, grep(toupper(rest_pred[1]), list.files(path="..")))
  idx_all_pred2 <- intersect(idx_all_pred, grep(toupper(rest_pred[2]), list.files(path="..")))
  
  enrich_pred <- read.csv(list.files(path="..", full.names=T)[idx_pred])
  
  if(sum(enrich_pred$padj<pval_cut)>0){
    
    enrich_all_pred1 <- read.csv(list.files(path="..", full.names=T)[idx_all_pred1])
    enrich_all_pred2 <- read.csv(list.files(path="..", full.names=T)[idx_all_pred2])
    
    enrich_pred$Predictor <- predictor
    enrich_all_pred1$Predictor <- rest_pred[1]
    enrich_all_pred2$Predictor <- rest_pred[2]
    
    enrich <- rbind(enrich_pred, enrich_all_pred1, enrich_all_pred2) %>%
      dplyr::select(-leadingEdge) %>% filter(padj<pval_cut)
    
    enrich_min <- enrich %>% 
      group_by(SetID) %>% 
      slice_min(padj) %>% 
      arrange(padj) %>%
      filter(Predictor==predictor) 
    
    if(nrow(enrich_min)>0){
      
      file_name <- gsub(
        projectname, 
        paste0(projectname, "_unique_enrichments"), 
        list.files(path="..")[idx_pred]
      )
      
      fwrite(enrich_min, file=file_name)
      
    } ## if(nrow(enrich_pred)>0){
    
  } ## if(sum(enrich_pred$padj<pval_cut)>0){
  
} ## unique_enrichments <- function(
