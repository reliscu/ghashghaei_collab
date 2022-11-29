library(tibble)
library(data.table)

module_enrichments <- function(
  projectname,
  model_df,
  fm_dir
){
  
  ## Subset to modules with at least one significant predictor:
  
  signif_mods <- apply(
    model_df[,grep("Pval", colnames(model_df))], 1, function(x){sum(x<.05)}
  )>0
  
  # signif_mods <- apply(
  #   model_df[,grep("Pval", colnames(model_df))], 1, function(x){sum(x<1)}
  # )>0
  
  if(sum(signif_mods)>0){
    
    model_df <- model_df[signif_mods,]
    
    ## Get top predictor per module:
    
    min_idx <- model_df %>% 
      dplyr::select(grep("Pval", colnames(.))) %>% apply(., 1, which.min)
    
    pred <- sapply(
      strsplit(unlist(lapply(1:length(min_idx), function(x){
        colnames(model_df)[grep("Est", colnames(model_df))][min_idx[x]]
      })), "_"), "[", 1)
    
    est <- unlist(lapply(1:length(min_idx), function(x){
      model_df[x,grep("Est", colnames(model_df))][,min_idx[x]]
    }))
    
    pval <- unlist(lapply(1:length(min_idx), function(x){
      model_df[x,grep("Pval", colnames(model_df))][,min_idx[x]]
    }))

    model_df <- model_df %>%
      dplyr::select(Network, Module) %>%
      mutate(
        Predictor=pred,
        Estimate=est,
        Pval=pval
      )
    
    ## Get signif. enrichments associated with top modules:
    
    networks <- unique(model_df$Network)
    
    model_df$MO_Sets <- NA
    model_df$Broad_Sets <- NA
      
    for(i in 1:length(networks)){
      
      ea_files <- list.files(
        path=file.path(fm_dir, networks[i]), 
        pattern="GSHyperG",
        full.names=T
      )
      
      for(j in 1:length(ea_files)){
        
        enrich <- fread(ea_files[j], data.table=F)
    
        modules <- model_df$Module[model_df$Network%in%networks[i]]
        bc_pval <- .05/nrow(enrich)
        
        sets_per_mod <- unlist(lapply(modules, function(mod){
          signif_sets <- enrich[,colnames(enrich)%in%mod]<bc_pval
          if(sum(signif_sets)>0){

            pvals <- signif(enrich[signif_sets, colnames(enrich)%in%mod], 2)
            set_list <- c(paste0(enrich$SetName[signif_sets], " (", pvals,")"))
            set_list <- paste(set_list[order(pvals)], collapse=" | ")
            
          } else {
            return(NA)
          }
        })) ## unlist(lapply(modules, function(mod){
        
        if(grepl("MY_SETS", ea_files[j])){
          model_df$MO_Sets[model_df$Network%in%networks[i]] <- sets_per_mod
        } else {
          model_df$Broad_Sets[model_df$Network%in%networks[i]] <- sets_per_mod
        }
        
      } ## for(j in 1:length(ea_files)){
      
    } ## for(i in 1:length(networks)){
    
    fwrite(model_df, file=paste(
      projectname, "enrichments", nrow(model_df), "signif_modules.csv", sep="_"
    ))
    
  } ## if(sum(signif_mods)>0){

} ## mod_enrichments <- function(

