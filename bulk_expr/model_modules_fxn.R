library(dplyr)
library(data.table)

model_modules <- function(projectname, fm_dir, sampleinfo, genotype_vec, region_vec){
  
  networks <- list.files(path=fm_dir, pattern="signum", full.names=T)
  idx <- unlist(lapply(networks, function(x) length(list.files(path=x))))>0
  networks <- networks[idx]
  
  df_list <- lapply(1:length(networks), function(i){
    
    print(networks[i])
    
    mod_eig <- fread(list.files(path=networks[i], pattern="eigengenes", full.names=T), data.table=F)
    mod_eig <- mod_eig[match(make.names(sampleinfo$FM_Label), make.names(mod_eig$Sample)),]
    
    if(!identical(make.names(sampleinfo$FM_Label), make.names(mod_eig$Sample))) {
      stop("!identical(sampleinfo$FM_Label, make.names(mod_eig$Sample))")
    }
    
    mod_stats <- fread(list.files(path=networks[i], pattern="statistics", full.names=T)[1], data.table=F)
    
    interaction_model <- lapply(2:ncol(mod_eig), function(x){
      lm(as.numeric(mod_eig[,x]) ~ genotype_vec + region_vec + genotype_vec:region_vec)
    })
    
    basic_model <- lapply(2:ncol(mod_eig), function(x){
      lm(as.numeric(mod_eig[,x]) ~ genotype_vec + region_vec)})
    
    network <- sapply(strsplit(networks[i], "/", fixed=T), function(x) x[length(x)])
    genotype_est <- sapply(interaction_model, function(x) summary(x)$coeff[2,1])
    genotype_pval <- sapply(interaction_model, function(x) summary(x)$coeff[2,4])
    region_est <- sapply(interaction_model, function(x) summary(x)$coeff[3,1])
    region_pval <- sapply(interaction_model, function(x) summary(x)$coeff[3,4])
    interaction_est <- sapply(interaction_model, function(x) summary(x)$coeff[4,1])
    interaction_pval <- sapply(interaction_model, function(x) summary(x)$coeff[4,4])
    basic_mdl_r2 <- sapply(basic_model, function(x) summary(x)$adj.r.squared)
    interaction_mdl_r2 <- sapply(interaction_model, function(x) summary(x)$adj.r.squared)
    
    return(data.frame(
      Network=network,
      Module=colnames(mod_eig)[2:ncol(mod_eig)],
      Genotype_Est=genotype_est,
      Genotype_Pval=genotype_pval,
      Region_Est=region_est,
      Region_Pval=region_pval,
      Interaction_Est=interaction_est,
      Interaction_Pval=interaction_pval,
      Basic_Model_R2=basic_mdl_r2,
      Interaction_Model_R2=interaction_mdl_r2,
      Mean_Expr=mod_stats$MeanExpr
    ))
    
  }) ## for(i in 1:length(networks)){
  
  df <- do.call(rbind, df_list)
  
  fwrite(df, file=paste(projectname, "modeling", nrow(model_df), "modules.csv", sep="_"))
  
} 

