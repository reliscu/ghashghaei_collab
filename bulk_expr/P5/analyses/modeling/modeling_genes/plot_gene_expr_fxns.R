library(ggplot2)
library(viridis)

plot_gene_expr_violin <- function(projectname, model_df, predictor, score_type, percentile_cut=0, gene_name, expr, meta_cols, gene_name_col, region_vec, genotype_vec){
  
  pred_col <- colnames(model_df)[grep(paste0(predictor, "_Est"), colnames(model_df))]
  
  if(score_type=="neg"){
    
    model_df <- model_df %>% arrange_at(pred_col)
    
  } else {
    
    model_df <- model_df %>% arrange_at(pred_col, desc)
    
  }

  model_df <- model_df[model_df$Mean_Expr_Percentile>=percentile_cut,]
  
  for(i in 1:5){
    
    gene_name <- model_df$Gene[i]
    expr_gene <- as.numeric(
      expr[is.element(toupper(expr[,gene_name_col]), toupper(gene_name)), -meta_cols]
    )
    
    df_gene <- data.frame(
      Gene=gene_name,
      Expr=expr_gene,
      Genotype=genotype_vec,
      Region=region_vec,
      Interaction=paste(
        genotype_vec, region_vec
      )
    )
    
    if(i==1){
      df <- df_gene
    } else {
      df <- rbind(df, df_gene)
    }
    
  } ## for(i in 1:3){
  
  p <- 
    ggplot(df, aes_string(
      x="Interaction", y="Expr", fill=predictor, group="Interaction"
    )) +
    geom_violin() + 
    geom_boxplot(width=.1, fill="white", show.legend=F) + 
    geom_jitter(color="grey", show.legend=F) +
    theme_bw() + 
    theme(
      plot.title=element_text(hjust=.5),
      legend.position="bottom",
      legend.title=element_text(hjust=.5, face="bold"),
      axis.title.x=element_blank(),
      axis.text.x=element_text(face="bold", size=10, angle=45, vjust=.5)
    ) +
    labs(
      title=paste(toupper(predictor), toupper(score_type))
    ) +
    ylab("TPM") +
    scale_fill_viridis(discrete=T) 
  
  pdf(file=paste0("figures/", projectname, "_", toupper(predictor), "_", toupper(score_type), "_expr_violin.pdf"), width=10, height=7)
  
  print(
    p + facet_grid(.~Gene) + theme(
      strip.text.x=element_text(size=12, face="bold"),
      strip.background.x=element_rect(fill="white")
    )
  )
  
  dev.off()
  
  
  
} ## plot_gene_expr_violin <- function(

plot_gene_expr_barplot <- function(projectname, model_df, predictor, score_type, percentile_cut=0, gene_name, expr, meta_cols, gene_name_col, region_vec, genotype_vec){
  
  pred_col <- colnames(model_df)[grep(paste0(predictor, "_Est"), colnames(model_df))]
  
  if(score_type=="neg"){
    model_df <- model_df %>% arrange_at(pred_col)
  } else {
    model_df <- model_df %>% arrange_at(pred_col, desc)
  }
  
  model_df <- model_df[model_df$Mean_Expr_Percentile>=percentile_cut,]
  
  for(i in 1:5){
    
    gene_name <- model_df$Gene[i]
    expr_gene <- as.numeric(
      expr[is.element(toupper(expr[,gene_name_col]), toupper(gene_name)), -meta_cols]
    )
    df_gene <- data.frame(
      Gene=gene_name,
      Expr=expr_gene,
      Genotype=genotype_vec,
      Region=region_vec,
      Interaction=paste(
        genotype_vec, region_vec
      ),
      Label=paste(
        genotype_vec, region_vec, seq(length(genotype_vec))
      )
    )
    if(i==1){
      df <- df_gene
    } else {
      df <- rbind(df, df_gene)
    }
    
  } ## for(i in 1:3){
  
  p <- 
    ggplot(df, aes_string(x="Label", y="Expr", fill="Interaction")) +
    geom_bar(stat="identity") + 
    theme_classic() + 
    theme(
      plot.title=element_text(hjust=.5),
      legend.position="bottom",
      legend.title=element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_line(size=0)
    ) +
    labs(
      title=paste(toupper(predictor), toupper(score_type))
    ) +
    ylab("TPM") +
    scale_fill_viridis(name="Condition", discrete=T) 
  
  pdf(file=paste0("figures/", projectname, "_", toupper(predictor), "_", toupper(score_type), "_expr_barplot.pdf"), width=10, height=7)
  
  print(
    p + facet_grid(.~Gene) + theme(
      strip.text.x=element_text(size=12, face="bold"),
      strip.background.x=element_rect(fill="white")
    )
  )
  
  dev.off()
  
} ## plot_gene_expr_barplot <- function(