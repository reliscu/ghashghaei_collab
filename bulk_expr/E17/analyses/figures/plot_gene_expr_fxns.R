library(ggplot2)
library(viridis)

source("/home/rebecca/code/misc/upper_first.R")

plot_gene_expr_violin <- function(projectname, expr, meta_cols, gene_name_col, genes, enrichment, predictor, score_type){

  expr[,gene_name_col] <- toupper(expr[,gene_name_col])
  genes <- 
    genes[is.element(toupper(genes), expr[,gene_name_col])]
  set_idx <- match(toupper(genes), expr[,gene_name_col])
  expr_gene <- expr[set_idx,] %>% 
    dplyr::select(-ENSEMBL)
  
  pdf(file=paste0("figures/", projectname, "_", toupper(predictor), "_", toupper(score_type), "_", enrichment, "_expr_violin.pdf"), width=8, height=7)
  
  for(i in 1:length(genes)){
    
    gene <- genes[i]
    expr_gene <- as.numeric(
      expr[is.element(expr[,gene_name_col], toupper(gene)), -meta_cols]
    )
    
    df <- data.frame(
      Gene=gene,
      Expr=expr_gene,
      Genotype=genotype_vec,
      Region=region_vec,
      Interaction=paste(
        genotype_vec, region_vec
      )
    )
    
    p <- 
      ggplot(df, aes_string(x="Interaction", y="Expr", fill=predictor, group="Interaction")) +
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
    
    p
    
    print(
      p + facet_grid(.~Gene) + theme(
        strip.text.x=element_text(size=12, face="bold"),
        strip.background.x=element_rect(fill="white")
      )
    )
  
    
  } ## for(i in 1:length(genes)){

  dev.off()
  
} ## plot_gene_expr_violin <- function(

