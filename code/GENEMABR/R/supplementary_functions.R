#' fisher_exact_test
#' one side? two side?
#' https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
#'
#' This function allows you to compute fish exact pvalue of gene list for selected  pathways
#'          genes_in_common      genes_in_reference
#'pathway          a                   c
#'no_pathway       b                   d
#' @param gene_input A vecor of genes to be annotated. It should have same ID types(Ensembl ID, HUGO gene symbol) as the genes in \emph{gene_pathway_matrix}.
#' @param gene_pathway_matrix A binary background matrix whose columns are the pathways/gene sets and 
#'whose rows are all the genes from pathways/gene sets . It could be in sparse matrix format ((inherit from class "sparseMatrix" as in package Matrix) to save memory.
#'For gene i and pathway j, the value of matrix(i,j) is 1 is gene i belonging to pathway j otherwise 0.
#' @param lambda We use glmnet function to do regression. \emph{lambda} is an argument in \strong{glmnet}. See \strong{glmnet} function for more details
#' Here we use default value 0.007956622 after preliminary study. It can be overridden by giving \emph{nlambda} and \emph{lambda.min.ratio arguments}.
#' @param alpha The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as
#'(1-α)/2||β||_2^2+α||β||_1.
#'alpha=1 is the lasso penalty, and alpha=0 the ridge penalty. Default value: 0.5.
#' @param ... Other paramaters for glmnet function.
#' @return  A list of three element: 
#' \itemize{
#'   \item selected_pathways_coef - Regression coefficients value for selected pathways
#'   \item selected_pathways_fisher_pvalue - Fisher exact pvalue for selected pathways
#'   \item selected_pathways_num_genes - The number of genes for selected pathways in background
#' }
#' @keywords 
#' @export
#' @examples
#' a=regression_selected_pathways(gene_input =c("TRPC4AP","CDC37","TNIP1","IKBKB","NKIRAS2","NFKBIA","TIMM50","RELB","TNFAIP3","NFKBIB","HSPA1A","NFKBIE","SPAG9","NFKB2","ERLIN1","REL","TNIP2","TUBB6","MAP3K8"),gene_pathway_matrix=gene_pathway_matrix,lambda=0.007956622,alpha=0.5)


fisher_exact_test=function(selected_pathways,gene_input,gene_pathway_matrix){ 
  
  all_genes=rownames(gene_pathway_matrix))
  module1_common_genes=intersect(all_genes,gene_input) 
  
  selected_pathways_fisher_pvalue=vector()
  selected_pathways_num_genes=vector()
  for(index in 1:length(selected_pathways)){
    fisher_pathway=selected_pathways[index]
    #          genes_in_common      genes_in_reference
    #pathway          a                   c
    #no_pathway       b                   d
    a=sum(gene_pathway_matrix[module1_common_genes,fisher_pathway])
    b=length(module1_common_genes)-a
    c=sum(gene_pathway_matrix[,fisher_pathway])-a
    d=length(all_genes)-a-c-b
    contigency_table=matrix(c(a,b,c,d),nrow = 2)
    fisher_result=fisher.test(contigency_table)
    selected_pathways_fisher_pvalue[fisher_pathway]=fisher_result[['p.value']]
    selected_pathways_num_genes[fisher_pathway]=sum(gene_pathway_matrix[,fisher_pathway])
  }
  return(list(selected_pathways_fisher_pvalue=selected_pathways_fisher_pvalue,selected_pathways_num_genes=selected_pathways_num_genes))
}
