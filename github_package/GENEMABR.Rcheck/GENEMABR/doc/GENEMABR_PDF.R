## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval= FALSE----------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE)){
#      install.packages("BiocManager")}
#  BiocManager::install("GENEMABR")

## --------------------------------------------------------------------------


## ---- eval=TRUE, message=FALSE---------------------------------------------
library(GENEMABR)

## load other required package to run GENEMABR
if(!require(glmnet)){
    install.packages("glmnet")
    library(glmnet)
}
if(!require(Matrix)){
    install.packages("Matrix")
    library(Matrix)
}
if(!require(igraph)){
    install.packages("igraph")
    library(igraph)
}


## --------------------------------------------------------------------------
hypergeometric_test_results=read.csv(file = "../data_raw/Daniel_S5_mod283.txt",header = T,sep = "\t")
## filterd results with only P.noncentral value less then 0.05
hypergeometric_test_results_filtered=hypergeometric_test_results[hypergeometric_test_results$P.noncentral<0.05,]

head(hypergeometric_test_results_filtered)
dim(hypergeometric_test_results_filtered)

## --------------------------------------------------------------------------
#Gene module from the paper
gene_list=c("TRPC4AP","CDC37","TNIP1","IKBKB","NKIRAS2","NFKBIA","TIMM50","RELB","TNFAIP3","NFKBIB","HSPA1A","NFKBIE","SPAG9","NFKB2","ERLIN1","REL","TNIP2","TUBB6","MAP3K8")
#help("regression_selected_pathways")

#Here use regression_selected_pathways with default gene pathway matrix and set the alpha value as 0.5
enrichment_results=regression_selected_pathways(gene_input=gene_list,gene_pathway_matrix="default",alpha=0.5)
enrichment_results

## --------------------------------------------------------------------------
## if you use the default pathway databases(GO Ontologyand REACTOME).
## After you extracted enriched pathways, you can use find_root_ids function to find thier GO sub-root or REACTOME roots(ID) to help you better understanding the biological meanings of pathways.
## Here we use GO sub-root instead of GO root nodes as there are only three roots in the GO ontology and it's not so specific 
GO_Reactome_root_id=find_root_ids(selected_pathways=names(enrichment_results$selected_pathways_coef))
GO_Reactome_root_id

# Or if you want to obatain root notes names instead of ID, you can use function from_id2name to get names from ids
GO_Reactome_root_id_names=from_id2name(GO_Reactome_root_id)
GO_Reactome_root_id_names

# Or you  can use function get_steps function to calculate the distance from selected pathways to GO or Reactome roots
step2root=get_steps(selected_pathways=names(enrichment_results$selected_pathways_coef))
step2root
# To view specic position of GO/REACOTEM pathways in ontology trees.
# You can use Visualization tool at https://www.ebi.ac.uk/QuickGO/ and https://reactome.org/PathwayBrowser/ 

## --------------------------------------------------------------------------
sessionInfo()

