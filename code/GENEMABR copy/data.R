#' Homo sapiens GO ontology tree 
#'@description A rds R object contains GO ontology relationships (tree structure)
#'@format Directed igraph format
#'@details  The igraph format tree was constructed by using data from http://geneontology.org/docs/download-ontology/ (May 2108)
#'It has three root notes representing Molecular Function,Cellular Component and Biological Process (http://geneontology.org/docs/ontology-documentation/)
#'It was saved as raw data in this package
#'The path to this data could be list by system.file("extdata", "go_ontology.rds", package = "GENEMABR")
#'And you can read the data by readRDS(system.file("extdata", "go_ontology.rds", package = "GENEMABR"))
#'"human_go_ontology"



#' Homo sapiens REACTOME ontology tree 
#'@description A rds R object contains Reactome ontology relationships (tree structure)
#'@format Directed igraph format
#'It has several root notes representing REACTOME pathway categories (https://reactome.org/PathwayBrowser/)
#'It was saved as raw data in this package
#'@details  The igraph format tree was constructed by using  data from https://reactome.org/download-data (May 2108)
#'It was saved as raw data in this package
#'The path to this data could be list by system.file("extdata", "human_reactome_ontology.rds", package = "GENEMABR")
#'And you can read the data by readRDS(system.file("extdata", "human_reactome_ontology.rds", package = "GENEMABR")
#'"human_reactome_ontology"



#' Homo sapiens GO ontology and REACTOME ontology gene-pathway realtionship
#'@description A rds R object contains GO ontology and REACTOME ontology gene-pathway realtionship
#'@format Formal class 'dgCMatrix' [package "Matrix"] 
#'A binary matrix whose columns are the pathways/gene sets from GO ontology and REATOME database and 
#'whose rows are all the genes(represented by gene HUGO gene symbols) from  GO ontology and REATOME database.
#'For gene i and pathway j, the value of matrix(i,j) is 1 is gene i belonging to pathway j otherwise 0
#'"gene_pathway_matrix"

