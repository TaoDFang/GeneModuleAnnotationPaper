library(ggplot2)
library(RColorBrewer)
setwd("../Figures/")

Daniel_results_name="dream_consensus_modules.functional_enrichment.txt"
Daniel_results=read.csv(Daniel_results_name,header = TRUE,sep = "\t",stringsAsFactors = FALSE)
Daniel_network_results=Daniel_results[Daniel_results$network=="PPI-STRING",]
GO_Rec_index=c(grep("GO",Daniel_network_results$termId),grep("REACTOME",Daniel_network_results$term))
Daniel_network_results=Daniel_network_results[GO_Rec_index,]
Daniel_network_results=Daniel_network_results[Daniel_network_results$P.noncentral.fdr<0.05,]

Tao_results= read.csv(file = "dream_consensus_modules_results_alpha0.5.csv",header = TRUE,stringsAsFactors = FALSE) # dream_consensus_modules_results_old.csv dream_consensus_modules_results_alpha0.5.csv
Tao_reactome_results=Tao_results[grep("HSA",Tao_results$pathway_id),]


## Fig 1A. boxlot
tao_module_names=unique(Tao_results$module)
hist_tao=data.frame(module=tao_module_names,pathway_num=rep(0,length(tao_module_names)),row.names = tao_module_names)  
for(i in 1:length(tao_module_names)){
  hist_tao[tao_module_names[i],"pathway_num"]=sum(Tao_results$module==tao_module_names[i]) 
}


daniel_module_names=unique(Daniel_network_results$module)
hist_daniel=data.frame(module=daniel_module_names,pathway_num=rep(0,length(daniel_module_names)),row.names=daniel_module_names)
for (i in 1:length(daniel_module_names)) {
  hist_daniel[daniel_module_names[i],"pathway_num"]=sum(Daniel_network_results$module==daniel_module_names[i])
}



combined_hist_frame=data.frame(name=c(rep('gerr',nrow(hist_tao)),rep('FET+FDR',nrow(hist_daniel))),
                               pathway_num=c(hist_tao$pathway_num,hist_daniel$pathway_num))

ggsave(filename = "Fig1A.jpg",
       plot = ggplot(combined_hist_frame,aes(x=name,y=pathway_num))+
         geom_boxplot()+
       xlab("Method")+
       ylab('Num of gene-sets')+
        theme(text = element_text(size=8)),
       width = 3,
       height = 4,
       dpi = 600)

###median=XXX, interquartile range/IQR=XXX
summary(combined_hist_frame[combined_hist_frame$name=="gerr","pathway_num"]) #9.000
IQR(combined_hist_frame[combined_hist_frame$name=="gerr","pathway_num"])     # 8.000


summary(combined_hist_frame[combined_hist_frame$name=="FET+FDR","pathway_num"]) #54
IQR(combined_hist_frame[combined_hist_frame$name=="FET+FDR","pathway_num"])     # 53



## Fig 1B scatter plot
common_module_names=intersect(tao_module_names,daniel_module_names)
scatter_plot=data.frame(module=common_module_names,row.names = common_module_names,
                        tao_pathway_num=hist_tao[common_module_names,]$pathway_num,
                        daniel_pathway_num=hist_daniel[common_module_names,]$pathway_num
                        )

ggsave(filename = "Fig1B.jpg",
       plot =ggplot(scatter_plot, aes(x=daniel_pathway_num, y=tao_pathway_num)) + geom_point(shape=1,size=1) +
         theme_bw()+
         geom_abline(intercept = 0, slope = 1,colour='red',linetype="dashed")+
         xlab("Num of gene-sets by FET+FDR")+
         ylab('Num of gene-sets by gerr')+
         theme(text = element_text(size=8)),
       width = 3,
       height = 4,
       dpi = 600)

### correlation
cor(scatter_plot$tao_pathway_num,scatter_plot$daniel_pathway_num,method = "pearson") ##0.70



 
## fig 1C  adn 1D
overlap_coefficient<-function(set1,set2){
  return(abs(length(intersect(set1,set2)))/min(length(set1),length(set2)))
}

jaccard_index<-function(set1,set2){
  return(abs(length(intersect(set1,set2)))/abs(length(union(set1,set2))))
}


### load gene information for gene sets
# library(data.table)
# go_gaf_filename="/Users/taofang/Documents/DreamChallengeModuleAnnotation/HPC_Disease_module_identification_DREAM_change/goa_human.gaf"
# go_gaf=fread(go_gaf_filename,header = FALSE,skip = 23,fill=TRUE, na.strings="NA",stringsAsFactors = FALSE)
# go_gaf_simple=go_gaf[,c(3,5)]
# go_gaf_simple=go_gaf_simple[which(go_gaf_simple[,1] !=""),]
# go_gaf_simple=as.data.frame(go_gaf_simple)
# all_go_ids=unique(go_gaf_simple[,2])
# go_list=vector("list",length = length(all_go_ids))
# names(go_list)=all_go_ids
# for (i in 1:length(all_go_ids)){
#   go_list[[i]]=go_gaf_simple[go_gaf_simple[,2]==all_go_ids[i],1]
# }
#saveRDS(go_list, "go_full_list.rds")

go_list=readRDS("go_full_list.rds")


Daniel_network_GO_results=Daniel_network_results[grepl("GO",Daniel_network_results$termId),]
daniel_GO_module_names=unique(Daniel_network_GO_results$module)
Tao_go_results=Tao_results[grep("GO",Tao_results$pathway_id),]
tao_go_module_names=unique(Tao_go_results$module)
common_module_names=intersect(tao_go_module_names,daniel_GO_module_names)

tao_go_id=unique(Tao_go_results$pathway_id)  #1949
daniel_go_id=unique(Daniel_network_GO_results$termId) #10838
total_overlap_go_id=intersect(tao_go_id,daniel_go_id)  #1872

overlap_go_terms_frame=data.frame(matrix(0,nrow = length(common_module_names),ncol = 7),row.names = common_module_names)
colnames(overlap_go_terms_frame)=c('tao_go_pathway_num','daniel_go_pathway_num','overlap_go_pathway_num',"overlapCoef_go_pathway_num","overlapCoef_go_gene_num","jaccardJndex_go_pathway_num","jaccardJndex_go_gene_num")

for(i in 1:length(common_module_names)){
  tao_selected_pathwayIDs=Tao_go_results[Tao_go_results$module==common_module_names[i],"pathway_id"]
  overlap_go_terms_frame[common_module_names[i],"tao_go_pathway_num"]=length(tao_selected_pathwayIDs)
  daniel_selected_pathwayIDs=Daniel_network_GO_results[Daniel_network_GO_results$module==common_module_names[i],"termId"]
  overlap_go_terms_frame[common_module_names[i],"daniel_go_pathway_num"]=length(daniel_selected_pathwayIDs)
  overlap_go_terms_frame[common_module_names[i],"overlap_go_pathway_num"]=length(intersect(tao_selected_pathwayIDs,daniel_selected_pathwayIDs))
  
  ##overlapp coefficient of gene sets/pathways
  overlap_go_terms_frame[common_module_names[i],"overlapCoef_go_pathway_num"]=overlap_coefficient(tao_selected_pathwayIDs,daniel_selected_pathwayIDs)
  
  ## jaccard index  of gene sets/pathways
  overlap_go_terms_frame[common_module_names[i],"jaccardJndex_go_pathway_num"]=jaccard_index(tao_selected_pathwayIDs,daniel_selected_pathwayIDs)
  
  
  ##overlapp coefficient of genes
  
  tao_selected_genes=unique(unlist(lapply(tao_selected_pathwayIDs, function(id){go_list[[id]]})))
  daniel_selected_genes=unique(unlist(lapply(daniel_selected_pathwayIDs, function(id){go_list[[id]]})))
  
  overlap_go_terms_frame[common_module_names[i],"overlapCoef_go_gene_num"]=overlap_coefficient(tao_selected_genes,daniel_selected_genes)
  
  ##jaccard index of genes
  overlap_go_terms_frame[common_module_names[i],"jaccardJndex_go_gene_num"]=jaccard_index(tao_selected_genes,daniel_selected_genes)
  
}


ggsave(filename = "Fig1C.jpg",
       plot = ggplot(overlap_go_terms_frame,aes(x=overlapCoef_go_pathway_num))+
         geom_histogram(
           aes(y=..density..),
           bins=30,
           colour="black", fill="white"
         )+
         geom_density(aes(y=..density..)) +
         theme_bw()+
         theme(text = element_text(size=8))+
         xlab("Overlap coefficients of gene-sets")+
         ylab('Density'),
       width = 3,
       height = 4,
       dpi = 600)

ggsave(filename = "Fig1C_jaccardIndex.jpg",
       plot = ggplot(overlap_go_terms_frame,aes(x=jaccardJndex_go_pathway_num))+
         geom_histogram(
           aes(y=..density..),
           bins=30,
           colour="black", fill="white"
         )+
         geom_density(aes(y=..density..)) +
         theme_bw()+
         theme(text = element_text(size=8))+
         xlab("Jaccard index of gene-sets")+
         ylab('Density'),
       width = 3,
       height = 4,
       dpi = 600)


ggsave(filename = "Fig1D.jpg",
       plot = ggplot(overlap_go_terms_frame,aes(x=overlapCoef_go_gene_num))+
         geom_histogram(
           aes(y=..density..),
           bins=30,
           colour="black", fill="white"
         )+
         geom_density(aes(y=..density..)) +
         theme_bw()+
         theme(text = element_text(size=8))+
         xlab("Overlap coefficients of genes")+
         ylab('Density'),
       width = 3,
       height = 4,
       dpi = 600)


ggsave(filename = "Fig1D_jaccardIndex.jpg",
       plot = ggplot(overlap_go_terms_frame,aes(x=jaccardJndex_go_gene_num))+
         geom_histogram(
           aes(y=..density..),
           bins=30,
           colour="black", fill="white"
         )+
         geom_density(aes(y=..density..)) +
         theme_bw()+
         theme(text = element_text(size=8))+
         xlab("Jaccard index of genes")+
         ylab('Density'),
       width = 3,
       height = 4,
       dpi = 600)

#tao_go_pathway_num_range=unique(sort(overlap_go_terms_frame$tao_go_pathway_num))


# Fig 1D





