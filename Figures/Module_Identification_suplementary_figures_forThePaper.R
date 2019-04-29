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


#####https://docs.google.com/document/d/1s6KyQH0M7WK8sGqXpxeD9UKjB_c3tYGoPF9vIpzcNtM/edit


#relationship between regression coefficients and Fisher’s exact p-values. Ideas:
#Distribution of (a) #significant terms (FDR corrected P-value from noncentral hypergeometric) and (b) 
#selected terms (regression). E.g., two boxplots or histograms. 
#Scatterplot: #significant terms (Fisher) vs. #selected terms for each module
#how that the selected terms are less redundant than the ones from Fisher: E.g., for each significant term, compute max Jaccard index across all other significant terms of that module. Show distribution of these max Jaccard indexes for all significant terms and modules. Do the same for the selected terms.


#Distribution of (a) #significant terms (FDR corrected P-value from noncentral hypergeometric) and (b) 
#selected terms (regression). E.g., two boxplots or histograms. 

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


combined_hist_frame=data.frame(name=c(rep('Regression-based method',nrow(hist_tao)),rep('FDR-corrected noncentral hypergeometric test',nrow(hist_daniel))),
                               pathway_num=c(hist_tao$pathway_num,hist_daniel$pathway_num))

ggsave(filename = "Hist_plot0.5.jpg",
       plot = ggplot(combined_hist_frame,aes(x=pathway_num))+
         geom_histogram(
           aes(y=..density..),
           bins=30,
           colour="black", fill="white"
         )+
         geom_density(aes(y=..density..)) +
         facet_wrap(~name,ncol=1)+
         theme_bw()+
         xlab("Selected pathways for each module")+
         ylab('Density'),
         #ggtitle("Distribution of num for selected terms "),
       dpi = 600)




##########################Scatterplot: #significant terms (Fisher) vs. #selected terms for each module
common_module_names=intersect(tao_module_names,daniel_module_names)
scatter_plot=data.frame(module=common_module_names,row.names = common_module_names,
                        tao_pathway_num=hist_tao[common_module_names,]$pathway_num,
                        daniel_pathway_num=hist_daniel[common_module_names,]$pathway_num
                        )

ggsave(filename = "scatter_plot0.5.jpg",
       plot =ggplot(scatter_plot, aes(x=daniel_pathway_num, y=tao_pathway_num)) + geom_point(shape=1) +
         theme_bw()+
         geom_abline(intercept = 0, slope = 1,colour='red',linetype="dashed")+
         xlab("FDR-corrected noncentral hypergeometric test")+
         ylab('Regression-based method'),
         #ggtitle("Scatter plot of num for selected terms "),
       dpi = 600)

## decrcit range 0-20 and 0-200 in the paper



###################################### 
#overlap GO selecetd terms between two methods
Daniel_network_GO_results=Daniel_network_results[grepl("GO",Daniel_network_results$termId),]
daniel_GO_module_names=unique(Daniel_network_GO_results$module)
Tao_go_results=Tao_results[grep("GO",Tao_results$pathway_id),]
tao_go_module_names=unique(Tao_go_results$module)
common_module_names=intersect(tao_go_module_names,daniel_GO_module_names)

tao_go_id=unique(Tao_go_results$pathway_id)  #1949
daniel_go_id=unique(Daniel_network_GO_results$termId) #10838
total_overlap_go_id=intersect(tao_go_id,daniel_go_id)  #1872

overlap_go_terms_frame=data.frame(matrix(0,nrow = length(common_module_names),ncol = 3),row.names = common_module_names)
colnames(overlap_go_terms_frame)=c('tao_go_pathway_num','daniel_go_pathway_num','overlap_go_pathway_num')

for(i in 1:length(common_module_names)){
  tao_selected_pathwayIDs=Tao_go_results[Tao_go_results$module==common_module_names[i],"pathway_id"]
  overlap_go_terms_frame[common_module_names[i],"tao_go_pathway_num"]=length(tao_selected_pathwayIDs)
  daniel_selected_pathwayIDs=Daniel_network_GO_results[Daniel_network_GO_results$module==common_module_names[i],"termId"]
  overlap_go_terms_frame[common_module_names[i],"daniel_go_pathway_num"]=length(daniel_selected_pathwayIDs)
  overlap_go_terms_frame[common_module_names[i],"overlap_go_pathway_num"]=length(intersect(tao_selected_pathwayIDs,daniel_selected_pathwayIDs))
}




ggsave(filename = "overlap_go_terms_scatter_plot0.5.jpg",
       plot =ggplot(overlap_go_terms_frame, aes(x=tao_go_pathway_num, y=overlap_go_pathway_num)) + geom_point(shape=1) +
         theme_bw()+
         geom_abline(intercept = 0, slope = 1,colour='red',linetype="dashed")+
         xlab("Num of pathways from regression-based method")+
         ylab('Num of overlapped pathways'),
         #ggtitle("Overlap GO terms between two methods "),
       dpi = 600)



