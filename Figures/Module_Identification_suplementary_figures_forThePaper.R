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
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        theme(text = element_text(size=10)),
       width = 75,
       height = 60,
       units = "mm",
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
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         geom_abline(intercept = 0, slope = 1,colour='red',linetype="dashed")+
         xlab("Num of gene-sets by FET+FDR")+
         ylab('Num of gene-sets by gerr')+
         theme(text = element_text(size=10)),
       width = 75,
       height = 60,
       units = "mm",
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
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         theme(text = element_text(size=10))+
         xlab("Overlap coefficients of gene-sets")+
         ylab('Density'),
       width = 75,
       height = 60,
       units = "mm",
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
         theme(text = element_text(size=10))+
         xlab("Jaccard index of gene-sets")+
         ylab('Density'),
       width = 85,
       height = 100,
       units = "mm",
       dpi = 600)


ggsave(filename = "Fig1D.jpg",
       plot = ggplot(overlap_go_terms_frame,aes(x=overlapCoef_go_gene_num))+
         geom_histogram(
           aes(y=..density..),
           bins=30,
           colour="black", fill="white"
         )+
         geom_density(aes(y=..density..)) +
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         theme(text = element_text(size=10))+
         xlab("Overlap coefficients of genes")+
         ylab('Density'),
       width = 75,
       height = 60,
       units = "mm",
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
       width = 85,
       height = 100,
       units = "mm",
       dpi = 600)



## fig1I
#Another interesting plot: ratio of GOI genes to non-GOI genes for:
#- gene sets identified by both methods
#- gene sets identified by your method
#- gene sets identified by FET+FDR
module_genesets=scan("dream_consensus_modules.gmt",what = "",sep = "\n")
module_genesets=strsplit(module_genesets,split = "\t")

module_names=sapply(module_genesets, function(x){x[1]})
network_names=lapply(module_names, function(x){
  strsplit(x,split = "_")[[1]][1]
})
unique_network_names=unique(network_names)

module_genesets=lapply(module_genesets,function(x){
  x_genes=x[-c(1,2)]
  return(x_genes)
})
names(module_genesets)=module_names
module_genesets=module_genesets[grep(names(module_genesets),pattern = "STRING")]


GOIgenes_frame=data.frame(matrix(0,nrow = length(common_module_names)*3,ncol = 3))
colnames(GOIgenes_frame)=c("module_names","method","ratio")
GOIgenes_frame$module_names=unlist(lapply(common_module_names,function(x)(rep(x,3))))

for(i in 1:length(common_module_names)){
  tao_selected_pathwayIDs=Tao_go_results[Tao_go_results$module==common_module_names[i],"pathway_id"]
  daniel_selected_pathwayIDs=Daniel_network_GO_results[Daniel_network_GO_results$module==common_module_names[i],"termId"]
  
  union_go_ids=union(tao_selected_pathwayIDs,daniel_selected_pathwayIDs)
  TL=length(union_go_ids)  #total length
  intsect_go_ids=intersect(tao_selected_pathwayIDs,daniel_selected_pathwayIDs)
  FET_ids=setdiff(daniel_selected_pathwayIDs,intsect_go_ids)
  gerr_ids=setdiff(tao_selected_pathwayIDs,intsect_go_ids)
  
  ## for GOI genes
  GOI_genes=module_genesets[[common_module_names[i]]]
  
  intsect_genes=unique(unlist(lapply(intsect_go_ids, function(id){go_list[[id]]})))
  FET_genes=unique(unlist(lapply(FET_ids, function(id){go_list[[id]]})))
  gerr_genes=unique(unlist(lapply(gerr_ids, function(id){go_list[[id]]})))
  
  intsect_ratio=length(intersect(GOI_genes,intsect_genes))/(length(setdiff(intsect_genes,GOI_genes))+1)
  FET_ratio=length(intersect(GOI_genes,FET_genes))/(length(setdiff(FET_genes,GOI_genes))+1)
  gerr_ratio=length(intersect(GOI_genes,gerr_genes))/(length(setdiff(gerr_genes,GOI_genes))+1)
  
  
  
  GOIgenes_frame[GOIgenes_frame$module_names==common_module_names[i],"method"]=c("FET+FDR only","gerr only","FET+FDR & gerr")
  GOIgenes_frame[GOIgenes_frame$module_names==common_module_names[i],"ratio"]=c(log2(FET_ratio),log2(gerr_ratio),log2(intsect_ratio))
  
}


GOIgenes_frame$method=factor(GOIgenes_frame$method,levels = c("FET+FDR only","FET+FDR & gerr","gerr only"))
ggsave(filename = "Fig1I.jpg",
       plot =ggplot(GOIgenes_frame, aes(x=ratio,y=..scaled..,fill=method)) +
         geom_density(alpha=0.4)+
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         theme(legend.title = element_text(size=8),legend.text=element_text(size=8), legend.key.size = unit(0.3, "cm"),
               text = element_text(size=10))+
         #theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
         theme(legend.position = "top",legend.title=element_blank())+
         theme(text = element_text(size=10)),
       width = 75,
       height = 60,
       units = "mm",
       dpi = 600)



## fig1J 
#I would suggest stacked barplot with three colors:
#  unique gene sets detected by one method, unique gene set by another and overlap.
#(perhaps normalized by total number of sets and sorted by number of unique sets detected by FET)
stackedBarplot_frame=data.frame(matrix(0,nrow = length(common_module_names)*3,ncol = 3))
colnames(stackedBarplot_frame)=c("module_names","method","ratio")
stackedBarplot_frame$module_names=unlist(lapply(common_module_names,function(x)(rep(x,3))))

for(i in 1:length(common_module_names)){
  tao_selected_pathwayIDs=Tao_go_results[Tao_go_results$module==common_module_names[i],"pathway_id"]
  daniel_selected_pathwayIDs=Daniel_network_GO_results[Daniel_network_GO_results$module==common_module_names[i],"termId"]
  
  union_go_ids=union(tao_selected_pathwayIDs,daniel_selected_pathwayIDs)
  TL=length(union_go_ids)  #total length
  intsect_go_ids=intersect(tao_selected_pathwayIDs,daniel_selected_pathwayIDs)
  FET_ids=setdiff(daniel_selected_pathwayIDs,intsect_go_ids)
  gerr_ids=setdiff(tao_selected_pathwayIDs,intsect_go_ids)
  
  ## for stackplot of gene sets
  stackedBarplot_frame[stackedBarplot_frame$module_names==common_module_names[i],"method"]=c("FET+FDR only","gerr only","FET+FDR & gerr")
  stackedBarplot_frame[stackedBarplot_frame$module_names==common_module_names[i],"ratio"]=c(length(FET_ids)/TL,length(gerr_ids)/TL,length(intsect_go_ids)/TL)
  
  
}

stackedBarplot_frame$method=factor(stackedBarplot_frame$method,levels = c("FET+FDR only","FET+FDR & gerr","gerr only"))
stackedBarplot_gerr_frame=stackedBarplot_frame[stackedBarplot_frame$method=="gerr only",]
stackedBarplot_frame$module_names=factor(stackedBarplot_frame$module_names,levels =stackedBarplot_gerr_frame$module_names[order(stackedBarplot_gerr_frame$ratio)] )
ggsave(filename = "Fig1J.pdf",
       plot =ggplot(stackedBarplot_frame, aes(x=module_names, y=ratio,fill=method)) +
         geom_bar(stat="identity", width=1)+
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         theme(legend.title = element_text(size=8),legend.text=element_text(size=8), legend.key.size = unit(0.3, "cm"),
               text = element_text(size=10))+
         theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
         theme(legend.position = "top",legend.title=element_blank())+
         xlab("Modules")+
         ylab('Proportion')+
         theme(text = element_text(size=10)),
       width = 75,
       height = 60,
       units = "mm",
       dpi = 600)

# run imagemagic
# convert -density 600 -median 3 Fig1J.pdf Fig1J.jpg

a=stackedBarplot_gerr_frame[order(stackedBarplot_gerr_frame$ratio),]
b=sapply(a$module_names, function(x){length(module_genesets[[x]])})

c=Daniel_network_GO_results[Daniel_network_GO_results$module=="PPI-STRING_Consensus_mod296",]
d=Tao_go_results[Tao_go_results$module=="PPI-STRING_Consensus_mod296",]

#tao_go_pathway_num_range=unique(sort(overlap_go_terms_frame$tao_go_pathway_num))


## Fig 1E F G 
#+taofang@ebi.ac.uk what are the ranks of the gene-sets identified by gerr in the FET+FDR procedure? 
#I guess they should on average rank higher than those not selected by gerr, not correct? Can be make a plot for that? 
#If we transform ranking to numeric values from 0 and 1, then we can use a density plot show that gerr selected gene-sets are generally ranking high.


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

gerr_FETFDR_list=c()
nonGerr_FETFDR_list=c()
for(i in 1:length(common_module_names)){
  tao_module_go_ids=Tao_go_results[Tao_go_results$module==common_module_names[i],"pathway_id"]
  
    
  daniel_module_go_ids=Daniel_network_GO_results[Daniel_network_GO_results$module==common_module_names[i],"termId"]
  daniel_module_go_ranks=rank(Daniel_network_GO_results[Daniel_network_GO_results$module==common_module_names[i],"P.noncentral.fdr"],
                              ties.method = "min")
  names(daniel_module_go_ranks)=daniel_module_go_ids
  daniel_module_go_norm_ranks=daniel_module_go_ranks/length(daniel_module_go_ranks)
  
  overlap_ids=intersect(tao_module_go_ids,daniel_module_go_ids)
  gerr_FETFDR=daniel_module_go_norm_ranks[overlap_ids]
  nonGerr_FETFDR=daniel_module_go_norm_ranks[setdiff(daniel_module_go_ids,overlap_ids)]
  
  gerr_FETFDR_list=c(gerr_FETFDR_list,gerr_FETFDR)
  nonGerr_FETFDR_list=c(nonGerr_FETFDR_list,nonGerr_FETFDR)
  
}


gerr_FETFDR_normRanks=data.frame(gerr_FETFFDR=gerr_FETFDR_list)
ggsave(filename = "Fig1E.jpg",
       plot = ggplot(gerr_FETFDR_normRanks,aes(x=gerr_FETFFDR))+
         geom_histogram(
           aes(y=..density..),
           bins=30,
           colour="black", fill="white"
         )+
         geom_density(aes(y=..density..)) +
         theme_bw()+
         theme(text = element_text(size=8))+
         xlab("Normalized ranks ")+
         ylab('Density'),
       width = 85,
       height = 100,
       units = "mm",
       dpi = 600)


##
nonGerr_FETFDR_normRanks=data.frame(nonGerr_FETFFDR=nonGerr_FETFDR_list)
ggsave(filename = "Fig1F.jpg",
       plot = ggplot(nonGerr_FETFDR_normRanks,aes(x=nonGerr_FETFFDR))+
         geom_histogram(
           aes(y=..density..),
           bins=30,
           colour="black", fill="white"
         )+
         geom_density(aes(y=..density..)) +
         theme_bw()+
         theme(text = element_text(size=8))+
         xlab("Normalized ranks ")+
         ylab('Density'),
       width = 85,
       height = 100,
       units = "mm",
       dpi = 600)




ggsave(filename = "Fig1G.jpg",
       plot = ggplot() + 
         geom_density(data =gerr_FETFDR_normRanks, aes(x = gerr_FETFFDR,fill = "FET+FDR & gerr "),color = "black",alpha=0.7) + 
         geom_density(data = nonGerr_FETFDR_normRanks, aes(x = nonGerr_FETFFDR,fill = "FET+FDR only"),color = "black", alpha = 0.7)+
         xlab('Normalized rank (gene-sets)')+ylab('density')+
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         theme(legend.title = element_text(size=8),legend.text=element_text(size=8), legend.key.size = unit(0.3, "cm"),
               text = element_text(size=10) )+
         theme(legend.position = "top",legend.title=element_blank()),
         #theme(legend.position = "top"),
         #scale_fill_manual(name = "", values = c("#E69F00", "#56B4E9"), labels = c("1" = "gerr", "2" = "nonGerr")) ,
         #scale_fill_manual(name = "", values = c("black", "red"), labels = c("1" = "gerr", "2" = "nonGerr")) ,
       width = 75,
       height = 60,
       units = "mm",
       dpi = 600)



## Fig 1H

gerr_FETFDR_genes_list=c()
nonGerr_FETFDR_genes_list=c()
for(i in 1:length(common_module_names)){
  print(i)
  tao_module_go_ids=Tao_go_results[Tao_go_results$module==common_module_names[i],"pathway_id"]
  tao_module_go_genes=unique(unlist(lapply(tao_module_go_ids, function(id){go_list[[id]]})))
  
  daniel_module_go_ids=Daniel_network_GO_results[Daniel_network_GO_results$module==common_module_names[i],"termId"]
  daniel_module_go_ranks=rank(Daniel_network_GO_results[Daniel_network_GO_results$module==common_module_names[i],"P.noncentral.fdr"],
                              ties.method = "min")
  names(daniel_module_go_ranks)=daniel_module_go_ids
  daniel_module_go_norm_ranks=daniel_module_go_ranks/length(daniel_module_go_ranks)
  

  daniel_module_go_genes=unique(unlist(lapply(daniel_module_go_ids, function(id){go_list[[id]]})))
  daniel_module_go_genes_norm_ranks=c()
  for(gene in daniel_module_go_genes){
    a=1
    for (id in daniel_module_go_ids) {
      if(gene %in% go_list[[id]] ){
        a=min(a,daniel_module_go_norm_ranks[id])
      }
    }
    daniel_module_go_genes_norm_ranks=c(daniel_module_go_genes_norm_ranks,a)
  }
  names(daniel_module_go_genes_norm_ranks)=daniel_module_go_genes
  
  
  overlap_genes=intersect(tao_module_go_genes,daniel_module_go_genes)
  gerr_FETFDR_genes=daniel_module_go_genes_norm_ranks[overlap_genes]
  nonGerr_FETFDR_genes=daniel_module_go_genes_norm_ranks[setdiff(daniel_module_go_genes,overlap_genes)]
  
  gerr_FETFDR_genes_list=c(gerr_FETFDR_genes_list,gerr_FETFDR_genes)
  nonGerr_FETFDR_genes_list=c(nonGerr_FETFDR_genes_list,nonGerr_FETFDR_genes)
  
}

saveRDS(gerr_FETFDR_genes_list,"gerr_FETFDR_genes_list.rds")
saveRDS(nonGerr_FETFDR_genes_list,"nonGerr_FETFDR_genes_list.rds")


gerr_FETFDR_genes_normRanks=data.frame(gerr_FETFFDR=gerr_FETFDR_genes_list)
nonGerr_FETFDR_genes_normRanks=data.frame(nonGerr_FETFFDR=nonGerr_FETFDR_genes_list)
ggsave(filename = "Fig1H.jpg",
       plot = ggplot() + 
         geom_density(data =gerr_FETFDR_genes_normRanks, aes(x = gerr_FETFFDR,fill = "FET+FDR & gerr"),color = "black",alpha=0.7) + 
         geom_density(data = nonGerr_FETFDR_genes_normRanks, aes(x = nonGerr_FETFFDR,fill = "FET+FDR only"),color = "black", alpha = 0.7)+
         xlab('Normalized rank (genes)')+ylab('density')+    
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
       theme(legend.title = element_text(size=8),legend.text=element_text(size=8), legend.key.size = unit(0.3, "cm"),
                                                                        text = element_text(size=10))+
       theme(legend.position = "top",legend.title=element_blank()),
       #scale_fill_manual(name = "", values = c("#E69F00", "#56B4E9"), labels = c("1" = "gerr", "2" = "nonGerr")) ,
       #scale_fill_manual(name = "", values = c("black", "red"), labels = c("1" = "gerr", "2" = "nonGerr")) ,
       width = 75,
       height = 60,
       units = "mm",
       dpi = 600)






