library(ggplot2)
library(RColorBrewer)
setwd("../Figures/")


ggsave(filename = "Fig1A.jpg",
       plot = ggplot(combined_hist_frame,aes(x=Method,y=pathway_num,fill=Method))+
         geom_boxplot()+
       xlab("Method")+
       ylab('Num of gene-sets')+
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        theme(text = element_text(size=10),legend.position = "none"),
       # width = 75,
       # height = 60,
       # units = "mm",
       width = 4,
       height = 4,
       units = "in",
       dpi = 600)






ggsave(filename = "Fig1B.jpg",
       plot =ggplot(scatter_plot, aes(x=daniel_pathway_num, y=tao_pathway_num)) + geom_point(shape=16,size=1,alpha=0.5) +
         #geom_density_2d()+
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         geom_abline(intercept = 0, slope = 1,colour='red',linetype="dashed")+
         xlab("Num of gene-sets by FET+FDR")+
         ylab('Num of gene-sets by gerr')+
         theme(text = element_text(size=10)),
       # width = 75,
       # height = 60,
       # units = "mm",
       width = 4,
       height = 4,
       units = "in",
       dpi = 600)



ggsave(filename = "Fig1C.jpg",
       plot = ggplot(overlap_go_terms_frame,aes(x=overlapCoef_go_pathway_num))+
         geom_histogram(
           aes(y=..density..),
           bins=30,
           colour="black", fill="lightblue"
         )+
         #geom_density(aes(y=..density..)) +
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         theme(text = element_text(size=10))+
         xlab("Overlap coefficients of gene-sets")+
         ylab('Density'),
       # width = 75,
       # height = 60,
       # units = "mm",
       width = 4,
       height = 4,
       units = "in",
       dpi = 600)


ggsave(filename = "Fig1D.jpg",
       plot = ggplot(overlap_go_terms_frame,aes(x=overlapCoef_go_gene_num))+
         geom_histogram(
           aes(y=..density..),
           bins=30,
           colour="black", fill="lightblue"
         )+
         #geom_density(aes(y=..density..)) +
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         theme(text = element_text(size=10))+
         xlab("Overlap coefficients of genes")+
         ylab('Density'),
       # width = 75,
       # height = 60,
       # units = "mm",
       width = 4,
       height = 4,
       units = "in",
       dpi = 600)





# figure 1E in the paper aturally 
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
       # width = 75,
       # height = 60,
       # units = "mm",
       width = 4,
       height = 4,
       units = "in",
       dpi = 600)



## Fig 1F in the paper  acturaly 

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
       # width = 75,
       # height = 60,
       # units = "mm",
       width = 4,
       height = 4,
       units = "in",
       dpi = 600)






