### Bar plots showing the number of cancer types that each TRP gene exhibited up-regulation or down-regulation. 
### The heat map at the bottom showing the fold-changes of TPRs in comparison between cancer and normal.
rm(list = ls())
##cancer with number of normal.samples>=5

cancer_names<-rev(c('HNSC',"GBM",'LUAD','LUSC','BRCA','THCA','COAD','ESCA','READ','STAD','CHOL','LIHC','BLCA','KICH','KIRC','KIRP','PRAD','UCEC'))
gruop_TRP<-read.table("TRPList.txt",header = T,as.is = T,sep = "\t",check.names=FALSE)
gruop_TRP<-gruop_TRP[order(gruop_TRP$`Approved symbol`),]
#return pValue
func_p<-function(v){return(wilcox.test(as.numeric(na.omit(v[1:n])),as.numeric(na.omit(v[(n+1):m])))$p.value)}
#return Fold Change
func_fc<-function(s){return(mean(na.omit(s[(n+1):m]))/mean(na.omit(s[1:n])))}
p_mat<-c()
fc_mat<-c()
logfc_mat<-c()
fdr_mat<-c()
#return fc, pValue, p.adjust(fdr), log2fc
for(i in 1:18){
  setwd("/TRP_FPKM")

  cancer_mat<-read.table(paste0("TCGA_",cancer_names[i],".txt"),row.names = 1,header = T,sep = "\t",
                         as.is=T,check.names=FALSE)
  cancer_mat<-cancer_mat[gruop_TRP$`Approved symbol`,]
  normal_sample<-read.table(paste0("/TRP_FPKM/",
                                   cancer_names[i],"_normal.txt"),header=F,sep="\t")
  normal_sample$V1 <- substr(normal_sample$V1,1,15)
  normal_sample<-unique(normal_sample$V1)
  tumor_sample<-read.table(paste0("/TRP_FPKM/",
                                  cancer_names[i],"_tumor.txt"),header=F,sep="\t")
  
  tumor_sample$V1 <- substr(tumor_sample$V1,1,15)
  tumor_sample<-unique(tumor_sample$V1)
  cancer_mat<-cancer_mat[,c(normal_sample,tumor_sample)]
  log_cancer_mat<-log2(cancer_mat+1)
  n<-length(normal_sample)
  m<-length(tumor_sample)
  pvalue<-apply(log_cancer_mat,1,func_p)
  fc<-apply(cancer_mat,1,func_fc)
  logfc<-log2(fc)
  fdr<-p.adjust(pvalue,method = "fdr",n = length(pvalue))
  p_mat<-cbind(p_mat,pvalue)
  fc_mat<-cbind(fc_mat,fc)
  logfc_mat<-cbind(logfc_mat,logfc)
  fdr_mat<-cbind(fdr_mat,fdr)
  cat(i,sep="\n")
}
colnames(p_mat)<-cancer_names
p_mat<-as.data.frame(p_mat)
colnames(fc_mat)<-cancer_names
fc_mat<-as.data.frame(fc_mat)
colnames(logfc_mat)<-cancer_names
logfc_mat<-as.data.frame(logfc_mat)
colnames(fdr_mat)<-cancer_names
fdr_mat<-as.data.frame(fdr_mat)
#
##########################
library(reshape2)
logfc_mat_pic<-melt(as.matrix(logfc_mat))
fdr_mat_pic<-melt(as.matrix(fdr_mat))
fc_mat_pic<-melt(as.matrix(fc_mat))
colnames(logfc_mat_pic)<-c("symbol","cancer","log2fc")
colnames(fdr_mat_pic)<-c("symbol","cancer","fdr")
colnames(fc_mat_pic)<-c("symbol","cancer","fc")
logfc_mat_pic<-as.data.frame(logfc_mat_pic)
fdr_mat_pic<-as.data.frame(fdr_mat_pic)
fc_mat_pic<-as.data.frame(fc_mat_pic)

pic_mat<-merge(logfc_mat_pic,merge(fdr_mat_pic,fc_mat_pic,by=c("symbol","cancer")),by=c("symbol","cancer"))
pic_mat[which(pic_mat$log2fc<(-2)),3]<--2
pic_mat[which(pic_mat$log2fc>2),3]<-2
pic_mat$log10fdr<-(-log10(pic_mat$fdr))
pic_mat<-as.data.frame(pic_mat)
pic_mat_sig<-pic_mat[which(pic_mat$fdr<0.05 & abs(pic_mat$log2fc)>1),]
colnames(pic_mat_sig)<-c("symbol","cancer", "log2fc","fdr","fc","log10fdr")

ggplot()+
  geom_tile(data=pic_mat,aes(symbol,cancer),fill="white",color=NA)+
  geom_point(data=pic_mat,aes(symbol,cancer,size=log10fdr,color=log2fc),
             shape=20)+scale_size_continuous(range = c(4,6))+
  scale_color_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,limits = c(-2, 2))+
  geom_point(data=pic_mat_sig,aes(symbol,cancer,size=log10fdr),alpha=2,
             shape=21)+scale_size_continuous(range = c(4,6))+
   theme(axis.text.x = element_text(angle = 90, hjust = 1))
###########################################
#number of cancer in which TRP gene up-regulated/down-regulated 
#up-regulated
up_gene<-pic_mat[which(pic_mat$fc>2 & pic_mat$fdr<0.05),]
up_gene_number<-as.data.frame(table(up_gene$symbol))
up_gene_number$lable<-"up"
#down-regulated
down_gene<-pic_mat[which(pic_mat$fc<0.5 & pic_mat$fdr<0.05),]
down_gene_number<-as.data.frame(table(down_gene$symbol))
down_gene_number$lable<-"down"
deg_num<-rbind(up_gene_number,down_gene_number)
colnames(deg_num)<-c("Symbol","Number","label")
### Plot

ggplot(deg_num)+
  geom_bar(aes(x=Symbol,y=Number,group=label,fill=label
               ),stat = 'identity',
           width = 1,position = position_dodge(width = 0),)+
  scale_alpha_discrete(range = c(1,0.5))+
  labs( y = 'Number of cancer') +
  scale_fill_manual(values = c("blue","red"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  
  scale_y_continuous(breaks = seq(-15, 15, 5),limits = c(-20,20))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
