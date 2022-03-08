### The spearman correlation coefficients between TRP gene expressions and drug IC50.
###Preprocessing pathway data
pathway<-read.table("drug_pathway.txt",
                    header = T,sep = "\t",as.is = T,quote = "")
drug_id<-pathway[,c('Identifier','Name')]
pathway<-pathway[,c('Identifier','Name','Targeted.process.pathway','Putative.Target')]
target_list <- strsplit(with(pathway,Putative.Target),", ")
count_n <- sapply(target_list,length)
target_pathway<-data.frame(Identifier=rep(with(target_pathway,Identifier),count_n),
                           Name=rep(with(target_pathway,Name),count_n),
                           #symbol=rep(with(target_pathway,symbol),count_n),
                           Targeted.process.pathway=rep(with(target_pathway,Targeted.process.pathway),count_n),
                           #label=rep(with(target_pathway,label),count_n),
                           Putative.Target=unlist(target_list))
target_pathway$Putative.Target<-gsub("\\(.*\\)","",target_pathway$Putative.Target)
write.table(target_pathway,"pathway_infomation.txt",sep = "\t",quote = F,col.names = T,row.names = F)


###drug
rm(list = ls())
###The expression profile of GDSC cell line 
cell_line_exp<-read.table("GDSCExp.txt",
                          header = T,sep = "\t",as.is = T,quote = "")
####TRP
gruop_TRP<-read.table("TRPList.txt",header = T,as.is = T,sep = "\t",check.names=FALSE)
##TRP symbol
gene<-c(gruop_TRP$`Approved symbol`)
##expression profile of TRP gene
cell_line_exp<-cell_line_exp[which(cell_line_exp$GENE_SYMBOLS %in% gene),]
rownames(cell_line_exp)<-cell_line_exp$GENE_SYMBOLS
cell_line_exp<-cell_line_exp[,-c(1,2)]
cell_line_symbol<-colnames(cell_line_exp)
cell_line_symbol<-substring(cell_line_symbol, 6)
colnames(cell_line_exp)<-cell_line_symbol
#############IC50,tableS4A
drug_data<-read.table("GDSCDrugIC50.txt",sep = "\t",as.is = T,
                      header = T,quote = "",check.names = F)

rownames(drug_data)<-drug_data$`Cell line cosmic identifiers`
drug_data<-drug_data[,-c(1,2)]
#####Cell lines with NA value greater than 30% were removed
func_na<-function(v){return(length(which(is.na(v)))/length(v))}
sample_na<-apply(drug_data,1,func_na)
index_30<-which(sample_na>0.3)
drug_data<-drug_data[-index_30,]
inter_sample<-intersect(colnames(cell_line_exp),rownames(drug_data))
cell_line_exp_inter<-cell_line_exp[,inter_sample]
cell_line_exp_inter<-as.matrix(cell_line_exp_inter)
drug_inter<-drug_data[inter_sample,]
drug_inter<-t(as.matrix(drug_inter))

####knn for make up the missing value
library(impute)
impute_drug<-impute.knn(drug_inter,k = 5, rowmax = 0.5,colmax = 0.8,maxp = 1500,rng.seed=362436069)
impute_drug<-as.matrix(impute_drug[["data"]])
###spearman
library(Hmisc)
cor_mat<-rcorr(cbind(t(cell_line_exp_inter),t(impute_drug)),type = "spearman")
###return p&r
###r
r_mat<-cor_mat[["r"]][-c(1:26),]
r_mat<-r_mat[,c(1:26)]
library(reshape2)
r_mat<-as.data.frame(melt(r_mat))
####p
p_mat<-cor_mat[["P"]][-c(1:26),]
p_mat<-p_mat[,c(1:26)]
#library(reshape2)
p_mat<-as.data.frame(melt(p_mat))
colnames(r_mat)<-c("drug","symbol","r")
colnames(p_mat)<-c("drug","symbol","p")
####fdr
p_mat$fdr<-p.adjust(p_mat$p,method = "fdr",n=length(p_mat$p))
####p<0.05,|r|>0.15
r_mat_0.15<-r_mat[which(abs(r_mat$r)>0.15),]
p_mat_sig<-p_mat[which(p_mat$fdr<0.05),]
pr_mat<-merge(r_mat,p_mat,by=c("drug",'symbol'))
colnames(pr_mat)<-c("Identifier","symbol","r",'pvalue','p.adjust')
pr_mat<-merge(pr_mat,drug_id,by="Identifier")
pr_mat<-pr_mat[,c("Identifier",'Name',"symbol","r",'pvalue','p.adjust')]
target_pathway<-read.table("pathway_infomation.txt",
                           header = T,as.is = T,sep = "\t")
drug_symbol<-merge(pr_mat,drug_id,by="Identifier")
drug_symbol<-drug_symbol[,c('Identifier','Name',"symbol","r",'pvalue','p.adjust')]
#all
#write.table(pr_mat,"TRP_Drug_spearman.txt",sep = "\t",quote = F,col.names = T,row.names = F)
###sig
pr_mat_sig<-merge(r_mat_0.15,p_mat_sig,by=c("drug",'symbol'))
pr_mat_sig<-pr_mat_sig[order(pr_mat_sig$symbol),]
pr_mat_sig$label<-NA
pr_mat_sig[which(pr_mat_sig$r<0),6]<-"negative"
pr_mat_sig[which(pr_mat_sig$r>0),6]<-"positive"
colnames(pr_mat_sig)<-c("Identifier","symbol","r",'pvalue','p.adjust','drug_label')
pr_mat_sig<-merge(pr_mat_sig,drug_id,by="Identifier")
pr_mat_sig<-pr_mat_sig[,c("Identifier",'Name',"symbol","r",'pvalue','p.adjust','drug_label')]
#write.table(pr_mat_sig,"gene_drug_significant.txt",sep = "\t",quote = F,col.names = T,row.names = F)
length(unique(pr_mat_sig$symbol))
length(unique(pr_mat_sig$Name))

r_mat$symbol<-factor(r_mat$symbol,levels=gene[order(gene)])
color_symbol=c(c(rep('#FF8C00',3)),c(rep('#EE82EE',3)),c(rep('#8470FF',1)),c(rep('#1E90FF',6)),c(rep('#3CB371',8)),
          c(rep('#FF6347',5)))
### Figures8A
ggplot(r_mat,aes(y=r,x=symbol,fill=symbol))+

  geom_jitter(aes(color=symbol),width =0.2,shape = 21,size=0.3)+ 
  scale_fill_manual(values = color_symbol)+  
  scale_color_manual(values=color_symbol)+ 
  theme_bw()+ 
  theme(legend.position="none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=c(0.15,-0.15), linetype="dotted",size=1,colour="black")+
  theme(axis.text.x=element_text(angle=90,hjust=0.5, vjust=0.5))+
  xlab("TRP symbol")+ylab("coefficient")
############################################################################################################
#for data of figure8B, figure8B were visualized with Cytoscape
rm(list = ls())
###spearman coefficient between TRP gene and IC50 from GDSC
###The expression profile of GDSC cell line 
pathway<-read.table("drug_pathway.txt",
                    header = T,sep = "\t",as.is = T,quote = "")
drug_id<-pathway[,c('Identifier','Name')]
cell_line_exp_all<-read.table("GDSCExp.txt",
                          header = T,sep = "\t",as.is = T,quote = "")
####TRP infomation
gruop_TRP<-read.table("TRPList.txt",header = T,as.is = T,sep = "\t",check.names=FALSE)
TRP_gene<-c(gruop_TRP$`Approved symbol`)
##TRP expression
cell_line_exp<-cell_line_exp_all[which(cell_line_exp_all$GENE_SYMBOLS %in% TRP_gene),]
rownames(cell_line_exp)<-cell_line_exp$GENE_SYMBOLS
cell_line_exp<-cell_line_exp[,-c(1,2)]
cell_line_symbol<-colnames(cell_line_exp)
cell_line_symbol<-substring(cell_line_symbol, 6)
colnames(cell_line_exp)<-cell_line_symbol
######
clinical_gene<-read.table("clinical_gene.txt",
                          header = T,sep = "\t",as.is = T,quote = "",check.names = F)

clinical_gene_exp<-cell_line_exp_all[which(cell_line_exp_all$GENE_SYMBOLS %in% clinical_gene$Gene),]

rownames(clinical_gene_exp)<-clinical_gene_exp$GENE_SYMBOLS
clinical_gene_exp<-clinical_gene_exp[,-c(1,2)]
clinical_symbol<-colnames(clinical_gene_exp)
clinical_symbol<-substring(clinical_symbol, 6)
colnames(clinical_gene_exp)<-clinical_symbol

cor_mat<-rcorr(cbind(t(clinical_gene_exp),t(cell_line_exp)),type = "spearman")
r_mat<-cor_mat[["r"]][-c(1:117),]
r_mat<-r_mat[,c(1:117)]
library(reshape2)
r_mat<-as.data.frame(melt(r_mat))
####p
p_mat<-cor_mat[["P"]][-c(1:117),]
p_mat<-p_mat[,c(1:117)]
#library(reshape2)
p_mat<-as.data.frame(melt(p_mat))
colnames(r_mat)<-c("TRP","clinical","r")
colnames(p_mat)<-c("TRP","clinical","p")
####Adjust pValue
p_mat$fdr<-p.adjust(p_mat$p,method = "fdr",n=length(p_mat$p))

color_symbol=c(c(rep('#FF8C00',3)),c(rep('#EE82EE',3)),c(rep('#8470FF',1)),c(rep('#1E90FF',6)),c(rep('#3CB371',8)),
               c(rep('#FF6347',5)))
gene<-c(gruop_TRP$`Approved symbol`)
r_mat$TRP<-factor(r_mat$TRP,levels=gene[order(gene)])
ggplot(r_mat,aes(y=r,x=TRP,fill=TRP))+ 
  geom_jitter(aes(color=TRP),width =0.2,shape = 21,size=0.5)+ 
  scale_fill_manual(values = color_symbol)+  
  scale_color_manual(values=color_symbol)+ 
  theme_bw()+ 
  theme(legend.position="none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=c(0.2,-0.2), linetype="dotted",size=1,colour="black")+
  theme(axis.text.x=element_text(angle=90,hjust=0.5, vjust=0.5))+
  xlab("TRP symbol")+ylab("coefficient")
####significant : p<0.05,|r|>0.2
r_mat_0.25<-r_mat[which(abs(r_mat$r)>0.2),]
p_mat_sig<-p_mat[which(p_mat$fdr<0.05),]
pr_mat<-merge(r_mat,p_mat,by=c("TRP","clinical"))
colnames(pr_mat)<-c("TRP.symbol","Target","spearman",'pvalue','p.adjust')
#write.table(pr_mat,"TRP_Target_spearman.txt",sep = "\t",quote = F,col.names = T,row.names = F)
pr_mat_sig<-merge(r_mat_0.25,p_mat_sig,by=c("TRP","clinical"))
colnames(pr_mat_sig)<-c("symbol","Putative.Target","r",'pvalue','p.adjust')
length(unique(pr_mat_sig$Putative.Target))
#####add labels
pr_mat_sig$target.label<-NA
pr_mat_sig[which(pr_mat_sig$r<0),6]<-"negative"
pr_mat_sig[which(pr_mat_sig$r>0),6]<-"positive"
##############
target_pathway<-read.table("pathway_infomation.txt",
                           header = T,as.is = T,sep = "\t")
pr_mat_sig_pair<-merge(target_pathway,pr_mat_sig,by="Putative.Target")
pr_mat_sig_pair<-pr_mat_sig_pair[-which(pr_mat_sig_pair$Targeted.process.pathway=="other"),]
drug_symbol<-read.table("target_pathway_clinical_target.txt",
                        sep="\t",as.is=T,header=T)

intersect(pr_mat_sig_pair$symbol,drug_symbol$symbol)
unique(drug_symbol$symbol)
unique(pr_mat_sig_pair$symbol)
drug_symbol_target<-merge(pr_mat_sig_pair,drug_symbol,by=c("symbol","Name","Targeted.process.pathway","Putative.Target"))

















