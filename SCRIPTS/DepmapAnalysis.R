rm(list=ls())
work_path='path/'
setwd(work_path)
depmappath=paste0(work_path,'Depmap/')

library(openxlsx)
library(tidyverse)
library(dplyr)
library(Hmisc)
library(ggplot2)

####################### import data #########################

#Read TRP information
TRPdata=read.xlsx('group-TRP.xlsx')

#CCLE:Drug data
dat_Drug=read.csv(paste0(depmappath,'CCLE_Drug_data.csv'),
                  header=T)

#CCLE:Expression
dat_exp=read.csv(paste0(depmappath,'CCLE_RNAseq_genes_rpkm_20180929.gct'),
                 sep='\t',header=T,skip=2,stringsAsFactors = F,check.names = F)
dat_exp$Name=gsub('\\.[0-9]+','',dat_exp$Name)

# sample information
dat_sample=read.csv(paste0(depmappath,'sample_info.csv'),
                    header=T)
########################################################################################


####################### Correlation between CCLE-Drug IC50 and TRPs expression in CCLE #########################
#---------- Drug data
dat_Drug_2=dat_Drug %>% dplyr::select(CCLE.Cell.Line.Name,Compound,IC50..uM.) %>%
  spread(CCLE.Cell.Line.Name,IC50..uM.)
head(dat_Drug_2)[,1:3]
rownames(dat_Drug_2)=dat_Drug_2$Compound
dat_Drug_2=dat_Drug_2 %>% dplyr::select(-Compound)
dat_Drug_4=t(dat_Drug_2)
dat_Drug_4=as.data.frame(dat_Drug_4)

#---------- TRP expression in CCLE
dat_exp_TRP=dat_exp %>% filter(Name %in% TRPdata$ENSG)
head(dat_exp_TRP)[,1:3]
dat_exp_TRP=dat_exp_TRP %>% arrange(Description)
rownames(dat_exp_TRP)=dat_exp_TRP$Description
dat_exp_TRP=dat_exp_TRP %>% dplyr::select(-Name,-Description)
dat_exp_TRP=t(dat_exp_TRP)

#---------- Intersection Cell line
sample_intersect=intersect(rownames(dat_exp_TRP),rownames(dat_Drug_4))
dat_exp_TRP=dat_exp_TRP[sample_intersect,]
dat_Drug_4=dat_Drug_4[sample_intersect,]
dat_exp_TRP=apply(dat_exp_TRP, 2, as.numeric)
dat_Drug_4=apply(dat_Drug_4, 2, as.numeric)


#---------- Spearman correlation analysis
drug_all=rcorr(dat_exp_TRP,dat_Drug_4,type='spearman')
drug_cor=as.data.frame(drug_all$r)
drug_cor2=drug_cor[which(rownames(drug_cor) %in% colnames(dat_exp_TRP)),
                   which(colnames(drug_cor) %in% colnames(dat_Drug_4))]
drug_P=as.data.frame(drug_all$P)
drug_P2=drug_P[which(rownames(drug_P) %in% colnames(dat_exp_TRP)),
               which(colnames(drug_P) %in% colnames(dat_Drug_4))]
drug_cor2=reshape2::melt(as.matrix(drug_cor2))
colnames(drug_cor2)=c('TRP','DRUG','R')
drug_P2=reshape2::melt(as.matrix(drug_P2))
colnames(drug_P2)=c('TRP','DRUG','Pvalue')
drug_spearman=merge(drug_cor2,drug_P2,by=c('TRP','DRUG'))
drug_spearman$Padjust=p.adjust(drug_spearman$Pvalue,method='fdr')
TRPdata2=TRPdata %>% dplyr::select(Approved.symbol,Family)
drug_spearman=merge(drug_spearman,TRPdata2,by.x='TRP',by.y='Approved.symbol')
write.table(drug_spearman,'CCLE_DRUG_TRP.txt',sep='\t',quote=F,row.names=F)
write.csv(drug_spearman,'CCLE_DRUG_TRP.csv',quote=F,row.names=F)

#---------- Draw boxplot
drug_spearman=read.csv('CCLE_DRUG_TRP.txt',sep='\t',header=T)
drug_spearman_1=drug_spearman %>% filter(Padjust < 0.05) %>% distinct()
drug_spearman_2=drug_spearman_1 %>% filter(abs(R) > 0.15)#filter(R > 0.15)

colortrp=c('TRPA'='#8470FF',
           'TRPC'='#1E90FF',
           'TRPM'='#3CB371',
           'TRPML'='#FF8C00',
           'TRPP'='#EE82EE',
           'TRPV'='#FF6347')

p=ggplot(data=drug_spearman,aes(x=TRP,y=R,color=Family))+
  geom_boxplot(position = position_dodge(0.9))+
  geom_jitter(aes(color=Family),size=1)+
  ylab('TRP_drug')+xlab('')+
  scale_color_manual(values = colortrp)+
  scale_fill_manual(values = colortrp)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=0,hjust = 1,colour = 'black',size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(color = 'black',
                                  size=10,
                                  hjust=0.5))+
  geom_hline(yintercept=-0.15, linetype=44,color='red')+#"dotted"
  geom_hline(yintercept=0.15, linetype=44,color='red')+#"dotted"
  guides(colour = 'none')

pdf('CCLE_drug_cor_boxplot.pdf')
p+coord_flip()
dev.off()

########################################################################################

####################### Correlation between Drug-targeted genes and TRPs expression in CCLE #######################
#---------- TRP expression
dat_exp_TRP=dat_exp %>% filter(Name %in% TRPdata$ENSG)
head(dat_exp_TRP)[,1:3]
dat_exp_TRP=dat_exp_TRP %>% arrange(Description)
rownames(dat_exp_TRP)=dat_exp_TRP$Description
dat_exp_TRP=dat_exp_TRP %>% dplyr::select(-Name,-Description)
dat_exp_TRP=t(dat_exp_TRP)

#---------- Drug-targeted genes expression
DRUG_TARGET=dat_Drug %>% dplyr::select(Compound,Target) %>% distinct()
targetgene=unique(as.character(DRUG_TARGET$Target))
dat_target=dat_exp %>% filter(Description %in% targetgene)
rownames(dat_target)=dat_target$Description
dat_target=dat_target%>% dplyr::select(-Name,-Description)
dat_target=t(dat_target)

#---------- Spearman correlation analysis
library(Hmisc)
drug_all=rcorr(dat_exp_TRP,dat_target,type='spearman')
drug_cor=as.data.frame(drug_all$r)
drug_cor2=drug_cor[which(rownames(drug_cor) %in% colnames(dat_exp_TRP)),
                   which(colnames(drug_cor) %in% colnames(dat_target))]
drug_P=as.data.frame(drug_all$P)
drug_P2=drug_P[which(rownames(drug_P) %in% colnames(dat_exp_TRP)),
               which(colnames(drug_P) %in% colnames(dat_target))]
drug_cor2=reshape2::melt(as.matrix(drug_cor2))
colnames(drug_cor2)=c('TRP','DRUG','R')
drug_P2=reshape2::melt(as.matrix(drug_P2))
colnames(drug_P2)=c('TRP','DRUG','Pvalue')
drug_spearman=merge(drug_cor2,drug_P2,by=c('TRP','DRUG'))
drug_spearman$Padjust=p.adjust(drug_spearman$Pvalue,method='fdr')
TRPdata2=TRPdata %>% dplyr::select(Approved.symbol,Family)
drug_spearman=merge(drug_spearman,TRPdata2,by.x='TRP',by.y='Approved.symbol')
write.table(drug_spearman,'CCLE_TARGET_TRP.txt',sep='\t',quote=F,row.names=F)
write.csv(drug_spearman,'CCLE_TARGET_TRP.csv',quote=F,row.names=F)


#---------- Draw boxplot

drug_spearman=read.csv('CCLE_TARGET_TRP.txt',sep='\t',header=T)
drug_spearman_1=drug_spearman %>% filter(Padjust < 0.05) %>% distinct()
drug_spearman_2=drug_spearman_1 %>% filter(abs(R) > 0.15)

colortrp=c('TRPA'='#8470FF',
           'TRPC'='#1E90FF',
           'TRPM'='#3CB371',
           'TRPML'='#FF8C00',
           'TRPP'='#EE82EE',
           'TRPV'='#FF6347')

p=ggplot(data=drug_spearman,aes(x=TRP,y=R,color=Family))+
  geom_boxplot(position = position_dodge(0.9))+
  geom_jitter(aes(color=Family),size=1)+
  ylab('TRP_drug')+xlab('')+
  scale_color_manual(values = colortrp)+
  scale_fill_manual(values = colortrp)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=0,hjust = 1,colour = 'black',size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(color = 'black',
                                  size=10,
                                  hjust=0.5))+
  geom_hline(yintercept=-0.15, linetype=44,color='red')+
  geom_hline(yintercept=0.15, linetype=44,color='red')+
  guides(colour = 'none')

pdf('CCLE_drugtarget_cor_boxplot.pdf')
p+coord_flip()
dev.off()

############################################################################################



####################### Connection diagram:TRP+Drug+Drug-targeted genes #######################

#Read correlation results
drug_spearman_DRUG=read.csv('CCLE_DRUG_TRP.txt',sep='\t',header=T)
drug_spearman_TARGET=read.csv('CCLE_TARGET_TRP.txt',sep='\t',header=T)
DRUG_TARGET=dat_Drug %>% dplyr::select(Compound,Target) %>% distinct()
drug_spearman_DRUG=merge(drug_spearman_DRUG,DRUG_TARGET,by.x='DRUG',by.y='Compound')

#filter significant results
drug_spearman_DRUG_sig=drug_spearman_DRUG %>% filter(Padjust < 0.05 & abs(R) > 0.15)
drug_spearman_TARGET_sig=drug_spearman_TARGET %>% filter(Padjust < 0.05 & abs(R) > 0.15)
drug_spearman_all_sig=merge(drug_spearman_DRUG_sig,drug_spearman_TARGET_sig,by.x=c('TRP','Target'),by.y=c('TRP','DRUG'))

drug_spearman_all_sig_drug=drug_spearman_all_sig%>% dplyr::select(TRP,DRUG,R.x,Family.x)
drug_spearman_all_sig_drug$DRUG_Type='DRUG'
drug_spearman_all_sig_target=drug_spearman_all_sig%>% dplyr::select(TRP,Target,R.y,Family.x)
drug_spearman_all_sig_target$DRUG_Type='TARGET'
colnames(drug_spearman_all_sig_drug)=c('TRP','DRUG','R','Family','DRUG_Type')
colnames(drug_spearman_all_sig_target)=c('TRP','DRUG','R','Family','DRUG_Type')
drug_spearman_all_sig2=rbind(drug_spearman_all_sig_drug,drug_spearman_all_sig_target)
drug_spearman_all_sig2$Type='Negative'
drug_spearman_all_sig2[which(drug_spearman_all_sig2$R > 0),'Type']='Positive'
drug_spearman_all_sig2=drug_spearman_all_sig2 %>% dplyr::select(TRP,DRUG,DRUG_Type,Type,Family)
write.table(drug_spearman_all_sig2,'TRP_DRUG_TARGET_cys0.15.txt',sep='\t',
            row.names = F,quote = F)

#Type of each term[DRUG,TARGET,TRP]
NodeType=drug_spearman_all_sig2 %>% dplyr::select(DRUG,DRUG_Type) %>% distinct()
TRP=cbind(as.character(unique(drug_spearman_all_sig2$TRP)),'TRP')
colnames(TRP)=c('DRUG','DRUG_Type')
NodeType=rbind(NodeType,TRP)
colnames(NodeType)=c('Node','Type')
write.table(NodeType,'TRP_DRUG_TARGET_NodeType0.15.txt',sep='\t',row.names = F,quote = F)

############################################################################################
####################### Correlation between IC50 and necessity score #######################
#---------- Drug data
dat_Drug_1 = dat_Drug %>% dplyr::select(CCLE.Cell.Line.Name,Compound,IC50..uM.)
dat_sample_1=dat_sample %>% dplyr::select(DepMap_ID,CCLE_Name)
dat_Drug_2=merge(dat_Drug_1,dat_sample_1,by.x='CCLE.Cell.Line.Name',by.y='CCLE_Name')
library(tidyverse)
dat_Drug_3 = dat_Drug_2 %>% dplyr::select(DepMap_ID,Compound,IC50..uM.) %>% 
spread(DepMap_ID,IC50..uM.)
rownames(dat_Drug_3)=dat_Drug_3$Compound
dat_Drug_3=dat_Drug_3 %>% dplyr::select(-Compound)
dat_Drug_4=t(dat_Drug_3)
dat_Drug_4=as.data.frame(dat_Drug_4)

#---------- necessity score
dat_gene_effect=read.csv(paste0('Cellline_gene_necessity_score_matrix.txt'),sep='\t',
					   header=T)
head(dat_gene_effect)[,1:3]
#change sample names
sample0=colnames(dat_gene_effect)
sample0=gsub('\\.','-',sample0)
colnames(dat_gene_effect)=sample0
dat_gene_effect_TRP=dat_gene_effect[intersect(TRPdata$Approved.symbol,rownames(dat_gene_effect)),]
dat_gene_effect_TRP=t(dat_gene_effect_TRP)
dat_gene_effect_TRP=as.data.frame(dat_gene_effect_TRP)

#---------- intersection of cell lines
sample_intersect=intersect(rownames(dat_gene_effect_TRP),rownames(dat_Drug_4))
dat_gene_effect_TRP=dat_gene_effect_TRP[sample_intersect,]
dat_Drug_4=dat_Drug_4[sample_intersect,]
dat_gene_effect_TRP=apply(dat_gene_effect_TRP, 2, as.numeric)
dat_Drug_4=apply(dat_Drug_4, 2, as.numeric)

#---------- Spearman correlation analysis
drug_all=rcorr(dat_gene_effect_TRP,dat_Drug_4,type='spearman')
drug_cor=as.data.frame(drug_all$r)
drug_cor2=drug_cor[which(rownames(drug_cor) %in% colnames(dat_gene_effect_TRP)),
				 which(colnames(drug_cor) %in% colnames(dat_Drug_4))]
drug_P=as.data.frame(drug_all$P)
drug_P2=drug_P[which(rownames(drug_P) %in% colnames(dat_gene_effect_TRP)),
			 which(colnames(drug_P) %in% colnames(dat_Drug_4))]
drug_cor2=reshape2::melt(as.matrix(drug_cor2))
colnames(drug_cor2)=c('TRP','DRUG','R')
drug_P2=reshape2::melt(as.matrix(drug_P2))
colnames(drug_P2)=c('TRP','DRUG','Pvalue')
drug_spearman=merge(drug_cor2,drug_P2,by=c('TRP','DRUG'))
drug_spearman$Padjust=p.adjust(drug_spearman$Pvalue,method='fdr')
TRPdata2=TRPdata %>% dplyr::select(Approved.symbol,Family)
drug_spearman=merge(drug_spearman,TRPdata2,by.x='TRP',by.y='Approved.symbol')
write.table(drug_spearman,'Necessity_score_DRUG_TRP.txt',sep='\t',quote=F,row.names=F)
write.csv(drug_spearman,'Necessity_score_DRUG_TRP.csv',quote=F,row.names=F)

#---------- Draw boxplot
drug_spearman=read.csv('Necessity_score_DRUG_TRP.txt',sep='\t',header=T)
colortrp=c('TRPA'='#8470FF',
		 'TRPC'='#1E90FF',
		 'TRPM'='#3CB371',
		 'TRPML'='#FF8C00',
		 'TRPP'='#EE82EE',
		 'TRPV'='#FF6347')
p=ggplot(data=drug_spearman,aes(x=TRP,y=R,color=Family))+
geom_boxplot(position = position_dodge(0.9))+
geom_jitter(aes(color=Family),size=1)+
ylab('TRP_drug')+xlab('')+
scale_color_manual(values = colortrp)+
scale_fill_manual(values = colortrp)+
theme_bw()+
theme(axis.text.x = element_text(angle=0,hjust = 1,colour = 'black',size = 10),
	  axis.text.y = element_text(size = 10),
	  axis.title.y = element_text(size = 10),
	  legend.position = 'none',
	  panel.grid.major = element_blank(),
	  panel.grid.minor = element_blank(),
	  plot.title = element_text(color = 'black',
								size=10,
								hjust=0.5))+
geom_hline(yintercept=-0.15, linetype=44,color='red')+
geom_hline(yintercept=0.15, linetype=44,color='red')+
guides(colour = 'none')
p+coord_flip()
pdf('Necessity_score_DRUG_cor_boxplot.pdf')
p+coord_flip()
dev.off()

############################################################################################
####################### knockdown of several TRP genes #######################
#---------- necessity proliferation score
dat_gene_effect=read.csv(paste0('Cellline_gene_necessity_score_matrix.txt'),sep='\t',
					   header=T)

head(dat_gene_effect)[,1:3]
sample0=colnames(dat_gene_effect)
sample0=gsub('\\.','-',sample0)
colnames(dat_gene_effect)=sample0
head(dat_gene_effect)[,1:3]


#---------- Calculate the background proliferation score of each sample
Background=colMeans(dat_gene_effect)


#---------- Extract the TRP gene we need and add the background proliferation score
dat_gene_effect_TRP=dat_gene_effect[c(intersect(TRPdata$Approved.symbol,rownames(dat_gene_effect))),]
dat_gene_effect_TRP_background=rbind(dat_gene_effect_TRP,Background)
rownames(dat_gene_effect_TRP_background)[nrow(dat_gene_effect_TRP_background)]='Background'


#---------- Draw boxplot
dat_gene_effect_TRP_background2=reshape2::melt(as.matrix(dat_gene_effect_TRP_background))
colnames(dat_gene_effect_TRP_background2)=c('Term','Sample','Proliferation_score')

TRPdata2=TRPdata %>% dplyr::select(Approved.symbol,Family)
dat_gene_effect_TRP_background2=merge(dat_gene_effect_TRP_background2,TRPdata2,by.x='Term',by.y='Approved.symbol',all.x = T)
dat_gene_effect_TRP_background2[is.na(dat_gene_effect_TRP_background2$Family),'Family']='Background'

colortrp=c('TRPA'='#8470FF',
		 'TRPC'='#1E90FF',
		 'TRPM'='#3CB371',
		 'TRPML'='#FF8C00',
		 'TRPP'='#EE82EE',
		 'TRPV'='#FF6347',
		 'Background'='#999999'
)

p=ggplot(data=dat_gene_effect_TRP_background2,aes(x=Term,y=Proliferation_score,color=Family))+#Sample_Type
stat_boxplot(geom="errorbar",width=0.15,aes(x=Term,color=Family),position = position_dodge(0.9))+
geom_boxplot(position = position_dodge(0.9),outlier.shape = NA)+
geom_jitter(aes(color=Family),size=0.05)+
ylab('Necessity score')+xlab('')+
scale_color_manual(values = colortrp)+
scale_fill_manual(values = colortrp)+
theme_bw()+
theme(axis.text.x = element_text(angle=90,hjust = 1,colour = 'black',size = 10),
	  axis.text.y = element_text(size = 10),
	  axis.title.y = element_text(size = 10),
	  legend.position = 'none',
	  panel.grid.major = element_blank(),
	  panel.grid.minor = element_blank(),
	  plot.title = element_text(color = 'black',
								size=10,
								hjust=0.5))+
geom_hline(yintercept=median(Background), linetype=44,color='black')+
ylim(-60,60)+
guides(colour = 'none')

pdf('Necessity_score_background_boxplot.pdf',width=10,height = 5)
print(p)
dev.off()


