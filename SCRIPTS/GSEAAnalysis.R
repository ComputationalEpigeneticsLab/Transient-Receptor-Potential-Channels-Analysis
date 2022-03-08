#GSEA
library(Hmisc)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(openxlsx)
library(ggplot2)
library(tidyverse)

#Read TRP information
Trimdata=read.xlsx(paste(workdir1,"TRIM.xlsx",sep = ""))

################## Calculate the correlation between TRP family and other genes #################
#load Expression data from TCGA
expdata_cancer_all=read.csv('expression/TCGA_Allcancer.txt',header=T,sep='\t',stringsAsFactors=F)
#delete genes with low expression
lowindex=which(apply(expdata_cancer_all,1,function(v){return((sum(v==0)/length(v))>=0.5)}))
if(length(lowindex)>0){
expdata_cancer_all2=expdata_cancer_all[-lowindex,]
}else{
expdata_cancer_all2=expdata_cancer_all
}
#transpose matrix
expdata_cancer_all3=t(expdata_cancer_all2)
#calculate the correlation
myresult=rcorr(expdata_cancer_all3, type='spearman')
#extract correlation results
myresult_r=myresult$r[which(rownames(myresult$r) %in% TRPdata$ENSG),]
myresult_P=myresult$P[which(rownames(myresult$P) %in% TRPdata$ENSG),]
myresult_r2=reshape2::melt(as.matrix(myresult_r))
myresult_P2=reshape2::melt(as.matrix(myresult_P))
colnames(myresult_r2)[3]='R'
colnames(myresult_P2)[3]='Pvalue'
myresult_all=merge(myresult_r2,myresult_P2,by=c('Var1','Var2'))
colnames(myresult_all)[1:2]=c('TRP','Gene')
#Save the result of correlation analysis
#write.table(myresult_all,paste0('COR/Allcancer_pearson.txt'),sep='\t',quote=F,row.names=F)
########################################################################################
################## Annotation gene name ##################
#Read annotation file
mart_export=read.csv('mart_export.txt',header=T)
#Remove version
mart_export$ENSG=gsub('\\.[0-9]+','',mart_export$Gene.stable.ID.version)
mart_export2=mart_export %>% select(ENSG,Gene.name,Gene.type,NCBI.gene..formerly.Entrezgene..ID) %>% distinct()
#Convert Ensembl ID to gene symbol
dat=merge(myresult_all,mart_export2,by.x='TRP',by.y='ENSG')
colnames(dat)[5]='TRP_symbol'
dat=dat %>% select(TRP_symbol,Gene,R,Pvalue)
mart_export3=mart_export2 %>% filter(Gene.type=='protein_coding')
dat=merge(dat,mart_export3,by.x='Gene',by.y='ENSG')
########################################################################################
################## Gene set enrichment analysis(GSEA) ##################
#Read hallmark gene set file
kegmt<-read.gmt("h.all.v7.4.symbols.gmt")

#TRP family
genetrp=unique(dat$TRP_symbol) %>% sort()
#separately calculate Gene set enrichment analysis for TRP members
for(j in 1:length(genetrp)){
print(genetrp[j])
#Select correlation results for TRP members
dat_1=dat %>% filter(TRP_symbol==genetrp[j] )
dat_1=dat_1 %>% select(Gene.name,R)

#Remove duplication
dupid = unique(dat_1$Gene.name[which(duplicated(dat_1$Gene.name))])
if(length(dupid)>0){
dat_2= dat_1[-which(dat_1$Gene.name %in% dupid),]
}else{
dat_2=dat_1
}

#Sort by correlation coefficient
dat_2$R=as.numeric(dat_2$R)
dat_2=dat_2[order(dat_2$R,decreasing = T),]
dat_2=dat_2[which(dat_2$Gene.name!=''),]
kgene=dat_2
geneList<-kgene[,2] 
names(geneList)=kgene[,1]
geneList=sort(geneList,decreasing = T) 

#GSEA
KEGG<-tryCatch({
clusterProfiler::GSEA(geneList,TERM2GENE = kegmt,pvalueCutoff = 1,eps=0,
				  minGSSize=5,maxGSSize=10000)
},error = function(e){
#Collect error reporting information
cat(paste0(genetrp[j],'_KEGGnot1\n'),file='notattend_all.txt',append=T)
}
)

if(!is.null(KEGG)){
if(nrow(KEGG@result)>0){
#Save GSEA results of corresponding TRP member and export TRP information successfully enriched
saveRDS(KEGG,paste0('GSEA/',genetrp[j],'_GSEA_all.Rds'))
cat(paste0(genetrp[j],'\n'),file='attend_all.txt',sep='\t',append=T)
}else{ 
#Collect information that are not enriched in the gene set
cat(paste0(genetrp[j],'_KEGGnot2\n'),file='notattend_all.txt',sep='\t',append=T)}
}
}
########################################################################################
################## Merge all enrichment results ##################
attend_all=read.table('attend_all.txt',header=F,sep='\t')
attend_all=attend_all[,1]
i=1
datG=readRDS(paste0('GSEA/',attend_all[i],'_GSEA_all.Rds'))
#Screening significant results
dat_all=datG@result %>% filter(qvalues<0.01) 
dat_all$TRP=attend_all[i]
dat_all= dat_all %>% select(TRP,ID,enrichmentScore,NES,pvalue,p.adjust,qvalues)
for(i in 2:length(attend_all)){
datG=readRDS(paste0('GSEA/',attend_all[i],'_GSEA_all.Rds'))
dat=datG@result# %>% filter(pvalue<0.05)
if(nrow(dat)>0){
  dat$TRP=attend_all[i]
  dat= dat %>% select(TRP,ID,enrichmentScore,NES,pvalue,p.adjust,qvalues)
  
  dat_all=rbind(dat_all,dat)
}
}
colnames(dat_all)[2]='Pathway'
dat_all$Pathway=gsub('HALLMARK_','',dat_all$Pathway)
#Save result
#write.table(dat_all,'TRP_GSEA_result.txt',sep='\t',quote=F,row.names=F)
########################################################################################
################## Drawing GSEA diagram ##################
#Extract the required functions from other packages
gsInfo <- getFromNamespace("gsInfo", "enrichplot")
ggtable <- getFromNamespace("ggtable", "enrichplot")
tableGrob2 <- getFromNamespace("tableGrob2", "enrichplot")
gseaScores <- getFromNamespace("gseaScores", "DOSE")

MyGSEAplot=function(dat,top1,condition){
#Set parameters
ES_geom='line'
ES_geom <- match.arg(ES_geom, c("line", "dot"))
geneList <- position <- NULL

#dat:Data of enrichment results
#top1:The name of the gene set you want to draw
#condition:Positive(NES>0) or negative(NES<0) correlation
title=paste(top1,condition,sep='<-')
dat_top1=dat %>% filter(ID == top1)
needTRP=dat_top1$TRP
needTRP=as.character(needTRP)

#Add information about Subfamily
needTRPfamily=TRPdata %>% filter(Approved.symbol %in% needTRP)
needTRPfamily2=needTRPfamily$Family
names(needTRPfamily2)=needTRPfamily$Approved.symbol
needTRPfamily2=needTRPfamily2[as.character(needTRP)]
colortrp=c('TRPA'='#8470FF',
		  'TRPC'='#1E90FF',
		  'TRPM'='#3CB371',
		  'TRPML'='#FF8C00',
		  'TRPP'='#EE82EE',
		  'TRPV'='#FF6347')
colortrp=colortrp[as.character(needTRPfamily2)]

#Select the gene set you want to draw
gsdata=data.frame()
x=data.frame()
for(i in 1:length(needTRP)){
  geneSetID=paste0('HALLMARK_',top1)
  x0=readRDS(paste0('GSEA/',needTRP[i],'_GSEA_all.Rds'))
  x1=x0@result %>% filter(Description == geneSetID)
  x1$Description=needTRP[i]
  x1$Family=as.character(needTRPfamily2[needTRP[i]])
  x1$Color=colortrp[as.character(needTRPfamily2[needTRP[i]])] %>% as.character()
  x=rbind(x,x1)
  gsdata0 <- gsInfo(x0, geneSetID)
  gsdata0$Description=needTRP[i]
  gsdata = rbind(gsdata,gsdata0)
}

#Draw the result diagram
library(ggplot2)
base_size=10
#Draw the base axis
p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
  theme_classic(base_size) +
  theme(panel.grid.major = element_line(colour = "grey92"),
		panel.grid.minor = element_line(colour = "grey92"),
		panel.grid.major.y = element_blank(),
		panel.grid.minor.y = element_blank()) +
  scale_x_continuous(expand=c(0,0))

#Add linear results about "Running Enrichment Score"
es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description),
					 size=1)
p.res <- p + es_layer +
  theme( legend.title = element_blank(),
		 legend.background = element_rect(fill = "transparent"))
#Modify theme style
p.res <- p.res + ylab("Running Enrichment Score") +
  theme(axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.line.x=element_blank(),
		plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))

#Add linerange
i <- 0
for (term in unique(gsdata$Description)) {
  idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
  gsdata[idx, "ymin"] <- i
  gsdata[idx, "ymax"] <- i + 1
  i <- i + 1
}
p2 <- ggplot(gsdata, aes_(x = ~x)) +
  geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
  xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
  theme(legend.position = "none",
		plot.margin = margin(t=-.1, b=0,unit="cm"),
		axis.ticks = element_blank(),
		axis.text = element_blank(),
		axis.line.x = element_blank()) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))


#Add Ranked list metric results about "Running Enrichment Score"
df2 <- p$data
df2$y <- p$data$geneList[df2$x]
p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
						 color="grey")
p.pos <- p.pos + ylab("Ranked List Metric") +
  xlab("Rank in Ordered Dataset") +
  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))


#Add title
if (!is.null(title) && !is.na(title) && title != "")
  p.res <- p.res + ggtitle(title)

#Modify color(Assign colors according to TRP family)
color=x$Color
# if (length(color) == length(geneSetID)) {
if(length(color) == length(unique(gsdata$Description))) {
  p.res <- p.res + scale_color_manual(values=color)
  if (length(color) == 1) {
	 p.res <- p.res + theme(legend.position = "none")
	 p2 <- p2 + scale_color_manual(values = "black")
  } else {
	 p2 <- p2 + scale_color_manual(values = color)
  }
}

#subplots:indicates which graphs are drawn,1, 2 and 3 correspond to the top linear graph, the middle dense graph and the lower rank graph respectively
subplots=1:3
#Combine graph
plotlist <- list(p.res, p2, p.pos)[subplots]
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] +
  theme(axis.line.x = element_line(),
		axis.ticks.x=element_line(),
		axis.text.x = element_text())

if (length(subplots) == 1)
  return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2,
												  l=.2, unit="cm")))
#Set the layout of three diagrams
rel_heights=c(1.5, 0.5, 1)

if (length(rel_heights) > length(subplots))
  rel_heights <- rel_heights[subplots]

#Plot the graph
aa=aplot::plot_list(gglist = plotlist, ncol=1, heights=rel_heights)
return(aa)
}

#Read the merged result file of GSEA
datafile='TRP_GSEA_pathway_all.txt'
dat_1=read.table(datafile,sep='\t',header=T)
dat_1$TRP=factor(dat_1$TRP,levels=sort(TRPdata$Approved.symbol))
#Filter results with qvalues less than 0.01
dat_2=dat_1 %>% filter(qvalues < 0.01) %>% distinct()

#plot Positive result
condition='Positive'
dat_H=dat_2 %>% filter(NES >0)
#Extract all significant genesets
top=names(sort(table(dat_H$ID),decreasing = T))  %>% sort()
for(k in 1:length(top)){
top1=top[k]
plotH=MyGSEAplot(dat_H,top1,condition)
dir.create(paste0('GSEAPLOT/',top1))
pdf(paste0('GSEAPLOT/',top1,'/',top1,'_',condition,'.pdf'),height = 6,width=9)
print(plotH)
dev.off()
}

#plot Negative result
condition='Negative'
dat_L=dat_2 %>% filter(NES < 0)
top=names(sort(table(dat_L$ID),decreasing = T)) %>% sort()
for(k in 1:length(top)){
top1=top[k]
plotL=MyGSEAplot(dat_L,top1,condition)
dir.create(paste0('GSEAPLOT/',top1))
pdf(paste0('GSEAPLOT/',top1,'/',top1,'_',condition,'.pdf'),height = 6,width=9)
print(plotL)
dev.off()
}
########################################################################################
################## Statistical enrichment results ##################
#Count the number of positive and negative gene sets
uptable=as.data.frame(table(dat_2$TRP[which(dat_2$NES > 0)]))
downtable=as.data.frame(table(dat_2$TRP[which(dat_2$NES < 0)]))
colnames(uptable)=c('TRP','UP')
colnames(downtable)=c('TRP','DOWN')

numdata=merge(uptable,downtable,by='TRP',all=T)
numdata[is.na(numdata)]=0
rownames(numdata)=numdata$TRP
numdata=numdata[,-1]
numdata=reshape2::melt(as.matrix(numdata))
colnames(numdata)=c('TRP','STATUS','Value')
pbar_cancer <- ggplot(numdata, aes(x=TRP, y=Value,group=STATUS))+
   geom_bar(stat="identity", position="dodge", aes(fill=STATUS))+
   scale_fill_manual(values=c("#CB181D","#1F78B4"))+
   theme_classic()+theme(axis.text = element_text(angle = 90, size = 10),
         legend.position = 'bottom',
         legend.title = element_text(size=10),
         legend.text = element_text(size=10),
         axis.title = element_text(size=10),
         panel.background = element_blank()
         )+
   geom_text(aes(label=Value), position=position_dodge(1),size=2,colour='grey20')

pdf('Numbers_TRP_pathway.pdf')
pbar_cancer
dev.off()
########################################################################################
################## Preparation documents of Circos ##################
#-----------------karyotype.txt
#set the frame of heatmap and barplot
len_pathway=length(unique(dat_2$ID))
len_pathway
len_TRP=length(TRPdata$Approved.symbol)
len_TRP
karyotype_1=c('chr','-','hs1','1',0,len_TRP*4,'black')
karyotype_2=c('chr','-','hs2','2',0,len_pathway*3,'black')
karyotype=rbind(karyotype_1,karyotype_2)
write.table(karyotype,'File/data/karyotype.txt',sep='\t',row.names = F,col.names = F,quote=F)

#-----------------Preparation documents of heatmap
#Set the position of the display value of each gene
library(dplyr)
library(tidyverse)
dat_1=read.table(datafile,sep='\t',header=T)
dat_1=dat_1 %>% filter(qvalues < 0.01)
dat2=dat_1 %>% select(ID,TRP,NES) %>% arrange(TRP)

#Add Na values if not enriched on the pathway
dat3=dat2 %>% spread(ID,NES)
rownames(dat3)=dat3$TRP
dat3=dat3[,-1]
dat3=reshape2::melt(as.matrix(dat3))
colnames(dat3)=c('TRP','ID','NES')
mygene=as.character(unique(dat3$TRP)) %>% sort()
dat3=cbind('hs1',dat3)
kk2=cbind(dat3,rep(seq(0,len_TRP*4-4,by=4),len_pathway))
kk2=cbind(kk2,rep(seq(4,len_TRP*4,by=4),len_pathway))
colnames(kk2)[5:6]=c('start','end')
kk2[is.na(kk2)]=0
#output file
samplenames=as.character(unique(kk2$ID)) %>% sort()
samplenames=cbind(samplenames,paste0('Sample_',1:len_pathway))
samplenames=as.data.frame(samplenames)
apply(samplenames,1,function(a){
   kk0=kk2 %>% filter(ID == a[1])
   kk0=kk0[order(as.character(kk0$TRP)),]
   kk0=kk0[,c(1,5,6,4)]
   write.table(kk0,paste0('File/data/',a[2],'.txt'),sep='\t',row.names = F,col.names = F,quote=F)
   print(a[2])
})
#-----------------KO_orthology_text.txt
#Show TRP Symbols
KO_orthology_1=cbind('hs1',seq(0+2,len_TRP*4-2,4),seq(0+2,len_TRP*4-2,4),mygene)
TRPdata2=TRPdata %>% select(Approved.symbol,Family)
KO_orthology=merge(KO_orthology_1,TRPdata2,by.x='mygene',by.y = 'Approved.symbol')

colortrp=c('TRPC'='trpc',#1E90FF
           'TRPP'='trpp',#PKD2L #EE82EE
           'TRPM'='trpm',#3CB371
           'TRPV'='trpv',#FF6347
           'TRPML'= 'trpml',#MCOLN #FF8C00
           'TRPA'= 'trpa'#8470FF
)
KO_orthology=merge(KO_orthology,colortrp,by.x='Family',by.y='row.names',sort = F)
KO_orthology[,c(4,5)]=apply(KO_orthology[,c(4,5)],2,as.numeric)
KO_orthology=KO_orthology[order(KO_orthology$V2),]
KO_orthology2=KO_orthology[,c(3,4,5,2,6)]
KO_orthology2[,5]=paste0('color=',KO_orthology2[,5])
write.table(KO_orthology2,'File/data/KO_orthology_text.txt',sep='\t',row.names = F,col.names = F,quote=F)

#-----------------class_text.txt
#Show gene sets names
dat=read.table(datafile,header=T,sep='\t')
dat_2=dat %>% filter(qvalues < 0.01)
len_pathway=length(unique(dat_2$ID))

mypathway=unique(dat_2$ID) %>% sort()


class_text=cbind('hs2',seq(0+1,len_pathway*3-3+1,3),seq(0+1,len_pathway*3-3+1,3),1:length(mypathway))
class_text=as.data.frame(class_text)
write.table(class_text,'File/data/class_text.txt',
            sep='\t',quote=F,row.names = F,col.names = F)

class_text=cbind('hs2',seq(0+1,len_pathway*3-3+1,3),seq(0+1,len_pathway*3-3+1,3),mypathway)
class_text=as.data.frame(class_text)
write.table(class_text,'File/data/class_text_2.txt',
            sep='\t',quote=F,row.names = F,col.names = F)

#-----------------links.txt
#Documents linking TRP and gene sets
dat=read.table(datafile,
               header=T,sep='\t')
dat_2=dat %>% filter(qvalues < 0.01 & NES >0)

dat_class=read.table('File/data/class_text_2.txt',
                     header=F,sep='\t')
colnames(dat_class)=c('chr','start_pathway',
                      'end_pathway','ID')

KO_orthology2=read.table('File/data/KO_orthology_text.txt',
                         sep='\t',header = F)
colnames(KO_orthology2)=c('chr','start','end','TRP','Color')

KO_orthology3=merge(KO_orthology2,dat_2,by='TRP')
dat_class2=merge(dat_class,KO_orthology3,by='ID')
dat_class2=dat_class2 %>% select(chr.y,start,end,
                                 chr.x,start_pathway,end_pathway,NES
)

dat_class2$color='color=blues-13-seq-12_a3'#  blues-13-seq-12_a3
dat_class2[which(dat_class2$NES > 0 ),'color']='color=lred'
dat_class2=dat_class2[order(dat_class2$start),]
dat_class2=dat_class2 %>% select(-NES)
write.table(dat_class2,'File/data/links.txt',
            sep='\t',quote=F,row.names = F,col.names = F)

#-----------------light_border.txt
#Add tick marks
len_pathway=length(unique(dat_2$ID))
dat=data.frame(chr=rep('hs2',20),
               start=c(0,0,len_pathway*3,len_pathway*3,
                       len_pathway*3-1,len_pathway*3-1,len_pathway*3-1,len_pathway*3-1,
                       len_pathway*3-1,len_pathway*3-1,len_pathway*3-1,len_pathway*3-1,
                       0,0,0,0,
                       0,0,0,0),
               end=c(0,0,len_pathway*3,len_pathway*3,
                     len_pathway*3,len_pathway*3,len_pathway*3,len_pathway*3,
                     len_pathway*3,len_pathway*3,len_pathway*3,len_pathway*3,
                     1,1,1,1,
                     1,1,1,1)
)
colormy0='fill_color=black,'
colormy=paste0(colormy0,'r0=1.0r,r1=1.5r')
colormy=c(colormy,paste0(colormy0,'r0=1.55r,r1=2.05r'))
colormy=c(colormy,paste0(colormy0,'r0=1.0r,r1=1.5r'))
colormy=c(colormy,paste0(colormy0,'r0=1.55r,r1=2.05r'))
colormy=c(colormy,paste0(colormy0,'r0=1.1r,r1=1.1r'))
colormy=c(colormy,paste0(colormy0,'r0=1.2r,r1=1.2r'))
colormy=c(colormy,paste0(colormy0,'r0=1.3r,r1=1.3r'))
colormy=c(colormy,paste0(colormy0,'r0=1.4r,r1=1.4r'))
colormy=c(colormy,paste0(colormy0,'r0=1.65r,r1=1.65r'))
colormy=c(colormy,paste0(colormy0,'r0=1.75r,r1=1.75r'))
colormy=c(colormy,paste0(colormy0,'r0=1.85r,r1=1.85r'))
colormy=c(colormy,paste0(colormy0,'r0=1.95r,r1=1.95r'))
colormy=c(colormy,paste0(colormy0,'r0=1.1r,r1=1.1r'))
colormy=c(colormy,paste0(colormy0,'r0=1.2r,r1=1.2r'))
colormy=c(colormy,paste0(colormy0,'r0=1.3r,r1=1.3r'))
colormy=c(colormy,paste0(colormy0,'r0=1.4r,r1=1.4r'))
colormy=c(colormy,paste0(colormy0,'r0=1.65r,r1=1.65r'))
colormy=c(colormy,paste0(colormy0,'r0=1.75r,r1=1.75r'))
colormy=c(colormy,paste0(colormy0,'r0=1.85r,r1=1.85r'))
colormy=c(colormy,paste0(colormy0,'r0=1.95r,r1=1.95r'))

dat=cbind(dat,colormy)
write.table(dat,'File/data/light_border.txt',sep='\t',
            quote=F,row.names = F,col.names = F)

#-----------------hist_GSH.txt and hist_GSL.txt
#Counting the number of TRPs enriched in each gene set 
dat_H=dat_2 %>% filter(NES>0)
hist_H=as.data.frame(table(dat_H$ID))

dat_L=dat_2 %>% filter(NES<0)
hist_L=as.data.frame(table(dat_L$ID))

class_text=read.table('File/data/class_text_2.txt',
                      sep='\t',header = F)

colnames(class_text)=c('chr','start','end','ID')

hist_H=merge(class_text,hist_H,by.x='ID',by.y='Var1',all.x = T)
hist_H[is.na(hist_H)]=0
hist_H$end=hist_H$end+1
hist_H=hist_H %>% select(-ID)
write.table(hist_H,'File/data/hist_GSH.txt',
            sep='\t',quote = F,row.names = F,col.names =F)

hist_L=merge(class_text,hist_L,by.x='ID',by.y='Var1',all.x = T)
hist_L[is.na(hist_L)]=0
hist_L$end=hist_L$end+1
hist_L=hist_L %>% select(-ID)
write.table(hist_L,'File/data/hist_GSL.txt',
            sep='\t',quote = F,row.names = F,col.names =F)
#-----------------hist_3.txt
#Set the tick mark of barplot
start=seq(0,len_pathway*3-3,3)
end=seq(3,len_pathway*3,3)
dat0=cbind(start,end)
dat0=as.data.frame(dat0)
dat0=cbind('hs2',dat0,'0')
write.table(dat0,'File/data/hist_3.txt',sep='\t',
            row.names = F,col.names = F,quote = F)

#-----------------highlight.txt
#Set barplot background color
KO_orthology2=read.table('File/data/KO_orthology_text.txt',sep='\t',header = F)
colnames(KO_orthology2)=c('chr','start','end','mygene','color')
TRPdata2=TRPdata %>% select(Approved.symbol,Family)
KO_orthology2=merge(KO_orthology2,TRPdata2,by.x='mygene',by.y = 'Approved.symbol')
unique(KO_orthology2$Family)

highlight=KO_orthology2 %>% group_by(Family) %>% summarise(minn=min(start),maxn=max(end))
highlight$minn=highlight$minn-1
highlight$maxn=highlight$maxn+1
highlight=as.data.frame(highlight)
highlight$color='fill_color=greys-5-seq-2,r0=0.87r,r1=0.94r'
write.table(highlight,'File/data/highlight.txt',sep='\t',
            row.names = F,col.names = F,quote = F)

#-----------------run circos.conf
circos –conf etc/circos.conf

########################################################################################

