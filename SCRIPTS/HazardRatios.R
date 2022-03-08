### The distribution of hazard ratios across different cancer types
#################################################################################################
rm(list = ls())
library(survival)
cancer_names<-c('HNSC','GBM','LGG','UVM','LUAD','LUSC','MESO','BRCA','ACC','PCPG','THCA','COAD','ESCA','READ',
                'STAD','CHOL','LIHC','PAAD','BLCA','KICH','KIRC','KIRP','PRAD','TGCT','CESC','OV','UCEC','UCS',
                'DLBC','LAML','THYM','SKCM','SARC')
sur_mat_res<-c()
func_mid<-function(v){
  return(median(v))
}
gene_HR_res<-c()
for(i in 1:33){
  clinical_road<-"/clinical/"
  clinical_cancer<-read.table(paste0(clinical_road,cancer_names[i],".txt"),header = T,as.is = T,sep = "\t",
                              row.names = 1,check.names = F,quote = "")
  if(length(which(is.na(clinical_cancer$OS.time)|is.na(clinical_cancer$OS)))>0){
    clinical_cancer<-clinical_cancer[-which(is.na(clinical_cancer$OS.time)|is.na(clinical_cancer$OS)),]
  }
  expression_road<-"/TRP_FPKM/"
  cancer_mat<-read.table(paste0(expression_road,"TCGA_",cancer_names[i],".txt"),header = T,as.is = T,sep = "\t",
                         row.names = 1,check.names = F)
  tumor_sample<-read.table(paste0("/clinical/",
                                  cancer_names[i],"_tumor.txt"),header=F,sep="\t")
  cancer_mat<-cancer_mat[,tumor_sample$V1]
  inter_sample<-intersect(colnames(cancer_mat),rownames(clinical_cancer))
  cancer_mat<-cancer_mat[,inter_sample]
  index0<-which(apply(cancer_mat,1,func_mid)==0)
  if(length(index0!=0)){
    cancer_mat<-cancer_mat[-index0,]
  }
  clinical_cancer<-clinical_cancer[inter_sample,]
  cancer_mat<-as.matrix(cancer_mat)
  clinical_cancer$sample<-rownames(clinical_cancer)
  gene_TRP<-"TRPA1"#TRPM2 for Figures5A, TRPA1 for Figures5B, Figures5C
  gene_exp<-cancer_mat[gene_TRP,]
  clinical_cancer$FPKM<-gene_exp
  clinical_cancer$OS.time<-as.numeric(clinical_cancer$OS.time)
  clinical_cancer$FPKM<-as.numeric(clinical_cancer$FPKM)
  res.cut <- surv_cutpoint(clinical_cancer, time = "OS.time", event = "OS",variables = "FPKM")
  res.cat <- surv_categorize(res.cut)
  res.cat$FPKM<-factor(res.cat$FPKM,levels=c("low","high"))
  data.survdiff <- survdiff(Surv(OS.time,OS)~FPKM,data = res.cat)
  ###Calculate P, HR, and confidence intervals
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 
  gene_HR<-c(gene_TRP,cancer_names[i],p.val,HR,up95,low95)
  gene_HR_res<-rbind(gene_HR_res,gene_HR)
  cat(i ,sep = "\n")
}
######
gene_HR_res<-as.data.frame(gene_HR_res)
colnames(gene_HR_res)<-c("Symbol","cancer","pvalue","HR","up95","low95")
gene_HR_res$pvalue<-as.numeric(gene_HR_res$pvalue)
gene_HR_res$HR<-as.numeric(gene_HR_res$HR)
gene_HR_res$low95<-as.numeric(gene_HR_res$low95)
gene_HR_res$up95<-as.numeric(gene_HR_res$up95)
####################Assign 0 to cancers with large confidence intervals, Different TRP gene will be different
gene_HR_res[which(gene_HR_res$cancer=="THCA"),]$up95<-0
gene_HR_res[which(gene_HR_res$cancer=="THCA"),]$low95<-0
gene_HR_res[which(gene_HR_res$cancer=="THYM"),]$up95<-0
gene_HR_res[which(gene_HR_res$cancer=="THYM"),]$low95<-0
gene_HR_res[which(gene_HR_res$cancer=="PCPG"),]$up95<-0
gene_HR_res[which(gene_HR_res$cancer=="PCPG"),]$low95<-0

gene_HR_res$cancer<-factor(gene_HR_res$cancer,levels = rev(cancer_names))
func_round<-function(v){return(round(v,3))}####Reserve three significant digits
#func_round(c(0.123456,0.4567889))
gene_HR_res$pvalue<-func_round(gene_HR_res$pvalue)
pdf("HR.pdf",width = 4,height=8)
ggplot(data=gene_HR_res, aes(y=cancer,color='red'),ylab=NA) + 
  geom_errorbar(data = gene_HR_res,aes(xmin = low95, xmax=up95, ylab=NA), #The error bar represents a 95% confidence interval
                width=0.2,#The width of the short horizontal line at the end of the error bar
                position=position_dodge(0.9), 
                alpha = 1,xlim=c(0,3),
                size=1)  +theme(axis.text.y = element_blank())+
  theme_bw()+
  #Plot the median as a dot plot
  geom_point(data = gene_HR_res,aes(x=HR, y=cancer),pch=19,position=position_dodge(0.9),size=4,ylab=NA,
             ylab=gene_HR_res$pvalue) +
  
  theme(
    panel.grid.major = element_blank(),   #Grid lines are not displayed
    panel.grid.minor = element_blank())+  #Grid lines are not displayed
  geom_vline(xintercept=1, linetype="dotted") #Adds a dotted line at the specified position
dev.off()
