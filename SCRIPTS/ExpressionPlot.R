### Expression of TRPM8 across cancers. 
### Expression of TRPA1 across cancers. 
### Expression of TPRM3 across cancers.

rm(list = ls())
library(ggplot2)
library(ggpubr)
library(ggsignif)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),          
                     median = median   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     quantile_25 = quantile(xx[[col]],0.25),
                     quantile_75 = quantile(xx[[col]],0.75)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  #datac <- rename(datac, c("median" = measurevar)) 
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

cancer_names<-c('HNSC',"GBM",'LUAD','LUSC','BRCA','THCA','COAD','ESCA','READ','STAD','CHOL','LIHC',
                'BLCA','KICH','KIRC','KIRP','PRAD','UCEC')
gruop_TRP<-read.table("TRPList.txt",header = T,as.is = T,sep = "\t",check.names=FALSE)
exp_mat_res<-c()
for(i in 1:18){

  setwd("/TRP_FPKM")
  cancer_mat<-read.table(paste0("TCGA_",cancer_names[i],".txt"),row.names = 1,header = T,as.is = T,sep="\t",check.names=FALSE)
  cancer_mat<-cancer_mat["TRPM2",]#TRPM8 for figures3A, TRPA1 for figures3B, TRPM3 for figures3C
  normal_sample<-read.table(paste0("/TRP_FPKM/",
                                   cancer_names[i],"_normal.txt"),header=F,sep="\t")
  normal_sample$V1 <- substr(normal_sample$V1,1,15)
  normal_sample<-unique(normal_sample$V1)
  tumor_sample<-read.table(paste0("/TRP_FPKM/",
                                  cancer_names[i],"_tumor.txt"),header=F,sep="\t")
  tumor_sample$V1 <- substr(tumor_sample$V1,1,15)
  tumor_sample<-unique(tumor_sample$V1)
  cancer_mat<-cancer_mat[,c(normal_sample,tumor_sample)]
  n<-length(normal_sample)
  m<-length(tumor_sample)
  normal_exp<-cancer_mat[,1:n]
  tumor_exp<-cancer_mat[,(n+1:m)]
  normal_exp_res<-cbind(rep(cancer_names[i],n),rep("normal",n),as.matrix(normal_exp)[1,])
  tumor_exp_res<-cbind(rep(cancer_names[i],m),rep("tumor",m),as.matrix(tumor_exp)[1,])
  exp_mat<-rbind(normal_exp_res,tumor_exp_res)
  exp_mat_res<-rbind(exp_mat_res,exp_mat)

}
mean(as.numeric(normal_exp[1,]))
mean(as.numeric(cancer_mat[1,]))
colnames(exp_mat_res)<-c("cancer","label","expression")
exp_mat_res<-as.data.frame(exp_mat_res)
exp_mat_res$expression<-as.numeric(exp_mat_res$expression)
exp_mat_res$log2exp<-log2(exp_mat_res$expression+1)
###############################################################################################
Data_summary <- summarySE(exp_mat_res, measurevar="log2exp", groupvars=c("cancer","label"))
colnames(Data_summary)[7:8]<-c('quantile_25','quantile_75')
Data_summary$cancer<-factor(Data_summary$cancer,levels=cancer_names)
setwd("E:\\work\\project\\TRP_pancancer\\result")
pdf("TRPM2_FPKM.pdf",height = 4,width = 10)#TRPM8 for figures3A, TRPA1 for figures3B, TRPM3 for figures3C
ggplot(data=Data_summary, aes(x=cancer,color=label)) + 
    geom_errorbar(data = Data_summary,aes(ymin = quantile_25, ymax=quantile_75,group=label), #The error bar represents a 95% confidence interval
                  width=0.2, #The width of the short horizontal line at the end of the error bar
                  position=position_dodge(0.9), 
                  alpha = 1,
                  size=1) + 
    theme_bw()+
    #Plot the median as a dot plot
    geom_point(data = Data_summary,aes(x=cancer, y=median),pch=19,position=position_dodge(0.9),size=4) +
    scale_color_manual(values = c('#4393C3','#D6604D') )+ 
    #Fill color : c("#56B4E9", "#E69F00") c("#1F78B4","#CB181D") c('#4393C3','#D6604D') 
    theme(
      panel.grid.major = element_blank(),   #Grid lines are not displayed
      panel.grid.minor = element_blank())+  #Grid lines are not displayed
    xlab("Cancer")+ylab("log2(FPKM+1) (TRPM2)")+ #title of X&Y,TRPM8 for figures3A, TRPA1 for figures3B, TRPM3 for figures3C
    geom_vline(xintercept=c(seq(1.5,17.5,by=1)), linetype="dotted") #Adds a dotted line at the specified position
dev.off()