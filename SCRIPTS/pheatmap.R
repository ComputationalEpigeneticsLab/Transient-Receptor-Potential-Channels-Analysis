### The mutation frequency profile of TRPs across cancer types. 

library(data.table)
library(pheatmap)
TCGAoreder <- read.table("tumorOrder.txt",sep="\t",header=T,stringsAsFactor=F)
mutation_number <- read.table("MutationRatio.txt",sep="\t",header=T,stringsAsFactor=F)
mutation_number1 <- mutation_number[,c(1,2,5)] 
mutation_matrix <- dcast(data=mutation_number1,gene ~ TCGA)
mutation_matrix1 <- mutation_matrix[,-1]
rownames(mutation_matrix1) <- mutation_matrix[,1]
mutation_matrix1 <- mutation_matrix1[-9,TCGAoreder[,1]]
mutation_matrix2 <- mutation_matrix1
mutation_matrix2[mutation_matrix2>0.05] <- 0.05

annotation_row <- data.frame(Subfamilies = factor(c(rep("TRPML",3),rep("TRPP",3),rep("TRPA",1),rep("TRPC",6),rep("TRPM",8),rep("TRPV",6))))
ann_colors <- list(Subfamilies = c(TRPML="#ff8c00",TRPP= "#ee82ee",TRPA= "#8470ff",TRPC= "#1e90ff",TRPM= "#3cb371",TRPV= "#ff6347"))

### mutaiton pheatmap
rownames(annotation_row) <- rownames(mutation_matrix2)
mutaiton_pic <- pheatmap(mutation_matrix2, show_colnames= T, show_rownames= T, scale= "none", fontsize= 6.5,
         cluster_rows = F,cluster_cols = F,	 
        annotation_colors= ann_colors,
		  annotation_row= annotation_row, 
        col = colorRampPalette(c("white","steelblue"))(10)
)