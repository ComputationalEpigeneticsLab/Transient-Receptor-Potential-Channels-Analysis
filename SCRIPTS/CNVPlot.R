### Heat map showing the frequency of CNV alterations for TRPs across cancer types. 
### The upper of each rectangle represents the frequency of CNV loss and the bottom shows the frequency for CNV gain. 

library(ComplexHeatmap)
library(circlize)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
cnv_number <- read.table("TCGA.ACC.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes_Ratio.txt",sep="\t",header=T,stringsAsFactor=F)
cnv_number_gain <- cnv_number[,c(1,2,5)]
cnv_number_loss <- cnv_number[,c(1,2,7)] 
cnv_number_gain1 <- dcast(data=cnv_number_gain,gene ~ TCGA)
cnv_number_loss1 <- dcast(data=cnv_number_loss,gene ~ TCGA)
cnv_number_gain2 <- cnv_number_gain1[,-1]
rownames(cnv_number_gain2) <- cnv_number_gain1[,1]
cnv_number_gain2 <- cnv_number_gain2[-9,TCGAoreder[,1]]
cnv_number_loss2 <- cnv_number_loss1[,-1]
rownames(cnv_number_loss2) <- cnv_number_loss1[,1]
cnv_number_loss2 <- cnv_number_loss2[-9,TCGAoreder[,1]]
UpColor <- colorRamp2(breaks = c(0, 1), colors = c("white","#AB221F"))
DnColor <- colorRamp2(breaks = c(0, 1), colors = c("white","#3878C1"))
Heatmap(cnv_number_gain2,
        column_title = "Copy number variation across cancer types", 
        rect_gp = gpar(type = "none"),  
        show_heatmap_legend = F, 
        cluster_rows = F, cluster_columns = F
        )
row_an <-  HeatmapAnnotation(type = c(rep("TRPML",3),rep("TRPP",3),rep("TRPA",1),rep("TRPC",6),rep("TRPM",8),rep("TRPV",6)), 
                                show_annotation_name = F,
                                col = list(type = c(TRPML="#ff8c00",TRPP= "#ee82ee",TRPA= "#8470ff",TRPC= "#1e90ff",TRPM= "#3cb371",TRPV= "#ff6347")), 
                                show_legend = T, 
                                annotation_legend_param = list(title = "Subfamilies",nrow = 1), 
                                which = "row" 
                             )
DiagFunc <- function(up, down){
    function(j, i, x, y, width, height, fill){
        grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                         unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                         gp = gpar(fill = DnColor(down[i, j]), col = "grey"))
        grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                         unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                         gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
    }
}
p1 <- Heatmap(cnv_number_gain2, column_title = "Copy number variation across cancer types",
              rect_gp = gpar(type = "none"),
              show_heatmap_legend = F,
              cluster_rows = F,
              cluster_columns = F, 
              left_annotation = row_an, 
              cell_fun = DiagFunc(up = cnv_number_gain2, down = cnv_number_loss2) 
              ); p1

lgd <- list(Legend(title = "CNV Gain", 
              col_fun = UpColor, 
              at = c(0,0.5,1), 
              direction = "horizontal" 
              ),
            Legend(title = "CNV Loss", 
              col_fun = DnColor, 
              at = c(0,0.5,1), 
              direction = "horizontal" 
              ) )
draw(p1, annotation_legend_list = lgd,
  annotation_legend_side = "bottom",
  heatmap_legend_side = "bottom",
  merge_legend = TRUE)
