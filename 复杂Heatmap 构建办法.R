library(devtools)
library(usethis)
devtools::install_github("jokergoo/ComplexHeatmap")
library(grid)

library(ComplexHeatmap)
a<-read.delim("clipboard",header=T,row.names=1)
x=as.matrix(a)
#自设置颜色
library(circlize)

#自定义颜色
mycol <- colorRamp2(c( 0.01, 0.05, 0.15,0.2), c("red","skyblue3","lightskyblue1","lightskyblue2"))
#自定义图例
heatmap_legend_param = list(
  title= "legend", title_position = "topcenter", 
  legend_height=unit(8,"cm"), legend_direction="vertical")



Heatmap(x,cluster_rows = FALSE,cluster_columns = FALSE,col =mycol
        , heatmap_width = unit(10, "cm"), heatmap_height = unit(17, "cm")
        ,heatmap_legend_param = list( title= "P_value", 
                                      title_position = "topcenter"
                                      , legend_height=unit(5,"cm")
                                      , legend_direction="vertical" ))

#其他参数设置参考
#https://jokergoo.github.io/ComplexHeatmap-reference/book/more-examples.html#visualize-cell-heterogeneity-from-single-cell-rnaseq
