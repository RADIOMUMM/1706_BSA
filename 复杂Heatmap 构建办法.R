library(devtools)
library(usethis)
devtools::install_github("jokergoo/ComplexHeatmap")
library(grid)

library(ComplexHeatmap)
a<-read.delim("clipboard",header=T,row.names=1)
x=as.matrix(a)
X=log10(x+1)
#自设置颜色
library(circlize)

#自定义颜色
mycol <- colorRamp2(seq(0.01, 2, length=5), c("royalblue","lightsteelblue","khaki", "orange","red"), space = "RGB")
#自定义图例
heatmap_legend_param = list(
  title= "legend", title_position = "topcenter", 
  legend_height=unit(8,"cm"), legend_direction="vertical")



Heatmap(X,cluster_rows = T,cluster_columns = T,col =mycol
        , heatmap_width = unit(18, "cm"), heatmap_height = unit(10, "cm")
        ,heatmap_legend_param = list( title= "P_value", 
                                      title_position = "topcenter"
                                      , legend_height=unit(5,"cm")
                                      , legend_direction="vertical" ))


#其他参数设置参考
#https://jokergoo.github.io/ComplexHeatmap-reference/book/more-examples.html#visualize-cell-heterogeneity-from-single-cell-rnaseq
