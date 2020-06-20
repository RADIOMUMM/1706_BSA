library(pheatmap)
mydata<-read.csv(file.choose())
a=as.matrix(mydata)
pheatmap(a)
pheatmap(a, cluster_row = FALSE, 
         legend_breaks = -1:4, legend_labels = c("0", "1e-4", "1e-3", "1e-2", "1e-1", "1"))
#默认参数下是对行列均进行聚类（可设置cluster_row = FALSE, cluster_col = FALSE不进行行列的聚类；
#如果进行聚类了，还可以通过设置treeheight_row=0, treeheight_col=0不显示dendrogram），
#矩阵没有进行标准化（标准化参数为scale，可选"none", "row", "column"），
#热图的每个小块之间以灰色隔开（参数border_color，如果不想要border可以设置为NA，当然也可以设置成其它颜色）
#legend显示在右上方（可设置legend = FALSE不显示legend）；热图的颜色可利用参数color调整


# 在热图格子里展示文本
pheatmap(a, display_numbers = TRUE)
pheatmap(a, display_numbers = TRUE, number_format = "\%.1e")
pheatmap(a, display_numbers = matrix(ifelse(a > 5, "*", ""), nrow(a)))#还可以自己设定要显示的内容；
#可设置参数display_numbers将数值显示在热图的格子中，
#可通过number_format设置数值的格式，
#较常用的有"%.2f"（保留小数点后两位），"%.1e"（科学计数法显示，保留小数点后一位），
#number_color设置显示内容的颜色：


#例子

library(pheatmap)
install.packages("viridis")
library(viridis)
# viridis 黄绿色 magma 黄紫色 plasma 黄蓝色 inferno黄黑色 cividis黄灰色，数字是颜色广度

library(ggplot2)
a<-read.delim("clipboard",header=T,row.names=1)
#识别剪切板，创建矩阵
pheatmap(a)
p=pheatmap(a, color = plasma(6),cluster_row = T, display_numbers = TRUE,fontsize = 6.8,
         legend_breaks = -1:4,
         legend_labels = c("0", "0", "1", "2", "3", "4"))#图例中的刻度修改
#默认参数下是对行列均进行聚类（可设置cluster_row = FALSE, cluster_col = FALSE不进行行列的聚类；
#如果进行聚类了，还可以通过设置treeheight_row=0, treeheight_col=0不显示dendrogram），
#矩阵没有进行标准化（标准化参数为scale，可选"none", "row", "column"），
#热图的每个小块之间以灰色隔开（参数border_color，如果不想要border可以设置为NA，当然也可以设置成其它颜色）
#legend显示在右上方（可设置legend = FALSE不显示legend）；热图的颜色可利用参数color调整





