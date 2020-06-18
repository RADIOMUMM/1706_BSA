install.packages("qqman")
library(qqman)
data("gwasResults")
head(gwasResults)
#第一列为SNP的名字，第二列CHR为所在染色体，第三列BP为染色体上所在位置。
#要注意如果你的CHR中存在X，Y这样的，
#需要给他们转化为数字如赋予23，24等，
#其中第一列SNP的名字是可选择的，后三列是必须提供的。
manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P" )
install.packages("viridis")
library(viridis)

manhattan(gwasResults, annotatePval = 0.01,annotateTop = T, col=c("orange","green4"))
