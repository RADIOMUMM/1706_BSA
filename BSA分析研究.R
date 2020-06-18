install.packages("ggplot2")
install.packages("devtools")
install.packages("farver"
                 )
library(usethis)
library(devtools)
devtools::install_local("D:/迅雷下载/QTLseqr-master.zip")
 library(QTLseqr)
library(ggplot2)
install.packages("vcfR")
library("vcfR")
devtools:: install_local("D:/迅雷下载/Yang2013data-master.zip")

library("Yang2013data")
#导入数据
rawData <- system.file(
  "extdata", 
  "Yang_et_al_2013.table",
  package = "Yang2013data",
  mustWork = TRUE)

#然后我们分别给高低池命名。写染色体编号。
HighBulk <-"SRR834931"

LowBulk <-"SRR834927"

Chroms <-paste0(rep("Chr",12),1: 12)


#snp数据读取
df <-
  importFromGATK(
    file = rawData,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  )
#对snp进行过滤
df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.20,
    minTotalDepth = 100,
    maxTotalDepth = 400,
    minSampleDepth = 5,#单个样品最小测序深度
    minGQ = 99
  )

#进行G运算
df_filt <- runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = 1e6,
  outlierFilter = "deltaSNP")

#Run QTLseq analysis
df_filt <- runQTLseqAnalysis(
  SNPset = df_filt,
  windowSize = 500000,#窗口大小
  popStruc = "F2",#群体类型
  bulkSize = c(25, 25),
  replications = 10000,
  intervals = c(95, 99)
)


par(mfrow=c(1,2))
#Plot
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)

#export summary CSV
getQTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")
