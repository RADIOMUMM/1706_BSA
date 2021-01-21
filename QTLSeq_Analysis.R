library(QTLseqr)
library("ggplot2")


HighBulk <- "H1706Z" ## 高值混池名�?
LowBulk <- "H1706F"  ## 低值混池名�?
Chroms <-paste0("arahy.Tifrunner.gnm2.Arahy.",c("01","02","03","04","05","06","07","08","09","10",
                                                "11","12","13","14","15","16","17","18","19","20"))
## SNP数据读取
df <- importFromGATK(file = choose.files() ,
    highBulk = HighBulk, # 指定高值混�?
    lowBulk = LowBulk, # 指定低至混池
    chromList = Chroms # 指定分析用染色体列表
    ) 


## 绘制深度分布�? 

p1 <- ggplot(data = df) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
  xlim(0,800)
  
pdf(file="SNP_depth.pdf")
p1
dev.off()

## 绘制ref等位基因频率分布�?
p2 <- ggplot(data = df) +
  geom_histogram(aes(x = REF_FRQ))
  
pdf(file="ref_allele_frequency.pdf")
p2
dev.off()


## 绘制高值混池SNP-index分布
p3 <- ggplot(data = df) +
  geom_histogram(aes(x = SNPindex.HIGH))
pdf(file="SNPindex.HIGH.dis.pdf")
p3
dev.off()

## 绘制低值混池SNP-index分布
p4 <- ggplot(data = df) +
  geom_histogram(aes(x = SNPindex.LOW))

pdf(file="SNPindex.LOW.dis.pdf")
p4
dev.off()

## SNP过滤
df <- subset(df, !is.na(SNPindex.LOW) & !is.na(SNPindex.HIGH))
df_filt <- filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.20, ## ref allele频率过滤�?0.2表示0.2~0.8之间
    minTotalDepth = 100, ## 最小深度过�?
    maxTotalDepth = 400, ## 最大深度过�?
    minSampleDepth = 5, ## 单个样品最小深�?
    minGQ = 99, ## genotype quality 过滤
    verbose = TRUE  ## 输出日志
  )

## 进行deltaSNPindex 计算
df_filt <- runQTLseqAnalysis(df_filt,
    windowSize = 3e6, ## 窗口大小
    popStruc = "F2",  ## 群体类型，F2 或者RIL
    bulkSize = c(25,25), ## 混样个数，第一个高值组
    replications = 10000, ## bootstrap次数
    intervals = c(95, 99) ## 置信区间
)


## SNP 沿染色体分布�?
p5 <- plotQTLStats(SNPset = df_filt, var = "nSNPs")
pdf(file="SNP_filter.window.pdf", width = 20, height = 4)
p5
dev.off()

## deltaSNPindex沿染色体分布
p6 <- plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
pdf(file="deltaSNPindex.pdf", width = 20, height = 4)
p6
dev.off()

## 提取显著性区�?
QTL <- getSigRegions(SNPset = df_filt, 
    method = "QTLseq", interval = 95)

## 输出到文�?
write.table(QTL[[1]],
    "QTLseq_result.SigRegions.table",
    sep="\t",quote=F)

results95 <- getQTLTable(
  SNPset = df_filt,  
  method = "QTLseq", 
  interval = 95,  ## 结果阈�?
  export = FALSE)

## 输出到文�?
write.table(results95,
            file="QTLseq_result_deltaSNPindex_CI95.table" ,
            sep="\t", quote = F)


## 提取QTLseq结果
results99 <- getQTLTable(
  SNPset = df_filt,  
  method = "QTLseq", 
  interval = 99,  ## 结果阈�?
  export = FALSE)

## 输出到文�?
write.table(results99,
    file="QTLseq_result_deltaSNPindex_CI99.table" ,
    sep="\t", quote = F)

