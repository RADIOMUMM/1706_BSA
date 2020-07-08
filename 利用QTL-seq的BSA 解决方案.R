HighBulk <-"R55"
LowBulk <-"R53"
Chroms1 <-paste0("Arahy.",c("01","02","03","04","05","06","07","08","09","10",
                            "11","12","13","14","15","16","17","18","19","20"))
library(QTLseqr)
library(vcfR)
vcf <- read.vcfR("C:/Users/Mumm/Desktop/America_Peanut.raw.filter.snp.anno.gatk.vcf")

chrom <- getCHROM(vcf)
pos <- getPOS(vcf)
ref <- getREF(vcf)
alt <- getALT(vcf)

ad <- extract.gt(vcf, "AD")
ref_split <- masplit(ad, record = 1, sort = 0)
alt_split <- masplit(ad, record = 2, sort = 0)
gt <- extract.gt(vcf, "GT")
dp= extract.gt(vcf,"DP")

df <- data.frame(CHROM = chrom,
                 POS = pos,
                 REF = ref,
                 ALT = alt,
                 AD_REF.R53 = ref_split[,3],
                 AD_ALT.R53 = alt_split[,3],
                 DP.R53=dp[,3],
                 AD_REF.R55 = ref_split[,4],
                 AD_ALT.R55 = alt_split[,4],
                 DP.R55=dp[,4]
               
)
mask <- which(gt[,"R01"] == "1/1" &  gt[,"R13"] == "0/0")
df <- df[mask,]
write.table(df, file = "rice.tsv", sep = "\t", row.names = F, quote = F)


df <- importFromTable("rice.tsv",
                      highBulk = HighBulk,
                      lowBulk = LowBulk,
                      chromList = Chroms1,
                      sep = "\t")

df <- subset(df, !is.na(SNPindex.LOW) & !is.na(SNPindex.HIGH))
#绘制深度分布图
library(ggplot2)
ggplot(data = df)+geom_histogram(aes(x=DP.HIGH+DP.LOW))+xlim(0,250)

ggplot(data = df)+geom_histogram(aes(x=REF_FRQ))

#开始分析
df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.20,
    minTotalDepth = 40,
    maxTotalDepth = 400,
    minSampleDepth = 5,#单个样品最小测序深度
    minGQ = 99
  )
#进行G分析
df_filt <- runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = 1e7/2.5,
  outlierFilter = "deltaSNP",
  filterThreshold =0.1)
#Run QTLseq analysis
df_filt <- runQTLseqAnalysis(
  SNPset = df_filt,
  windowSize = 1e7/2.5,
  popStruc = "F2",
  bulkSize = c(25, 25),
  replications = 10000,
  intervals = c(95, 99)
)
#SNP沿染色体分布图
plotQTLStats(SNPset = df_filt,var = "nSNPs")

plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)

plotQTLStats(SNPset =df_filt,var ="Gprime",plotThreshold =TRUE,q =0.01)




