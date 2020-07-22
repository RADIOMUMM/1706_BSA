library(vcfR)
vcf <- read.vcfR("C:/Users/Mumm/Desktop/R07_H1314WvsR06_H1314T.mutmap.vcf")

chrom <- getCHROM(vcf)
pos <- getPOS(vcf)
ref <- getREF(vcf)
alt <- getALT(vcf)

ad <- extract.gt(vcf, "AD")
ref_split <- masplit(ad, record = 1, sort = 0)
alt_split <- masplit(ad, record = 2, sort = 0)
gt <- extract.gt(vcf, "GT")
dp= extract.gt(vcf,"DP")
pl= extract.gt(vcf,"PL")

df= data.frame(CHROM = chrom,
                POS = pos,
                REF = ref,
                ALT = alt,
                GT.R06=gt[,1],
               AD.R06= ad[,1],
               PL.R06=pl[,1],
                DP.R06=dp[,1],
               GT.R07 =gt[,2],
               AD.R07 = ad[,2],
                DP.R07=dp[,2],
               PL.RO7=pl[,2])

m1=which(gt[,"R06"]!=gt[,"R07"])
DF=df[m1,]
write.table(DF, file = "R07_H1314WvsR06_H1314T.table", sep = "\t", row.names = F, quote = F)
