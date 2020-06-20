install.packages("devtools")
install.packages("vcfR")
library(devtools)
devtools::install_github("xuzhougeng/binmapr")
library("vcfR")
library("binmapr")

vcf <- read.vcfR(choose.files())
gt=extract.gt(vcf)
ad=extract.gt(vcf,"AD")
head(gt)
head(ad)

#R01为“1/1” R13为“0/0” 的情况下
mask=which(gt[,"R01"] == "1/1"& gt[,"R13"] == "0/0")
ad_flt =ad[mask,c("R53","R55")]
colnames(ad_flt) <- c("T_Pool","S_Pool")
freq=calcFreqFromAd(ad_flt,min.depth = 8,max.depth = 100)
freq2=freq[Matrix::rowSums(is.na(freq)) == 0, ]


par(mfrow = c(4,5))

for (i in paste0("Arahy.", formatC(1:20, width = 2, flag=0)) ){
  freq_flt <- freq2[grepl(i,row.names(freq2)), ]
  pos <- as.numeric(substring(row.names(freq_flt), 10))
  plot(pos, freq_flt[,1] - freq_flt[,2], ylim = c(-1,1),
       pch = 20, cex = 0.2,
       xlab = i,
       ylab = expression(paste(Delta, " " ,"SNP index")))
}



