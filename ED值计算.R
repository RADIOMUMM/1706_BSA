install.packages("vcfR")
library("vcfR")
library("binmapr")

vcf <- read.vcfR(choose.files())
gt=extract.gt(vcf)
ad=extract.gt(vcf,"AD")

mask <- which(gt[,"R01"] == "0/0" &  gt[,"R13"] == "1/1")

ad_flt <- ad[mask,c("R53", "R55")]

ED_list <- apply(ad_flt, 1, function(x){
  count <- as.numeric(unlist(strsplit(x, ",",fixed = TRUE,useBytes = TRUE)))
  depth1 <- count[1] + count[2]
  depth2 <- count[3] + count[4]
  
  ED <- sqrt((count[3] / depth2 - count[1] / depth1)^2 + 
               (count[4] / depth2- count[2] /depth1)^2)
  return(ED^5)
  
})

par(mfrow = c(3,4))

for (i in paste0("Arahy.", formatC(1:20, width = 2, flag=0)) ){
  ED_flt <- ED_list[grepl(i,names(ED_list))]
  pos <- as.numeric(substring(names(ED_flt), 10))
  plot(pos, ED_flt,
       pch = 20, cex = 0.2,
       xlab = i,
       ylab = "ED")
}
