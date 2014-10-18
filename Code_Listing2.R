# 6.2	Code Listing 2: R code to generate Figure 8, the copy number distribution of CNVs of CEU and Melanoma
# Copyright Sebastian Burgstaller-Muehlbacher, 2014. Code is distributed under GPLv3 or later, without warranty

library("sm")

getCNVSizes <- function(CNVs, sizeFilter){
  cnvSizes <- matrix(data=0, ncol=1)
  for(i in 1:nrow(CNVs)){
    if(CNVs[i,2] != 24 & CNVs[i,2] != 25 & sizeFilter == 0 & (CNVs[i,5]<(-0.1) | CNVs[i,5]>0.1)){
      cnvSizes <- rbind(cnvSizes, CNVs[i,4] - CNVs[i,3])
    
    }else if(CNVs[i,2] != 24 & CNVs[i,2] != 25 & sizeFilter > 0 & (CNVs[i,5]<(-0.1) | CNVs[i,5]>0.1)){
      size <- CNVs[i,4] - CNVs[i,3]
      
      if(size <= sizeFilter)
        cnvSizes <- rbind(cnvSizes, CNVs[i,4] - CNVs[i,3])
    }
  }
  return(cnvSizes)
}

#open CEU CNV file and melanoma CNV file

#first, do it unfiltered
abb_CEU <- read.table("regions_CEU.xls",stringsAsFactors=F, header=T)
abb_MEL <- read.table("myCSV_regions.xls",stringsAsFactors=F, header=T)

maxCNVSize <- 1000000

sizeDistrCEU <- getCNVSizes(abb_CEU,maxCNVSize)
sizeDistrMEL <- getCNVSizes(abb_MEL,maxCNVSize)

tmp1<- cbind(sizeDistrCEU, sizeDistrCEU)
tmp1[,1] <- 1

tmp2 <- cbind(sizeDistrMEL, sizeDistrMEL)
tmp2[,1] <- 2

cnv_dens <- rbind(tmp1,tmp2)
sm.density.compare(cnv_dens[,2],cnv_dens[,1] ,  positive=T, xlab="Size in bp")
title(main="CNV Size Distribution CEU vs Melanoma")

grp.f <- factor(cnv_dens[,1], levels=c(1,2), labels=c("CEU","Melanoma"))

colfill<-c(2:(2+length(levels(grp.f))))
legend(locator(1), levels(grp.f), fill=colfill)
