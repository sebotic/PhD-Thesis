# Code Listing 4: CNVs and epigenetics. Annotates CNVs with regulatory elements and TFBS.
# Copyright Sebastian Burgstaller-Muehlbacher, 2014. Distributed under GPLv3 or later, without warranty.
# evaluation of CNVs for epigenetic marks and enhancers from ENCODE

library(xlsx)

new.cnvs <- read.xlsx("new_cnvs.xlsx", sheetIndex=1, stringsAsFactors=F, header=T)

genome.segmentations <- read.table("wgEncodeAwgSegmentationCombinedHelas3.bed", stringsAsFactors=F)


for(i in 1:nrow(new.cnvs)){
  
  new.cnvs[i, 10] <- paste0(searchRegulatoryElements(new.cnvs[i, 5], new.cnvs[i, 8], new.cnvs[i, 9])[,4], collapse=" ")
  
}
rm(genome.segmentations)



tfbs.161 <- read.table("wgEncodeRegTfbsClusteredV3.bed", stringsAsFactors=F)

for(i in 1:nrow(new.cnvs)){
  
  new.cnvs[i, 11] <- paste0(searchTFBS(new.cnvs[i, 5], new.cnvs[i, 8], new.cnvs[i, 9])[,4], collapse=" ")
  
}

searchRegulatoryElements <- function(chr, cnv_from, cnv_to){
  chr <- paste0("chr", chr, collapse="")
  
  temp.gs <- subset(genome.segmentations, genome.segmentations[, 1] == chr)
  
  upstream <- temp.gs[, 8] - cnv_from
  upstream[upstream < 0] <- NA
  
  
  downstream <- cnv_to - temp.gs[, 7]
  upstream[downstream < 0] <- NA
  
  
  temp.gs <- temp.gs[!is.na(upstream), ]
  
  return(temp.gs)
  
}

searchTFBS <- function(chr, cnv_from, cnv_to){
  chr <- paste0("chr", chr, collapse="")
  
  temp.tfbs <- subset(tfbs.161, tfbs.161[, 1] == chr)
  
  upstream <- temp.tfbs[, 3] - cnv_from
  upstream[upstream < 0] <- NA
  
  
  downstream <- cnv_to - temp.tfbs[, 2]
  upstream[downstream < 0] <- NA
  
  
  temp.tfbs <- temp.tfbs[!is.na(upstream), ]
  
  
  return(temp.tfbs)
}
