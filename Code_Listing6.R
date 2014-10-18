# Code Listing 6: Generates annotation of all CNVs found (Table 16)
# Copyright Sebastian Burgstaller-Muehlbacher, 2014. Distributed under GPLv3 or later, without warranty.

# list all genes, affected by a cnv, irrespective of DGV
# list all cnvs present in multiple samples + genes
# intergentic and regulatory regions affected?

library("xlsx")

all.cnvs <- data.frame()

for (i in 1:10) {
  tmp <- getAllCNV(i)
  
  #replace chip name with patient study code
  tmp[, 1] <- patient_names[i]
  
  
  all.cnvs <- rbind(all.cnvs, tmp)
}

all.cnvs <- subset(all.cnvs, all.cnvs[, 4] != "Y")


cnv.gained <- subset(all.cnvs, all.cnvs[, 3] == "Gain")
cnv.lost <- subset(all.cnvs, all.cnvs[, 3] == "Loss")


cnv.gained <- uniqueCNVs(cnv.gained)
cnv.lost <- uniqueCNVs(cnv.lost)

genes.gained <- getGenesInCNVList(cnv.gained)
genes.lost <- getGenesInCNVList(cnv.lost)


write.xlsx(genes.gained, "CNV-annotation.xlsx", sheetName="genes.gained", row.names=F)
write.xlsx(genes.lost, "CNV-annotation.xlsx", sheetName="genes.lost", row.names=F, append=T)


getGenesInCNVList <- function(cnv.list) {
  col.index.genes <- ncol(cnv.list) + 1
  

  for(i in 1:nrow(cnv.list)){
    genes <- mapToGene(cnv.list[i, 4], cnv.list[i, 11], cnv.list[i, 12])
    if(length(genes) == 0){
      genes <- getNeighbourGenes(cnv.list[i, 4], cnv.list[i, 11], cnv.list[i, 12])
    }
    cat(genes)
    cnv.list[i, col.index.genes] <- paste0(genes, collapse= " ")
    
  }
  
  colnames(cnv.list)[col.index.genes] <- "Genes"
  return(cnv.list)
}

