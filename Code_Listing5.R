# Code Listing 5: Join expression levels from melanocytic cell line and add expression levels to annotation.
# Copyright Sebastian Burgstaller-Muehlbacher, 2014. Distributed under GPLv3 or later, without warranty.
# evaluate gene expression of genes affected by CNVs in "Melano" melanomcyte line from ENCODE

#load 4 gene expression replicates from the "Melano" human melanocyte cell line and generate a mean expression of those 4 datasets
#Gene expression data was generated on the Affymetrix Exon 1.0 ST arrays

#Gene expression omnibus accession numbers: GSM993582, GSM993583, GSM993584, GSM993585

melano.expr1 <- read.table("./expression/GSM993582_Melano_B1_expression.bed.gz", stringsAsFactors = F)
colnames(melano.expr1) <- c("chr", "start", "end", "gene", "score", "strand", "signalValue", "exonCount", "constituitiveExons")
melano.expr2 <- read.table("./expression/GSM993583_Melano_B2_expression.bed.gz", stringsAsFactors = F)
melano.expr3 <- read.table("./expression/GSM993584_Melano_B3_expression.bed.gz", stringsAsFactors = F)
melano.expr4 <- read.table("./expression/GSM993585_Melano_B4_expression.bed.gz", stringsAsFactors = F)


expression.mean <- c()
expression.stdev <- c()


for(i in 1:nrow(melano.expr1)) {
  expression.mean <- c(expression.mean, mean(c(melano.expr1[i, 7], melano.expr2[i, 7], melano.expr3[i, 7], melano.expr4[i, 7])))
  expression.stdev <- c(expression.stdev, sd(c(melano.expr1[i, 7], melano.expr2[i, 7], melano.expr3[i, 7], melano.expr4[i, 7])))
}

melano.expr1[, 7] <- expression.mean
melano.expr1[, 10] <- expression.stdev
colnames(melano.expr1)[10] <- "expr.stdev"


for(i in 1:nrow(new.cnvs)){
  
  new.cnvs[i, 19] <- getExpressionLevel(new.cnvs[i, 16])
  
}

colnames(new.cnvs)[19] <- "gene.expr"

write.xlsx(new.cnvs, "CNV-annotation_encode+expression.xlsx")


getExpressionLevel <- function(cnv.genes) {
  
  gene.list <- c()
  if(grepl("\\|", cnv.genes)) {
    gene.list <- unlist(strsplit(x=cnv.genes, " \\| "))
  } else {
    gene.list <- unlist(strsplit(x=cnv.genes, " "))
  }
    
  gene.subset <- melano.expr1[melano.expr1[, 4] %in% gene.list, ]
  
  #round expression levels to 2 digits
  gene.subset[, 7] <- round(gene.subset[, 7], digits = 2)
  
  #return a string with the expression levels only
  return(paste0(gene.subset[, 7], collapse=" "))
  
}


###########################
#load encode annotated multi-case cnvs and add gene expression level
multi.case.new <- read.xlsx(file="new_cnvs_encode.xlsx", stringsAsFactors=F, sheetIndex=1)

for(i in 1:nrow(multi.case.new)){
  
  multi.case.new[i, 13] <- getExpressionLevel(multi.case.new[i, 3])
  
}

colnames(multi.case.new)[13] <- "gene.expr"

write.xlsx(multi.case.new, file="new_cnvs_encode+expression.xlsx", row.names=F)

