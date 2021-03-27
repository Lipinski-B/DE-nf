#!/usr/bin/Rscript
library(DESeq2)
library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(ggbeeswarm)
library(genefilter)
library(biomaRt)
library(stringr)
library(ggplot2)
library(NMF)
library(tidyverse)
install.packages(as.character(commandArgs(TRUE)[3]), repos = NULL, type="source")
library(grid)
install.packages(as.character(commandArgs(TRUE)[4]), repos = NULL, type="source")
library(gridExtra)

remplace <- function(name){name <- gsub("-", ".", name)} 
split <- function(name){name <- str_split(name,"_")[[1]][1]}

#matrix <- "/home/boris/Bureau/DE/DE-nf/output/finale.txt" 
matrix = as.character(commandArgs(TRUE)[1])
#metadata <- "/home/boris/Bureau/DE/DE-nf/data/Metadata.xls"
metadata = as.character(commandArgs(TRUE)[2])

###############################################################################################################################################
## -- DESeq2 annalyse -- ##
###########################
DEA_analysis <- function(matrix, metadata){
  #load
  countData <- read.table(matrix, row.names = 1, header = T, sep = "\t")
  countData <- head(countData,-5)
  #countData <- countData[1:61]
  
  #cutoff
  cutoff <- as.integer(cpm(10, mean(colSums(countData))))
  keep <- rowSums(cpm(countData)>cutoff) >= 2
  countData <- countData[keep,]
  
  #rename
  colnames(countData) <- apply(as.matrix(colnames(countData)), 1, split)
  countData <- countData[,order(colnames(countData))]
  
  #process metadata
  colData <- readxl::read_xls(path = metadata)
  colData[1] <- apply(colData[1], 1, remplace)
  colData <- colData[order(colData[,1]),]
  
  #intersect with data already processed
  colData <- colData[which(colData$`Name of sample` %in% intersect(colData$`Name of sample`, colnames(countData))),]
  #colData <- colData[which(colData$`Sex`=='Male'),]                            #Pour ne garder QUE les males
  #colData <- colData[which(colData$`Sex`=='Female'),]                          #Pour ne garder QUE les females
  
  #colData <- colData[which(colData$`Modality`!="Gnotobiotic Microbacterium oxydans"),]                        #Pour supprimer les Gnotobiotic Microbacterium oxydans
  #colData <- colData[which(colData$`Modality`!="Gnotobiotic Carnobacterium maltaromatycum"),]                 #Pour supprimer les Gnotobiotic Carnobacterium maltaromatycum
  #colData <- colData[which(colData$`Modality`!="Gnotobiotic Oerskovia"),]                                     #Pour supprimer les Gnotobiotic Oerskovia
  #colData <- colData[which(colData$`Modality`!="Conventional"),]                                              #Pour supprimer les Conventional
  #colData <- colData[which(colData$`Modality`!="Axenic"),]                                                    #Pour supprimer les Axenic
  
  countData <- countData[,which(colnames(countData) %in% intersect(colData$`Name of sample`, colnames(countData)))]
  
  #DESeq2 analysis
  dataset <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~Modality)
  dataset <- estimateSizeFactors(dataset) #+ dataset@colData$sizeFactor
  
  dataset <- DESeq(dataset)

  return(dataset)
}
dataset <- DEA_analysis(matrix, metadata)
rld <- rlog(dataset, blind = T)


# Top Genes Results
DEA_result <- function(dataset){
  result <- results(dataset, pAdjustMethod = "BH")
  result <- result[complete.cases(result),]
  result <- result[order(result$padj),]
  
  #Tri des gènes les plus différentiellement exprimé
  topResults <- as.matrix(subset(result , padj < 0.05))
  
  return(topResults)
}
result <- DEA_result(dataset)
topGenes <- assay(rld)[rownames(result),]



###############################################################################################################################################
## -- Figure time -- ##
#######################
pdf("/mnt/DE/result/R/Result.pdf", height = 12, width = 12)

# Résultat
grid.text("Résultats des top gènes les plus différentiellement exprimés :", x = 0.05, hjust = 0.05, vjust = -20, gp = gpar(fontsize = 25))
grid.table(result)

# Dendrogramme
rlog.norm.counts <- assay(rld)
distance.m_rlog  <- as.dist(1 - cor(rlog.norm.counts , method = "pearson" ))
plot(hclust(distance.m_rlog), labels = colnames(rlog.norm.counts),main = "rlog  transformed  read  counts\ndistance: Pearson  correlation")

# PCA
plotPCA(rld, intgroup = "Modality") + theme_bw() + ggtitle("Rlog transformed counts")
plotPCA(rld, intgroup = "Sex")+ theme_bw() + ggtitle("Rlog transformed counts")

# aheatmap
aheatmap(topGenes, Colv = TRUE, hclustfun = "average", scale = "row", main="By modality") 
aheatmap(topGenes, Colv = TRUE, hclustfun = "average", scale = "row", main="By sex") 
#top40Genes <- head(assay(rld)[rownames(topGenes),], 40)
#aheatmap(top40Genes, Colv = TRUE, hclustfun = "average", scale = "row")

# Plot count
topGene <- rownames(result)[which.min(result[,"padj"])]
els <- plotCounts(dataset, gene = topGene, intgroup="Modality", returnData=TRUE)
els$sex = dataset@colData$Sex
ggplot(els,aes(x = Modality, y = count)) + ggtitle(topGene) + geom_point(aes(color = sex))

# PlotMA
plotMA(topGenes, main = "MA plot: group1 vs group2 (padj < 0.05)", ylim = c(-4,4))
abline(h = c(-1,1), col = "blue", lty = 2)
mtext(c(paste("-2 fold"), paste("+ 2 fold")), side = 4, at = c(-1, 1), cex = 0.8, line = 0.5, col = "blue")

dev.off()