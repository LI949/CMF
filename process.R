library(Seurat)
library(limma)
library(dplyr)
library(magrittr)
library(rmarkdown)
library(tidyverse)
library(patchwork)
Ze= read.table("data/GSE60361_Zeisel.txt.gz",header = T,sep = "\t",quote = "",fill = T)
rt2=as.matrix(Ze)
rownames(rt2)=rt2[,1]
Zeiselexp=rt2[,2:ncol(rt2)]
Zeisel=avereps(Zeiselexp)

sc <- CreateSeuratObject(counts = Zeiselexp,project = "seurat", 
                           min.cells = 3, min.features = 1, names.delim = "_",)
#min.features规定每个细胞需要被检测到的最少基因数,将会筛掉低质量的细胞，
sc
QCZeisel=as.matrix(sc[["RNA"]]@counts)
write.csv(QCZeisel,"D:\\Zeisel.csv")
rowsum = apply(QCZeisel,1,sum)
NolmaZel=sweep(QCZeisel,1,rowsum,"/")*median(rowsum)
log10Zel=log10(NolmaZel+1)
write.csv(NolmaZel,"D:\\NolmaZel.csv")
