## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE---------------
#rm(list = ls())
rm(list = ls())
gc()
library(Seurat)
#lpsdata<-read.table("data/uso.txt",header=T,row.names=1,sep="\t",check.names=F)
lpsdata<-read.csv("data/iPSCdata.csv")
rt2=as.matrix(lpsdata)
rownames(rt2)=rt2[,1]
exp=rt2[,2:ncol(rt2)]
lpsdata=exp
colnames(lpsdata)=c(1:length(colnames(lpsdata)))
#lpsdata<-read.csv("data/Usoskindata.csv")
seurat <- CreateSeuratObject(counts = lpsdata,project = "seurat", min.cells = 3,min.features = 1)
seurat
data1 <- as.matrix(seurat@assays$RNA@data) 
#H=apply(data1,2,sum)
#data2= sweep(data1,2,H,"/")*median(H)
#data3=log10(data2+1)

#class.label<- read.table("data/us_label.txt", header=T,sep="\t",check.names=F) 
class.label<- read.csv("data/iPSC_label.csv")
class.label<-as.matrix(class.label)
data2<-read.csv("data/iPSC_CMF.csv",header = FALSE)
lpsdata1=as.matrix(data2)

#lpsdata=log10(lpsdata+1)
yan1<-lpsdata1
label<-class.label[,2]
## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE---------------
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png')

## ---- message=FALSE, warning=FALSE-----------------------------------------
library(SingleCellExperiment)
library(SC3)
library(scater)
sce <- SingleCellExperiment(
    assays = list(
    counts = as.matrix(yan1),
    #logcounts = log10(as.matrix(yan1) + 1)
    logcounts = yan1
  ), 
   colData = label
    #as.vector(as.matrix(label[,1]))
)

# 定义 feature names 在 feature_symbol 列
rowData(sce)$feature_symbol <- rownames(data1)
# 去除feature中重复的基因名
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

# define spike-ins
#isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)

## --------------------------------------------------------------------------
K=5  #uso 4 loh 8 iPSC 5 Zeisel 11
sce <- sc3(sce,gene_filter = FALSE, ks =K, biology = TRUE)

#sc3_n_clusters = "sc3_4_clusters"
#plotPCA(sce, colour_by =sc3_n_clusters)
preditlabel=sce@colData@listData[["sc3_5_clusters"]]
#write.csv(preditlabel,file = "sc3_results_label_loh.csv",row.names = T)
adjustedRandIndex(preditlabel,class.label[,2])

T=table(preditlabel)
print(T)
Ttrue=table(label)
print(Ttrue)
N=length(colnames(data1))
test.entropy <- function(d){
  res <- 0
  for(i in 1:length(d))
  {
    if(d[i]!=0)
      res <- res + d[i]*log(d[i])
  }
  return (-res)
}
#loh
#d1=c(87 ,75 ,64 ,60, 45, 55  ,22 ,21  )  #随数据集更改
#d2=c(87, 61, 55, 28, 45, 64 ,75, 14)
#uso
#d1=c(234, 165 ,162 ,61)  #随数据集更改
#d2=c(233, 81 ,169 ,139)
#iPSC
d1=c(83 ,88, 80, 49,15 )  #随数据集更改
d2=c(83, 77, 60, 41, 54 )
HU=test.entropy(d=d1/N)
HV=test.entropy(d=d2/N)

U= matrix(data=c(1:K*K),nrow=K,ncol=K,byrow=FALSE,dimnames=NULL)
for (i in 1:K)
  for (j in 1:K)
  {
    U[i,j]=length(intersect(which(preditlabel==i),which(label==j)))
  }
#U列联表,U[i,j]=|Ui交Vj|
mu_infor <- function(u,v,U){
  res <- 0
  for(i in 1:length(u))
    for (j in 1:length(v))
    {
      if(U[i,j]!=0)
        res <- res + U[i,j]/N*log(N*U[i,j]/u[i]/v[j])
    }
  return (res)
}
muinfor=mu_infor(d1,d2,U)
NMI=2*muinfor/(HU+HV)
print(NMI)
## --------------------------------------------------------------------------

#sc3_export_results_xls(sce)

