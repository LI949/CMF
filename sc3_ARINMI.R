
lpsdata<-read.csv("data/uso_CMF.csv",header = FALSE)
#colnames(lpsdata)=c(1:622)
data2= as.matrix(lpsdata)

class.label<- read.table("data/us_label.txt", header=T,sep="\t",check.names=F) 
#class.label<- read.csv("data/iPSC_label.csv")
class.label<-as.matrix(class.label)

yan1<-as.matrix(data2)
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
    logcounts = log10(as.matrix(yan1) + 1)
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
K=4
sce <- sc3(sce,gene_filter = FALSE, ks =K, biology = TRUE)

#sc3_n_clusters = "sc3_4_clusters"
#plotPCA(sce, colour_by =sc3_n_clusters)
preditlabel=sce@colData@listData[["sc3_4_clusters"]]
adjustedRandIndex(preditlabel,class.label[,1])

label=class.label[,2]
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

d1=c(234, 161 ,167 ,60)  #随数据集更改
d2=c(233 ,81 ,169,139)
d1=c(964,   2, 690, 506,  37, 317, 188 ,246,  55)  #随数据集更改
d2=c(290, 390, 948, 820,  98, 175, 198,  26,  60)
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
