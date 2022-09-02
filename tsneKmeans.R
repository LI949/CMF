rm(list = ls())
library(devtools)
library(mclust)
library(Rtsne)
library(ggplot2)
# before imputation

#lpsdata<-read.table("data/uso.txt",header=T,row.names=1,sep="\t",check.names=F)
lpsdata<-read.csv("data/iPSC_CMF.csv", header=FALSE)  #uso_CMF0.01
#class.label<- read.table("data/us_label.txt", header=T,sep="\t",check.names=F) 
class.label<-read.csv("data/iPSC_label.csv")
class.label<-as.matrix(class.label)

K=5
set.seed(2017)#207 209
tSNE.resul=Rtsne(t(as.matrix(lpsdata)))
kmeans.resul=kmeans(tSNE.resul$Y, centers = K)
adjustedRandIndex(kmeans.resul$cluster, class.label[,1])

pca.rmd <- prcomp(t(lpsdata))  #PCA函数
cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)#Kmeans聚类
ARI <- adjustedRandIndex(class.label[,1], cl.rmd$cluster)

preditlabel=kmeans.resul$cluster
label=class.label[,2]
T=table(preditlabel)
print(T)
Ttrue=table(label)
print(Ttrue)
N=length(colnames(lpsdata))
test.entropy <- function(d){
  res <- 0
  for(i in 1:length(d))
  {
    if(d[i]!=0)
      res <- res + d[i]*log(d[i])
  }
  return (-res)
}

d1=c(38, 75, 64, 31, 87, 45, 60, 29 )  #随数据集更改
d2=c(87, 61, 55, 28, 45, 64, 75, 14)#loh
#d2=c(233 ,81 ,169,139)  #uso
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
# after imputation
#adjustedRandIndex(kmeans(Rtsne(t(as.matrix(X.imp)))$Y, centers = 4)$cluster, class.label)

cell = class.label[,1]
cluster = kmeans.resul[["cluster"]]
tSNEdata<-data.frame(tSNE.resul$Y,cell,cluster)
names(tSNEdata) <- c("tSNE_1", "tSNE_2", "cell","cluster")
tSNEdata$cluster=as.character(tSNEdata$cluster)
dev.off() 
ggplot(tSNEdata,mapping = aes(x=tSNE_1,
                              y=tSNE_2,
                              col=cell))+
  geom_point()
