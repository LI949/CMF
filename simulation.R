## ---- include = FALSE------------------------------------------------------
rm(list = ls())
gc()
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup-----------------------------------------------------------------
library("splatter")
library("scater")
library("ggplot2")


# Five groups
## ----nGenes----------------------------------------------------------------
# Set the number of genes to 1000
params = newSplatParams()

params = setParams(params, list(batchCells = 500,
                                nGenes =1000,
                                group.prob = c(0.3, 0.3,0.4),
                                de.prob = c(0.05, 0.08, 0.01),#基因进行差异表达的概率
                                de.facLoc = 0.5,
                                de.facScale = 0.8)
)

# Set up the vector of dropout.mid
#dropout_mid = c(4, 5, 5.5)  ,c(0,0,0)

# determine if it is a good parameter


# Generate the simulation data using Splatter package
sim = splatSimulateGroups(params,
                          dropout.shape =c(-0.4,-0.4,-0.4),
                          dropout.mid = c(4, 5, 5.5),
                          dropout.type = "group", #确定要模拟的dropout效果的类型
                          )
#dropout.mid参数控制概率等于0.5的点
#dropout.shape控制概率如何随表达式的增加而变化

sim <- logNormCounts(sim)
#sim <- runPCA(sim)
#plotPCA(sim, colour_by = "Group")
sim <- runTSNE(sim,ncomponents = 2)
plotTSNE(sim, colour_by = "Group")

simtrue <- as.matrix(sim@assays@data@listData[["TrueCounts"]])
simtrue2 <- log10(simtrue+1)

X <- as.matrix(assays(sim)$count)
X.log <- log10(X+ 1)
simlabel<-sim$Group
write.csv(X.log,file = "sim-0.4.csv")
write.csv(simtrue2,file = "sim-0.4true.csv")
write.csv(simlabel,file = "sim-0.4_label.csv")

tsne_re<-Rtsne(t(simtrue2),perplexity=10,check_duplicates = FALSE)
tSNEdata<-data.frame(tsne_re$Y,simlabel)
names(tSNEdata) <- c("tSNE_1", "tSNE_2","cluster")

tSNEdata$cluster=as.character(tSNEdata$cluster)
dev.off() 
ggplot(tSNEdata,mapping = aes(x=tSNE_1,
                              y=tSNE_2,
                              col=cluster))+
  geom_point()


tsne_re<-Rtsne(t(simim0.4),perplexity=10,check_duplicates = FALSE)
tSNEdata<-data.frame(tsne_re$Y,simlabel)
names(tSNEdata) <- c("tSNE_1", "tSNE_2","cluster")

tSNEdata$cluster=as.character(tSNEdata$cluster)
dev.off() 
ggplot(tSNEdata,mapping = aes(x=tSNE_1,
                              y=tSNE_2,
                              col=cluster))+
  geom_point()

#PCC SSE
simim0.4 <- read.csv("sim-0.4_CMF.csv")
simim0.4 <- simim0.4[,2:501]
simim0.4 <- as.matrix(simim0.4)
X.log = read.csv("sim-0.4.csv")
X.log <- X.log[,2:501]
X.log <- as.matrix(X.log)
#EE1=t(matrix(colSums(simtrue2)/1000,nrow = 500,ncol = 1000))
#EE2=t(matrix(colSums(simim0.4)/1000,nrow = 500,ncol = 1000))
A=c(simtrue2)-mean(c(simtrue2))
B=c(simim0.4)-mean(c(simim0.4))
C=c(X.log)-mean(c(X.log))
AB= sum(B*c(simtrue2))

#s=y1%*%AB
s=AB
PCC=s/norm(A,c("2"))/norm(B,c("2"))#0.9351

AC=sum(C*c(simtrue2))
PCCraw=s/norm(C,c("2"))/norm(B,c("2"))#0.8239

Pearson1 <- round(cor(simim0.4,method = c("pearson")),2)
mean(Pearson1)
Pearson2 <- round(cor(X.log,method = c("pearson")),2)
mean(Pearson2)
rowE=rowSums((simtrue2-X.log)^2)  #行和
y1=matrix(1:1,nrow = 1,ncol = 1000)#新建一个一行1000列的单位矩阵
SSE1 <- (y1%*%rowE)^(1/2) #386

rowEim=rowSums((simtrue2-simim0.4)^2)  #行和
SSEim <- (y1%*%rowEim)^(1/2) #194
