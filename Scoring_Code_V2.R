###ASSUME YOU ALREADY CLUSTERED EXPERMINETAL DATA AND YOU HAVE 4 CLUSTERS,
#THE CLUSTER SIZES ARE IN THE VECTOR clSizeExp
#EXPERIMENTAL DATA IS IN dataM, 
#dataM HAS CELLS AS ROWS AND GENES/TFs AS COLUMNS

dataMt<-read.csv(file ="~/Desktop/Capstone/processesedData (2).csv")
lrgdatasetnormcolgene <- readRDS("~/Desktop/Capstone/lrgdatasetnormcolgene")
#fix row name problem
rNames <- dataMt[,1]
rownames(dataMt) <- rNames
dataMt <- dataMt[,-1]
dataM <- t(dataMt)
dataM <- dataM[,colnames(dataM) %in% colnames(lrgdatasetnormcolgene)]
clSizeExp<-c(200,200,200,200,200,200)

#Find the centers of the 6 clusters

cumClSize<-cumsum(clSizeExp)
m1<-apply(dataM[1:cumClSize[1],],2,median)
m2<-apply(dataM[(cumClSize[1]+1):cumClSize[2],],2,median)
m3<-apply(dataM[(cumClSize[2]+1):cumClSize[3],],2,median)
m4<-apply(dataM[(cumClSize[3]+1):cumClSize[4],],2,median)
m5<-apply(dataM[(cumClSize[4]+1):cumClSize[5],],2,median)
m6<-apply(dataM[(cumClSize[5]+1):cumClSize[6],],2,median)

#Find radius of cluster 1
dcall<-dist(rbind(m1,dataM[1:cumClSize[1],]),method="euclidian")
dc<-dcall[1:cumClSize[1]]
sdc<-sort(dc)
threshDist<-0.95*clSizeExp[1] #you can replace 0.95 by 0.97 if you want a bigger radius
r1<-sdc[threshDist]

#Find radius of cluster 2
dcall<-dist(rbind(m2,dataM[(cumClSize[1]+1):(cumClSize[2]),]),method="euclidian")
dc<-dcall[1:clSizeExp[2]]
sdc<-sort(dc)
threshDist<-0.95*clSizeExp[2]
r2<-sdc[threshDist]

#Find radius of cluster 3
dcall<-dist(rbind(m2,dataM[(cumClSize[2]+1):(cumClSize[3]),]),method="euclidian")
dc<-dcall[1:clSizeExp[3]]
sdc<-sort(dc)
threshDist<-0.95*clSizeExp[3]
r3<-sdc[threshDist]

#Find radius of cluster 4
dcall<-dist(rbind(m4,dataM[(cumClSize[3]+1):(cumClSize[4]),]),method="euclidian")
dc<-dcall[1:clSizeExp[4]]
sdc<-sort(dc)
threshDist<-0.95*clSizeExp[4]
r4<-sdc[threshDist]

#Find radius of cluster 5
dcall<-dist(rbind(m5,dataM[(cumClSize[4]+1):(cumClSize[5]),]),method="euclidian")
dc<-dcall[1:clSizeExp[5]]
sdc<-sort(dc)
threshDist<-0.95*clSizeExp[5]
r5<-sdc[threshDist]

#Find radius of cluster 6
dcall<-dist(rbind(m4,dataM[(cumClSize[5]+1):(cumClSize[6]),]),method="euclidian")
dc<-dcall[1:clSizeExp[6]]
sdc<-sort(dc)
threshDist<-0.95*clSizeExp[6]
r6<-sdc[threshDist]

#PUT THE CLUSTER CENTERS IN A MATRIX, CALL IT cenMedRef
#AND PUT THE RADII IN A VECTOR, CALL IT cutOffM

cenMedRef<-rbind(rbind(m1,m2,m3,m4,m5,m6))
cutOffM<-c(r1,r2,r3,r4,r5,r6)
numbClusters<-dim(cenMedRef)[1] 

#FIND THE SCORE OF EACH SIMULATED MODEL
#ASSUME THE SIMULATED DATA IS IN A MATRIX simulM
#EACH SIMULATED MODEL IS A ROW, GENES/TFs ARE COLUMNS

lrgdatasetnormcolgene <- readRDS("~/Desktop/Capstone/lrgdatasetnormcolgene")

simulM<-lrgdatasetnormcolgene

sceGeneNames<-colnames(dataM)
simGeneNames<-colnames(simulM)
permN<-rep(0,length(sceGeneNames))
for(i in 1:length(permN))
{
  permN[i]<-which(simGeneNames==sceGeneNames[i])[1]
}
shSimM<-simulM[,permN]

M<-dim(shSimM)[1] #this is number of simulated models
#M<-10000

clusterScoreSimM<-matrix(0, nrow=M, ncol=2)

#clusterScoreSimM will contain on 2nd column the score of each simulated model
#in the first column we will have the cluster membership of each sim. model


for (t in 1:M)
{
  #For each model calculate the distances to cluster centers
  
  dmedMod<-dist(rbind(shSimM[t,],cenMedRef), method="euclidean")
  dIntm<-dmedMod[1:numbClusters]
  
  #Divide the distances by the corresponding radii of the clusters
  ratDm<-dIntm/cutOffM
  
  #Find the cluster that minimizes the ratios 
  indClm<-which(ratDm==min(ratDm))[1] #assign the cluster membership
  
  if(ratDm[indClm]<=3)
  {
    clusterScoreSimM[t,1]<-indClm  
  }
  
  #If minimum ratio is bigger than 1 assign it to a cluster "0", that is
  #treat it as noise
  
  if(ratDm[indClm]>3)
  {
    clusterScoreSimM[t,1]<-0
  }
  
  clusterScoreSimM[t,2]<-ratDm[indClm] #assign the score to model t
  
}

scoresOut<-rep(0,7) #this will contain a bunch of scores

indZe<-which(clusterScoreSimM[,1]==0)
scoresOut[1]<-length(indZe)#this is number of "noise" simulated models


clSizeSimul<-rep(0,numbClusters)#this will contain number of simulated models assigned to each cluster
avgCl<-rep(0,numbClusters)#this will contain median score for each cluster
chScore<-FALSE #boolean value that indicates if we should impose a penalty to the score

for(j in 1:numbClusters)
{
  indCAv<-which(clusterScoreSimM[,1]<=j)
  avgCl[j]<-median(clusterScoreSimM[indCAv,2])
  clSizeSimul[j]<-length(indCAv)
  
  #See if need to impose a penalty for not having a good enough distribution
  #of the simulated models
  
  if((clSizeSimul[j]/sum(clSizeSimul))<(0.5*clSizeExp[j]/sum(clSizeExp)))
  {
    chScore<-TRUE
  }
  
}

scoresOut[2]<-(M/numbClusters)*sum(avgCl)
scoresOut[3]<-sum(clusterScoreSimM[indZe,2])
#scoresOut[4]<-sum(clusterScoreSimM[indZe,3])
#scoresOut[5]<-sum(clusterScoreSimM[indZe,4])
#scoresOut[6]<-sum(clusterScoreSimM[indZe,5])
scoresOut[7]<-scoresOut[2]+scoresOut[3]
#+scoresOut[4]+scoresOut[5]+scoresOut[6]

if(chScore==TRUE)
{
 scoresOut[7]<-scoresOut[7]+50000
}

####scoreOut[7]  IS THE FINAL SCORE OF THE SCENIC NETWORK

####IF YOU WANT TO WRITE A FUNCTION THAT AUTOMATICALLY CALCULATES THE SCORE OF A 
####SCENIC-GENERATED NETWORK, YOU WILL NEED  A TOPOLOGY OF A NETWORK PRODUCED
###BY SCENIC, CALL IT topol_scenic
###THEN YOU CAN WRITE A FUNCTION LIKE THIS

score_scenic<-function(clSizeExp, cenMedRef,cutOffM, topol_scenic)
{
 M<-10000
 racSCE<-sracipeSimulate(circuit = topol_scenic, numModels = M, plots = FALSE, 
                        integrateStepSize = 0.1, simulationTime = 100)

 racNSCE<-sracipeNormalize(racSCE)
 dataRowGeneSCE<-assay(racNSCE,1)
 shSimM<-t(as.matrix(dataRowGeneSCE))
 
 ###COPY CODE FROM ABOVE
 clusterScoreSimM<-matrix(0, nrow=M, ncol=5)

 #clusterScoreSimM will contain on 2nd column the score of each simulated model
 #in the first column we will have the cluster membership of each sim. model


 for (t in 1:M)
 {
  ...
 }

 ...

 if(chScore==TRUE)
 {
   scoresOut[7]<-scoresOut[7]+50000
 }

return(scoresOut[7])
}#end function

