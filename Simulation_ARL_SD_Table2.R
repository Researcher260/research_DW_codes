library(pracma)
library(matrixStats)
library(DiscreteWeibull)
library(writexl)
tic()
clear()
# Definition of parameters n, q, beta
target=200 #ARL0
#In-control 
q = 0.5
beta = 0.5
set.seed(123)
#out-of-control - Table 2
Q1 = c(0.5,0.6,0.8,0.5,0.5,0.4,0.55)
Beta1 = c(0.5,0.5,0.5,0.4,0.3,0.4,0.4)

N=c(1,2,3,5,7,10,30,50,100,300) #Sample size used
n=300 #Sample size.
Ref=which(N == n) #Position of n in N 
FREQ=c(rep(10000000,5),rep(1000000,2),600000,300000,100000)
freq=FREQ[Ref]
Result=matrix(0,7,2)
for (j in 1:7){
  cat('Out-of-Control. Completed in 7: ',j,"\n")
  q1=Q1[j]
  beta1=Beta1[j]

Simulate=50000 #Monte Carlo - mixed geometric distribution
runs<-1000 #number of simulations for Monte Carlo
Saida=matrix(0,runs,1)
for (i in 1:runs){
cat('Runs in 1000: ',i,"\n")

SC0<-(matrix(rdweibull(n*freq,q,beta,TRUE),n,freq)) 
SCmean0<-colMeans(SC0) 

#Limits obtained by experimentation to reach target ARL0
LI=c(57,38.5,30.5,22.5,18.5,15.5,9.5,7.5,5.5,4.5)
LS=c(59,40.5,32.5,24.5,20.5,17.5,11.5,9.5,7.5,6.5)

li=LI[Ref] 
ls=LS[Ref]
D=seq(li,ls,0.01)
DD=size(D)
DD=DD[2]

ARLM=matrix(0,DD,3)
for (jj in 1:DD){
 LL=D[jj]
 Prob1=mean(SCmean0<=LL)
 ARL0=1/(1-Prob1)
 ARLM[jj,1]=ARL0
 ARLM[jj,2]=LL
}

ARLM[,3]=abs(ARLM[,1]-target)
ARLM=sortrows(ARLM,3)
ARL=1
cont=1
while (ARL<target){
ARL=ARLM[cont,1]
LSC=ARLM[cont,2]
cont=cont+1
}
#Distribution out-of-control
SC1<-(matrix(rdweibull(n*freq,q1,beta1,TRUE),n,freq)) 
SCmean1<-colMeans(SC1) 
Prob1=mean(SCmean1<=LSC)
ARL1=1/(1-Prob1)

Saida[i,1]=ARL1

}
SaidaA=matrix(0,Simulate,runs)
for (ii in 1:runs){ 
  pp=1/Saida[ii,1]
  R=rgeom(Simulate,pp)+1
  SaidaA[,ii]=R
  cat('Simulate mixed geometric in 1000: ',ii,"\n")
}
Result[j,1]=mean(SaidaA)
Vari=sum((SaidaA-mean(SaidaA))^2)/(runs*Simulate-1)

Result[j,2]=Vari^0.5

}
Result=data.frame(Result)
writexl::write_xlsx(Result,"Result.xlsx")
toc()

