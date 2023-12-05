library(tidyverse)
library(DiscreteWeibull)
library(expm)
library(pracma)
library(writexl)
tic()
clear()
# Definition of parameters n, q, beta
n <- 1 #Sample size. 
#In-control average
q <- 0.5
beta <- 0.5
Alvo<-200
cumulative_probability=0.999
#n=5 and cumulative_probability=0.999 implies an elapsed time of 2149 seconds

#out-of-control - Table 2

Q1 = c(0.5,0.6,0.8,0.5,0.5,0.4,0.55)
Beta1 = c(0.5,0.5,0.5,0.4,0.3,0.4,0.4)

# Computing k and x
k <- qdweibull(cumulative_probability, q, beta, zero = TRUE)
x <- 0:(k*n - 1) 
# Computing P(X = x)
# x = {0, 1, ..., k*n - 1}
prob <- pmap_dbl(.l= list(x=x, q=list(q), beta=list(beta),
                          zero = list(TRUE)), .f = ddweibull) 
# Creating of matrix
creating_matrix <- function(prob){
  prob_last <- 1-sum(prob)
  k_n <- length(prob) + 1  
  prob_matrix <- matrix(0, nrow = k_n,ncol = k_n)
  prob_matrix[1,] <- c(prob,prob_last)  
  for(i in 2:k_n){
    prob_matrix[i,] <- c(0, prob_matrix[i-1,1:(k_n-2)],
                         sum(prob_matrix[i-1,-(1:(k_n-2))]))
  }
  return(prob_matrix)
}
(matrix <- creating_matrix(prob))    
# Creating of transition matrix Q and 
# Determination of the sum distribution
# of X by Markov Chain approach
Q_exp_n <- matrix %^% n
vetor_b <- c(1,rep(0, length(prob)))
v <- t(vetor_b) %*% Q_exp_n
acum_matriz <- cumsum(t(vetor_b) %*% Q_exp_n)
p_x <- seq(0, length(prob), 1)
MC_tab_probs <- data.frame(n,q,beta,p_x,acum_matriz)
MC_tab_probs
A=abs(MC_tab_probs[,5]-1/Alvo)
MC_tab_probs=cbind(MC_tab_probs,A)

MC_tab_probs=data.matrix(MC_tab_probs)
MC_tab_probs=sortrows(MC_tab_probs,6)
ARL0=1
i=1
while (ARL0< Alvo){
ARL0=1/(1-MC_tab_probs[i,5])
UCL=MC_tab_probs[i,4]/n
i=i+1
}

ResultCM=matrix(0,7,4)
for (j in 1:7){
  q1=Q1[j]
  beta1=Beta1[j]
q <- q1
beta <- beta1
# Computing k and x
k <- qdweibull(cumulative_probability, q, beta, zero = TRUE)
x <- 0:(k*n - 1) 
# Computing P(X = x)
# x = {0, 1, ..., k*n - 1}
prob <- pmap_dbl(.l= list(x=x, q=list(q), beta=list(beta),
                          zero = list(TRUE)), .f = ddweibull) 
# Creating of matrix
creating_matrix <- function(prob){
  prob_last <- 1-sum(prob)
  k_n <- length(prob) + 1  
  prob_matrix <- matrix(0, nrow = k_n,ncol = k_n)
  prob_matrix[1,] <- c(prob,prob_last)  
  for(i in 2:k_n){
    prob_matrix[i,] <- c(0, prob_matrix[i-1,1:(k_n-2)],
                         sum(prob_matrix[i-1,-(1:(k_n-2))]))
  }
  return(prob_matrix)
}
(matrix <- creating_matrix(prob))    
# Creating of transition matrix Q and 
# Determination of the sum distribution
# of X by Markov Chain approach
Q_exp_n <- matrix %^% n
vetor_b <- c(1,rep(0, length(prob)))
v <- t(vetor_b) %*% Q_exp_n
acum_matriz <- cumsum(t(vetor_b) %*% Q_exp_n)
p_x <- seq(0, length(prob), 1)
MC_tab_probs <- data.frame(n,q,beta,p_x,acum_matriz)
MC_tab_probs
ARL1=1/(1-MC_tab_probs[(UCL*n+1),5])
ResultCM[j,3]=ARL1
pp=1/ARL1
ResultCM[j,4]=((1-pp)/(pp^2))^0.5
ResultCM[j,1]=q1
ResultCM[j,2]=beta1
}

colnames(ResultCM)<-c("q1","beta1", "ARL", "SDARL")
ResultCM=data.frame(ResultCM)
writexl::write_xlsx(ResultCM,"ResultCMa.xlsx")

toc()