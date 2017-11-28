source('sbr.R')
library(mvtnorm)
library(MCMCpack)
library(doParallel)

################# Toy example with 3 data sources #################

### GENERATION OF DATA ###
## sample size and number of predictors
n<-50 
p1<-10
p2<-100
p3<-500

## generation of covariance matrices and feature matrices
S1<-riwish(p1,diag(p1))
S2<-riwish(p2,diag(p2))
S3<-riwish(p3,diag(p3))
X1<-matrix(rmvnorm(n*p1,rep(0,p1),S1),n,p1)
X2<-matrix(rmvnorm(n*p2,rep(0,p2),S2),n,p2)
X3<-matrix(rmvnorm(n*p3,rep(0,p3),S3),n,p3)

## sparsity and generation of betas
s2<-p2*0.3
s3<-p3*0.01
non.zero2<-sample(1:p2,s2)
non.zero3<-sample(1:p3,s3)

b1<-rnorm(10,0,2.5)
b2<-rep(0,p2)
b2[non.zero2]<-rnorm(s2)
b3<-rep(0,p3)
b3[non.zero3]<-rnorm(s3)

## generation of responce
mu<-X1%*%b1+X2%*%b2+X3%*%b3
y<-rnorm(n,mu,sd=0.5)

### STANDARIZING & GRAM MATRIX COMPUTATION ###
y<-scale(y)
X1<-scale(X1)
X2<-scale(X2)
X3<-scale(X3)

## first a test
T0<-X3%*%t(X3)                         # usual matrix multiplication                          
T1<-gram(X3,block=TRUE,block.size=100) # block matrix multiplication

cores<-detectCores()-1
cores<-makeCluster(cores)
registerDoParallel(cores)
T2<-gram.parallel(X3,cl=cores,block=FALSE)               # matrix multiplication in parallel
T3<-gram.parallel(X3,cl=cores,block=TRUE,block.size=100) # block matrix multiplication in parallel
stopCluster(cores)

all.equal(T0,T1)
all.equal(T0,T2)
all.equal(T0,T3)
rm(T0,T1,T2,T3)

## calculation of gram matrices
G1<-gram(X1); G2<-gram(X2); G3<-gram(X3)

## make lists
G<-list(G1,G2,G3)
X<-list(X1,X2,X3)

### RUN SBR/SSBR ###

## SBR with the ML lambda-estimator

model1<-sbr(y=y,X=X,G=G,estimator='ML') 

## SSBR with the ML lambda-estimator using block-matrix 
## computations for the variances of X3 (since p3=500)

model2<-sbr(y=y,X=X,G=G,estimator='ML',sparsify=TRUE,p.threshold=100,cov.blocks=100) 

## SSBR with the ML lambda-estimator using block-matrix 
## computations for the variances of X3 (since p3=500)
## and control equal to log(n) for the effect of sample size

model3<-sbr(y=y,X=X,G=G,estimator='ML',sparsify=TRUE,p.threshold=100,cov.blocks=100,sparse.control=log(n))

## parallel computing for the configuration of model3

cores<-detectCores()-1
cores<-makeCluster(cores)
registerDoParallel(cores)
model4<-sbr(y=y,X=X,G=G,estimator='ML',parallel=TRUE,cl=cores,sparsify=TRUE,p.threshold=100,cov.blocks=50,sparse.control=log(n))
stopCluster(cores)

### EXTRACTING OUTPUT FROM A MODEL ###

model3@BetaSBR  # SBR coefficients
model3@BetaSSBR # SSBR coefficients
model3@Lambda   # vector of lambdas  
model3@Sigma2   # error variance  
model3@Duration # runtime