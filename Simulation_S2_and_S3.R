# Load libraries
library(mixdist)
library(qvalue)
library(mixtools)
library(data.table)

# Set the number of bits for integer compression
BITS<- 16

# Initialize matrices for storage of results
FDR_ARRAY <- array(NA,c(100,4,2))
TPR_ARRAY <- array(NA,c(100,4,2))
for (ii in 1:100) {
  # Set sample size
  N <- 100000
  # Sample labels according to declared proportions
  components <- sample(1:2,prob=c(0.8,0.2),size=N,replace=TRUE)
  # Set distribution means
  mus <- c(0,2)
  # Set distribution variances
  sds <- sqrt(c(1,1))
  # Sample AR chains
  CHAIN0 <- mus[1] + c(arima.sim(list(ar=c(-0.5)),n=100000,sd=sqrt(1-0.5^2)))
  CHAIN1 <- mus[2] + c(arima.sim(list(ar=c(-0.5)),n=100000,sd=sqrt(1-0.5^2)))
  # Compute test statistics
  TVALUES <- c()
  TVALUES[which(components==1)] <- CHAIN0[which(components==1)]
  TVALUES[which(components==2)] <- CHAIN1[which(components==2)]  
  # Compute p-values
  PVALUES <- 1-pnorm(TVALUES,0,1)
  # Round data according to integer encoding scheme
  a=data.table(Value=PVALUES)
  a[,merge:=Value]
  b=data.table(Value=seq(0,1,length.out = 2^BITS))
  b[,merge:=Value]
  setkeyv(a,c('merge'))
  setkeyv(b,c('merge'))
  encoded=b[a,roll='nearest']
  # Compute z-scores
  zscores <- qnorm(1-encoded$Value)
  ZSCORE <- zscores
  # Determine histogram binning
  HISTDAT <- mixgroup(ZSCORE,breaks=hist(ZSCORE,breaks='Sturges')$breaks)
  HISTDAT[1,2] <- length(which(ZSCORE<=HISTDAT[1,1]))
  HISTDAT[dim(HISTDAT)[1],2] <- HISTDAT[dim(HISTDAT)[1],2]+sum(ZSCORE==Inf)
  # Fit mixture model
  MMIX <- mix(HISTDAT,mixparam(c(0,2),c(1,1),c(.8,.2)),
              emsteps = 5,iterlim=100,steptol=1e-3)

  # Obtain list of component values and round p-values
  CLASS <- components - 1
  ROUNDED <- PVALUES
  
  # Conduct FDR control via EB model
  OUTPUT <- unlist(c(MMIX$parameters[1,c(1:3)],MMIX$parameters[2,2:3]))
  TAU <- OUTPUT[1]*dnorm(ZSCORE,OUTPUT[2],OUTPUT[3])/
    (OUTPUT[1]*dnorm(ZSCORE,OUTPUT[2],OUTPUT[3])+(1-OUTPUT[1])*dnorm(ZSCORE,OUTPUT[4],OUTPUT[5]))
  TAU[which(is.na(TAU)&(ZSCORE==-Inf))] <- 1
  TAU[which(is.na(TAU)&(ZSCORE==Inf))] <- 0
  CUMTAU <- cumsum(sort(TAU))/(1:N)
  REJECT <- rep(0,N)
  REJECT[which(TAU<=CUMTAU[which(CUMTAU>=0.05)[1]-1])] <- 1
  FDR_ARRAY[ii,1,1] <-   sum((REJECT==1)&(CLASS==0))/max(1,sum(REJECT==1))
  TPR_ARRAY[ii,1,1] <-   sum((REJECT==1)&(CLASS==1))/sum(CLASS==1)
  
  REJECT[which(TAU<=CUMTAU[which(CUMTAU>=0.10)[1]-1])] <- 1
  FDR_ARRAY[ii,1,2] <-   sum((REJECT==1)&(CLASS==0))/max(1,sum(REJECT==1))
  TPR_ARRAY[ii,1,2] <-   sum((REJECT==1)&(CLASS==1))/sum(CLASS==1)
  
  # Conduct FDR control via BH
  REJECT <- rep(0,N)
  REJECT[which(p.adjust(ROUNDED,method='BH')<=0.05)] <- 1
  FDR_ARRAY[ii,2,1] <- sum((REJECT==1)&(CLASS==0))/max(1,sum(REJECT==1))
  TPR_ARRAY[ii,2,1] <- sum((REJECT==1)&(CLASS==1))/sum(CLASS==1)
  
  REJECT <- rep(0,N)
  REJECT[which(p.adjust(ROUNDED,method='BH')<=0.10)] <- 1
  FDR_ARRAY[ii,2,2] <- sum((REJECT==1)&(CLASS==0))/max(1,sum(REJECT==1))
  TPR_ARRAY[ii,2,2] <- sum((REJECT==1)&(CLASS==1))/sum(CLASS==1)
  
  # Conduct FDR control via BY
  REJECT <- rep(0,N)
  REJECT[which(p.adjust(ROUNDED,method='BY')<=0.05)] <- 1
  FDR_ARRAY[ii,3,1] <- sum((REJECT==1)&(CLASS==0))/max(1,sum(REJECT==1))
  TPR_ARRAY[ii,3,1] <- sum((REJECT==1)&(CLASS==1))/sum(CLASS==1)
  
  REJECT <- rep(0,N)
  REJECT[which(p.adjust(ROUNDED,method='BY')<=0.10)] <- 1
  FDR_ARRAY[ii,3,2] <- sum((REJECT==1)&(CLASS==0))/max(1,sum(REJECT==1))
  TPR_ARRAY[ii,3,2] <- sum((REJECT==1)&(CLASS==1))/sum(CLASS==1)

  # Conduct FDR control via q-values
  REJECT <- rep(0,N)
  REJECT <- qvalue(ROUNDED,0.05)$significant
  FDR_ARRAY[ii,4,1] <- sum((REJECT==1)&(CLASS==0))/max(1,sum(REJECT==1))
  TPR_ARRAY[ii,4,1] <- sum((REJECT==1)&(CLASS==1))/sum(CLASS==1)

  REJECT <- rep(0,N)
  REJECT <- qvalue(ROUNDED,0.10)$significant
  FDR_ARRAY[ii,4,2] <-sum((REJECT==1)&(CLASS==0))/max(1,sum(REJECT==1))
  TPR_ARRAY[ii,4,2] <- sum((REJECT==1)&(CLASS==1))/sum(CLASS==1)
  
  # Print diagnostics
  print(ii)
  print(FDR_ARRAY[ii,,])
  print(TPR_ARRAY[ii,,])
}

# Return table of results
cbind(apply(FDR_ARRAY,c(2,3),mean),apply(TPR_ARRAY,c(2,3),mean))

