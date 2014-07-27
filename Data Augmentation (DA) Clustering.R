#required packages
library(mvtnorm)
library(MCMCpack)
library(calibrate)

data=read.csv("health1.csv")
data1 = data_com[1,,][,2:12]
mu=rep(0,11)
for (i in 1:11){
  mu[i]=mean(data1[,i])
}
n=nrow(data1)   #46

z=rep(0,n)
mu1=matrix(rep(0,451),41,11)
mu2=matrix(rep(0,451),41,11)

sigma1=array(0,dim=c(41,11,11))
sigma2=array(0,dim=c(41,11,11))

lambda1=c(1/2,rep(0,40))
lambda2=c(1/2,rep(0,40))  # initials lambda1=lambda2=lambda3=1/3
n1=rep(0,40)
n2=rep(0,40)


mu1[1,]=mu-1
mu2[1,]=mu

sigma1[1,,]=cov(data1)
sigma2[1,,]=cov(data1) #set initial values

for(t in 1:20){  #10 iterations
  for(i in 1:n){
    a=lambda1[t]*dmvnorm(data1[i,], mean=mu1[t,], sigma=sigma1[t,,])
    b=lambda2[t]*dmvnorm(data1[i,], mean=mu2[t,], sigma=sigma2[t,,])
    z[i]=sample(c(1,2),size=1,prob=c(a/(a+b),b/(a+b)))
  }
  
  n1[t]=sum(z==1)
  n2[t]=sum(z==2)
  rdir=rdirichlet(1,alpha=c(n1[t]+1,n2[t]+1))
  lambda1[t+1]=rdir[1]
  lambda2[t+1]=rdir[2]
  y1_bar=apply(data1[z==1,],2,mean)
  y2_bar=apply(data1[z==2,],2,mean)
  s1=cov(data1[z==1,])
  s2=cov(data1[z==2,])
  
  
  sigma1[t+1,,]=riwish(v=n1[t]-1, S=s1)
  sigma2[t+1,,]=riwish(v=n1[t]-1, S=s2)

  mu1[t+1,]=rmvnorm(1, mean = y1_bar, sigma = sigma1[t+1,,]/n1[t])
  mu2[t+1,]=rmvnorm(1, mean = y2_bar, sigma = sigma2[t+1,,]/n2[t])
}

#### result ####
z
z1=c()
z2=c()
for (i in 1:46) {
  if (z[i]==1) {z1[i]=i}
  if (z[i]==2) {z2[i]=i}
}
z1=z1[!is.na(z1)] # group1
z2=z2[!is.na(z2)] # group2
z1;z2  

z1=c(8,12,15,19,28,33,34,36,38,43,44,46)
z2=c(1,2,3,4,5,6,7,9,10,11,13,14,16,17,18,20,21,22,23,24,25,26,27,29,30,31,32,35,37,39,40,41,42,45)

i=1
g1=lm(MR[,i][z1]~Pneumonia[,i][z1]+Diarrhoea[,i][z1]+Prematurity[,i][z1]
            +Neonatal[,i][z1]+Birth[,i][z1]+Water[,i][z1]+Health[,i][z1]
            +Antibiotic[,i][z1]+Vitamin[,i][z1]+Underweight[,i][z1])
g2=lm(MR[,i][z2]~Pneumonia[,i][z2]+Diarrhoea[,i][z2]+Prematurity[,i][z2]
            +Neonatal[,i][z2]+Birth[,i][z2]+Water[,i][z2]+Health[,i][z2]
            +Antibiotic[,i][z2]+Vitamin[,i][z2]+Underweight[,i][z2])

search1=step(g1, trace=0)
search1$anova
summary(search1)

search2=step(g2, trace=0)
search2$anova
summary(search2)