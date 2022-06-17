rm(list=ls())
source('Subsampling.R')
library(cluster)
energy = read.csv('energydata_complete.csv')
colnames(energy);head(energy)
summary(energy$Appliances)
hist(log(energy$Appliances),freq=FALSE)
plot(log(energy$RH_2),log(energy$Appliances))
dat = as.matrix(energy[,c('Appliances','RH_1','RH_2','RH_3','RH_4','RH_5','RH_6','RH_7','RH_8','RH_9')])
unique(which(is.na(dat),arr.ind=TRUE)[,1])

n = dim(dat)[1]

y = log(dat[,'Appliances']); X = cbind(1, log(dat[,c('RH_1','RH_2','RH_3')]))

re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop  <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),dim(X)[2]*2,re))
vari.Aop  <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),2,re))
p.Aop  <- p.Lop  <- p.Aop.beta  <- p.Uni <- array(0,dim = c(length(subsize),2,re))
sig.est.Aop  <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),dim(X)[2]*2 + 3,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <- matrix(0, re, length(subsize))
time.full <- c()
set.seed(20210320)
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = X[ind,]; unif.y = y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(y, X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops(y, X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta(y, X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(y, X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta); beta.Aop[j,,i] = c(fit.Aop$beta)
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta); beta.Lop[j,,i] = c(fit.Lop$beta)  
      vari.Uni[j,,i] = fit.uni$sigma; vari.Aop[j,,i] = fit.Aop$sigma
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p; p.Aop[j,,i] = fit.Aop$p 
      p.Aop.beta[j,,i] = fit.Aop.beta$p; p.Lop[j,,i] = fit.Lop$p
      time.Uni[i,j] = fit.uni$time[3]; time.Aop[i,j] = fit.Aop$time[3]
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
}
set.seed(20210320)
ind = sample(1:n, 500, T); unif.x = X[ind,]; unif.y = y[ind] 
init = kmeans(unif.y, centers = 2, nstart = 25);
ind.c = order(init$centers)
ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
Full = FULL(y, X, beta.ini, vari.ini = vari.ini)
time.full = Full$time[3]
save.image("energy.RData")
#####
ind = order(Full$p)
beta = Full$beta[,ind]; sigma = Full$sigma[ind]; p = Full$p[ind]
mse1 <- mse2 <- mse3 <- mse4 <- matrix(0,re,4)
for(i in 1:re){
   mse1[i,] = colSums((t(beta.Aop[,,i]) - c(beta))^2) + colSums((t(sqrt(vari.Aop[,,i]))  - c(sqrt(sigma)))^2) + colSums((t(p.Aop[,-1,i])  - c(p)[-1])^2)
   mse2[i,] = colSums((t(beta.Aop.beta[,,i]) - c(beta))^2) + colSums((t(sqrt(vari.Aop.beta[,,i]))  - c(sqrt(sigma)))^2) + colSums((t(p.Aop.beta[,-1,i])  - c(p)[-1])^2)
   mse3[i,] = colSums((t(beta.Lop[,,i]) - c(beta))^2) + colSums((t(sqrt(vari.Lop[,,i]))  - c(sqrt(sigma)))^2) + colSums((t(p.Lop[,-1,i])  - c(p)[-1])^2)
   mse4[i,] = colSums((t(beta.Uni[,,i]) - c(beta))^2) + colSums((t(sqrt(vari.Uni[,,i]))  - c(sqrt(sigma)))^2) + colSums((t(p.Uni[,-1,i])  - c(p)[-1])^2)
}
colMeans(mse1)
colMeans(mse2)
colMeans(mse3)
colMeans(mse4)

library(xtable)
A = matrix(round(rowMeans(matrix(beta.Aop[3,,],8,1000)),3),4,2)[,ind]
Lbeta=matrix(round(rowMeans(matrix(beta.Aop.beta[3,,],8,1000)),3),4,2)[,ind]
L = matrix(round(rowMeans(matrix(beta.Lop[3,,],8,1000)),3),4,2)[,ind]
uni = matrix(round(rowMeans(matrix(beta.Uni[3,,],8,1000)),3),4,2)[,ind]
xtable(cbind(beta[,ind],A,Lbeta,L,uni),digits=3)


round(rowMeans(matrix(vari.Aop[3,,],2,1000)),3)
round(rowMeans(matrix(vari.Aop.beta[3,,],2,1000)),3)
round(rowMeans(matrix(vari.Lop[3,,],2,1000)),3)
round(rowMeans(matrix(vari.Uni[3,,],2,1000)),3)
Full$sigma
round(apply(p.Aop[3,,],1, mean),3)
round(apply(p.Aop.beta[3,,],1, mean),3)
round(apply(p.Lop[3,,],1, mean),3)
round(apply(p.Uni[3,,],1, mean),3)
Full$p

time.Uni[,3]
mean(time.Aop[,3])
mean(time.Aop.beta[,3] )
mean(time.Lop[,3])
time.full
