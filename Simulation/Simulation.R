###############################################################################
################# Model 1 with sigma =(1,1) and p = (0.5,0.5)  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.5, .5)
beta <- matrix(c( 1, 1, 1, 1, 4, 4, 4, 4), 4, 2) ; sigma <- c(1,1); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop1 <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
      beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); 
      beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; 
      vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; 
      vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; 
      sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; 
      sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3];
      time.Aop[i,j] = fit.Aop$time[3];
      time.Aop.beta[i,j] = fit.Aop.beta$time[3];
      time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case1.RData")
##############################################################################
################# Model 1 with sigma =(1,1) and p = (0.8,0.2)  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.8, .2)
beta <- matrix(c( 1, 1, 1, 1, 4, 4, 4, 4), 4, 2) ; sigma <- c(1,1); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
      beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); 
      beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; 
      vari.Aop[j,,i] = fit.Aop$sigma;
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; 
      vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; 
      sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; 
      sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; 
      time.Aop[i,j] = fit.Aop$time[3]; 
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; 
      time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case2.RData")
###############################################################################
################# Model 2 with sigma =(1,1) and p = (0.5,0.5)  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.5, .5)
beta <- matrix(c( -4, -4, -4, -4, 1, 1, 1, 1), 4, 2) ; sigma <- c(1,1); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
      beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); 
      beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; 
      vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; 
      vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; 
      sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; 
      sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; 
      time.Aop[i,j] = fit.Aop$time[3]; 
      time.Aop.beta[i,j] = fit.Aop.beta$time[3];
      time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case3.RData")
###############################################################################
################# Model 2 with sigma =(1,1) and p = (0.8,0.2)  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.2, .8)
beta <- matrix(c( -4, -4, -4, -4, 1, 1, 1, 1), 4, 2) ; sigma <- c(1,1); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
      beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; 
      vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; 
      vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; 
      sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; 
      time.Aop[i,j] = fit.Aop$time[3]; 
      time.Aop.beta[i,j] = fit.Aop.beta$time[3];
      time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case4.RData")
###############################################################################
############## Model 3 with sigma =(1,1,1) and p = (1/3,1/3,1/3)  #############
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(1/3, 1/3, 1/3)
beta <- matrix(c(-4, -4, -4, -4, 1, 1, 1, 1, 4, 4, 4, 4), 4, 3) ; 
sigma <- c(1,1,1); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop  <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 3, nstart = 25);
      ind.c = order(init$centers);init$centers
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]);
      ind.t = which(init$cluster == ind.c[3]);
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      t.lm = lm(unif.y[ind.t]~unif.x[ind.t,]-1)
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual,
                   sum(t.lm$residuals^2)/t.lm$df.residual)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients, t.lm$coefficients);beta.ini
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
      beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; 
      vari.Aop[j,,i] = fit.Aop$sigma;
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; 
      vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; 
      sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; 
      sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; 
      time.Aop[i,j] = fit.Aop$time[3]; 
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; 
      time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case5.RData")
###############################################################################
############# Model 3 with sigma =(1,1,1) and p = (0.25,0.5,0.25)  ############
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(0.25, 0.5, 0.25)
beta <- matrix(c(-4, -4, -4, -4, 1, 1, 1, 1, 4, 4, 4, 4), 4, 3) ; 
sigma <- c(1,1,1); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 3, nstart = 25);
      ind.c = order(init$centers);init$centers
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]);
      ind.t = which(init$cluster == ind.c[3]);
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      t.lm = lm(unif.y[ind.t]~unif.x[ind.t,]-1)
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual,
                   sum(t.lm$residuals^2)/t.lm$df.residual)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients, t.lm$coefficients);beta.ini
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; time.Aop[i,j] = fit.Aop$time[3]; 
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case6.RData")
###############################################################################
############### Model 1 with sigma =(0.5,0.5) and p = (0.5,0.5)  ##############
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.5, .5)
beta <- matrix(c( 1, 1, 1, 1, 4, 4, 4, 4), 4, 2) ; sigma <- c(0.5,0.5); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
      beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; 
      vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; 
      vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; 
      sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; 
      sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; 
      time.Aop[i,j] = fit.Aop$time[3]; 
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; 
      time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case1-1.RData")
###############################################################################
############### Model 1 with sigma =(0.5,0.5) and p = (0.8,0.2)  ##############
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.8, .2)
beta <- matrix(c( 1, 1, 1, 1, 4, 4, 4, 4), 4, 2) ; sigma <- c(0.5,0.5); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
      beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; 
      vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; 
      sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; 
      sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; 
      time.Aop[i,j] = fit.Aop$time[3]; 
      time.Aop.beta[i,j] = fit.Aop.beta$time[3];
      time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case2-1.RData")
###############################################################################
############### Model 2 with sigma =(0.5,0.5) and p = (0.5,0.5)  ##############
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.5, .5)
beta <- matrix(c( -4, -4, -4, -4, 1, 1, 1, 1), 4, 2) ; sigma <- c(0.5,0.5); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
      beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; 
      vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; 
      sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; 
      sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; 
      time.Aop[i,j] = fit.Aop$time[3]; 
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; 
      time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case3-1.RData")
###############################################################################
############### Model 2 with sigma =(0.5,0.5) and p = (0.2,0.8)  ##############
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.2, .8)
beta <- matrix(c(-4, -4, -4, -4, 1, 1, 1, 1), 4, 2) ; sigma <- c(0.5,0.5); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
      beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; 
      vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; 
      vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; 
      sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; 
      sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; 
      time.Aop[i,j] = fit.Aop$time[3]; 
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; 
      time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case4-1.RData")
###############################################################################
########### Model 3 with sigma =(0.5,0.5,0.5) and p = (1/3,1/3,1/3)  ##########
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(1/3, 1/3, 1/3)
beta <- matrix(c(-4, -4, -4, -4, 1, 1, 1, 1, 4, 4, 4, 4), 4, 3) ; 
sigma <- c(0.5,0.5,0.5); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 3, nstart = 25);
      ind.c = order(init$centers);init$centers
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]);
      ind.t = which(init$cluster == ind.c[3]);
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      t.lm = lm(unif.y[ind.t]~unif.x[ind.t,]-1)
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual,
                   sum(t.lm$residuals^2)/t.lm$df.residual)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients, t.lm$coefficients);beta.ini
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
      beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); 
      beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; 
      vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; 
      vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; 
      sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; 
      sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; 
      time.Aop[i,j] = fit.Aop$time[3]; 
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; 
      time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case5-1.RData")
###############################################################################
########## Model 3 with sigma =(0.5,0.5,0.5) and p = (0.25,0.5,0.25)  #########
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(0.25, 0.5, 0.25)
beta <- matrix(c(-4, -4, -4, -4, 1, 1, 1, 1, 4, 4, 4, 4), 4, 3) ; 
sigma <- c(0.5,0.5,0.5); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop  <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 3, nstart = 25);
      ind.c = order(init$centers);init$centers
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]);
      ind.t = which(init$cluster == ind.c[3]);
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      t.lm = lm(unif.y[ind.t]~unif.x[ind.t,]-1)
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual,
                   sum(t.lm$residuals^2)/t.lm$df.residual)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients, t.lm$coefficients);beta.ini
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
      beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); 
      beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; 
      sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; 
      sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; 
      time.Aop[i,j] = fit.Aop$time[3];
      time.Aop.beta[i,j] = fit.Aop.beta$time[3];
      time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case6-1.RData")
###############################################################################
################# Model 1 with sigma =(2,2) and p = (0.5,0.5)  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.5, .5)
beta <- matrix(c( 1, 1, 1, 1, 4, 4, 4, 4), 4, 2) ; sigma <- c(2,2); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; time.Aop[i,j] = fit.Aop$time[3]; 
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case1-2.RData")
###############################################################################
################# Model 1 with sigma =(2,2) and p = (0.8,0.2)  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.8, .2)
beta <- matrix(c( 1, 1, 1, 1, 4, 4, 4, 4), 4, 2) ;sigma <- c(2,2); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; time.Aop[i,j] = fit.Aop$time[3];
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case2-2.RData")
###############################################################################
################# Model 2 with sigma =(2,2) and p = (0.5,0.5)  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.5, .5)
beta <- matrix(c( -4, -4, -4, -4, 1, 1, 1, 1), 4, 2) ; sigma <- c(2,2); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop  <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; vari.Aop[j,,i] = fit.Aop$sigma;
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; sig.est.Aop[j,,i] = fit.Aop$sig.est;
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; time.Aop[i,j] = fit.Aop$time[3];
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case3-2.RData")
###############################################################################
################# Model 2 with sigma =(2,2) and p = (0.2,0.8)  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.2, .8)
beta <- matrix(c( -4, -4, -4, -4, 1, 1, 1, 1), 4, 2) ;sigma <- c(2,2); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop  <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop  <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop  <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 2, nstart = 25);
      ind.c = order(init$centers)
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; time.Aop[i,j] = fit.Aop$time[3]; 
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case4-2.RData")
###############################################################################
################# Model 3 with sigma =(2,2,2) and p = (1/3,1/3,1/3)  ##########
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(1/3, 1/3, 1/3)
beta <- matrix(c(-4, -4, -4, -4, 1, 1, 1, 1, 4, 4, 4, 4), 4, 3) ; 
sigma <- c(2,2,2); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop  <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop  <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop  <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop  <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop  <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop  <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 3, nstart = 25);
      ind.c = order(init$centers);init$centers
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]);
      ind.t = which(init$cluster == ind.c[3]);
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      t.lm = lm(unif.y[ind.t]~unif.x[ind.t,]-1)
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual,
                   sum(t.lm$residuals^2)/t.lm$df.residual)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients, t.lm$coefficients);beta.ini
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; vari.Aop[j,,i] = fit.Aop$sigma;
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; time.Aop[i,j] = fit.Aop$time[3];
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case5-2.RData")
###############################################################################
########## Model 3 with sigma =(2,2,2) and p = (0.25,0.5,0.25)  #########
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(0.25, 0.5, 0.25)
beta <- matrix(c(-4, -4, -4, -4, 1, 1, 1, 1, 4, 4, 4, 4), 4, 3) ; 
sigma <- c(2,2,2); corr = 0.5
set.seed(20220613)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 1000; subsize <- c(500, 1000, 1500, 2000) 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <- array(0,dim = c(length(subsize),length(c(beta)),re))
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- array(0,dim = c(length(subsize),length(sigma),re))
p.Aop <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <- array(0,dim = c(length(subsize),length(pi),re))
iter.Aop <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <-array(0,dim = c(length(subsize),2,re))
sig.est.Aop <- sig.est.Lop  <- sig.est.Aop.beta  <- sig.est.Uni <- 
   array(0,dim = c(length(subsize),length(c(beta)) + length(sigma) + length(pi) - 1,re))
time.Aop <- time.Aop.beta <- time.Lop <- time.Uni <-matrix(0, re, length(subsize))
Full <- matrix(0, re, length(c(beta))); time.full <- iter.full <- c()
for(i in 1:re)
{
   for(j in 1:length(subsize))
   {
      ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
      init = kmeans(unif.y, centers = 3, nstart = 25);
      ind.c = order(init$centers);init$centers
      ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]);
      ind.t = which(init$cluster == ind.c[3]);
      f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
      t.lm = lm(unif.y[ind.t]~unif.x[ind.t,]-1)
      vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual,
                   sum(t.lm$residuals^2)/t.lm$df.residual)
      beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients, t.lm$coefficients);beta.ini
      fit.uni = Uni(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop = A.ops( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      fit.Lop = L.ops(dat$Y, dat$X, r = subsize[j], pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
      beta.Uni[j,,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); beta.Aop[j,,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); beta.Lop[j,,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
      vari.Uni[j,,i] = fit.uni$sigma; vari.Aop[j,,i] = fit.Aop$sigma; 
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma; vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est; sig.est.Aop[j,,i] = fit.Aop$sig.est; 
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est; sig.est.Lop[j,,i] = fit.Lop$sig.est; 
      p.Uni[j,,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
      p.Aop[j,,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
      p.Aop.beta[j,,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
      p.Lop[j,,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
      prop.Aop[j,,i] = colMeans(dat$ind[c(fit.Aop$ind,ind),])
      prop.Aop.beta[j,,i] = colMeans(dat$ind[c(fit.Aop.beta$ind,ind),])
      prop.Lop[j,,i] =  colMeans(dat$ind[c(fit.Lop$ind,ind),])
      prop.Uni[j,,i] = colMeans(dat$ind[fit.uni$ind,])
      iter.Aop[j,,i] =  c(fit.Aop$first.iter, fit.Aop$iter)
      iter.Aop.beta[j,,i] = c(fit.Aop.beta$first.iter, fit.Aop.beta$iter)
      iter.Lop[j,,i] =  c(fit.Lop$first.iter, fit.Lop$iter)
      iter.Uni[j,,i] = fit.uni$iter
      time.Uni[i,j] = fit.uni$time[3]; time.Aop[i,j] = fit.Aop$time[3];
      time.Aop.beta[i,j] = fit.Aop.beta$time[3]; time.Lop[i,j] = fit.Lop$time[3]
      print(j)
   }
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]; iter.full[i] = Full$iter
}
save.image("simul_Li_case6-2.RData")



