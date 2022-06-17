###############################################################################
################# Model 1 with sigma =(0.5,1) and p = (0.5,0.5)  ##############
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.5, .5)
beta <- matrix(c( 1, 1, 1, 1, 4, 4, 4, 4), 4, 2) ; sigma <- c(0.5,1); corr = 0.5
set.seed(20210320)
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
      beta.Uni[j,,i] = c(fit.uni$beta); beta.Aop[j,,i] = c(fit.Aop$beta); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta); beta.Lop[j,,i] = c(fit.Lop$beta)  
      vari.Uni[j,,i] = fit.uni$sigma
      vari.Aop[j,,i] = fit.Aop$sigma
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma  
      vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est
      sig.est.Aop[j,,i] = fit.Aop$sig.est
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est
      sig.est.Lop[j,,i] = fit.Lop$sig.est
      p.Uni[j,,i] = fit.uni$p
      p.Aop[j,,i] = fit.Aop$p
      p.Aop.beta[j,,i] = fit.Aop.beta$p
      p.Lop[j,,i] = fit.Lop$p
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
save.image("simul_Li_case1-1add.RData")
###############################################################################
################# Model 1 with sigma =(0.5,1) and p = (0.8,0.2)  ##############
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.8, .2)
beta <- matrix(c( 1, 1, 1, 1, 4, 4, 4, 4), 4, 2) ; sigma <- c(0.5,1); corr = 0.5
set.seed(20210320)
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
      beta.Uni[j,,i] = c(fit.uni$beta); beta.Aop[j,,i] = c(fit.Aop$beta); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta); beta.Lop[j,,i] = c(fit.Lop$beta)  
      vari.Uni[j,,i] = fit.uni$sigma
      vari.Aop[j,,i] = fit.Aop$sigma
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma  
      vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est
      sig.est.Aop[j,,i] = fit.Aop$sig.est
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est
      sig.est.Lop[j,,i] = fit.Lop$sig.est
      p.Uni[j,,i] = fit.uni$p
      p.Aop[j,,i] = fit.Aop$p
      p.Aop.beta[j,,i] = fit.Aop.beta$p
      p.Lop[j,,i] = fit.Lop$p
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
save.image("simul_Li_case2-1add.RData")
###############################################################################
################# Model 1 with sigma =(1,2) and p = (0.5,0.5)  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.5, .5)
beta <- matrix(c( 1, 1, 1, 1, 4, 4, 4, 4), 4, 2) ; sigma <- c(1,2); corr = 0.5
set.seed(20210320)
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
      beta.Uni[j,,i] = c(fit.uni$beta); beta.Aop[j,,i] = c(fit.Aop$beta); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta); beta.Lop[j,,i] = c(fit.Lop$beta)  
      vari.Uni[j,,i] = fit.uni$sigma
      vari.Aop[j,,i] = fit.Aop$sigma
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma  
      vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est
      sig.est.Aop[j,,i] = fit.Aop$sig.est
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est
      sig.est.Lop[j,,i] = fit.Lop$sig.est
      p.Uni[j,,i] = fit.uni$p
      p.Aop[j,,i] = fit.Aop$p
      p.Aop.beta[j,,i] = fit.Aop.beta$p
      p.Lop[j,,i] = fit.Lop$p
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
save.image("simul_Li_case1-2add.RData")
###############################################################################
################# Model 1 with sigma =(1,2) and p = (0.8,0.2)  ################
###############################################################################
n <- 10^5; pi <- c(.8, .2)
beta <- matrix(c( 1, 1, 1, 1, 4, 4, 4, 4), 4, 2) ;sigma <- c(1,2); corr = 0.5
set.seed(20210320)
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
      beta.Uni[j,,i] = c(fit.uni$beta); beta.Aop[j,,i] = c(fit.Aop$beta); 
      beta.Aop.beta[j,,i] = c(fit.Aop.beta$beta); beta.Lop[j,,i] = c(fit.Lop$beta)  
      vari.Uni[j,,i] = fit.uni$sigma
      vari.Aop[j,,i] = fit.Aop$sigma
      vari.Aop.beta[j,,i] = fit.Aop.beta$sigma  
      vari.Lop[j,,i] = fit.Lop$sigma
      sig.est.Uni[j,,i] = fit.uni$sig.est
      sig.est.Aop[j,,i] = fit.Aop$sig.est
      sig.est.Aop.beta[j,,i] = fit.Aop.beta$sig.est
      sig.est.Lop[j,,i] = fit.Lop$sig.est
      p.Uni[j,,i] = fit.uni$p
      p.Aop[j,,i] = fit.Aop$p
      p.Aop.beta[j,,i] = fit.Aop.beta$p
      p.Lop[j,,i] = fit.Lop$p
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
save.image("simul_Li_case2-2add.RData")
