###############################################################################
################# 10^5  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^5; pi <- c(.5, .5)
beta <- matrix(c( 1, 1, 1, 1,1,1,1,1,1,1, 4, 4, 4, 4,4,4,4,4,4,4), 10, 2) ; sigma <- c(2,2); corr = 0.5
set.seed(20210320)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 100; subsize <- 2000 
beta.Aop <- beta.Lop  <- beta.Aop.beta  <- beta.Uni <-matrix(0,length(c(beta)),re)
vari.Aop <- vari.Lop  <- vari.Aop.beta  <- vari.Uni <- matrix(0,length(sigma),re)
p.Aop <- p.Aop1  <- p.Lop  <- p.Aop.beta  <- p.Uni <- 
   prop.Aop <- prop.Aop1 <- prop.Lop  <- prop.Aop.beta  <- prop.Uni <-  matrix(0,length(pi),re)

time.Aop <- time.Aop1 <- time.Aop.beta <- time.Lop <- time.Uni <- time.full <- c()
for(i in 1:re)
{
   ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
   init = kmeans(unif.y, centers = 2, nstart = 25);
   ind.c = order(init$centers)
   ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
   f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
   beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
   vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
   fit.uni = Uni(dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Aop = A.ops( dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Lop = L.ops(dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   beta.Uni[,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
   beta.Aop[,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
   beta.Aop.beta[,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); 
   beta.Lop[,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
   vari.Uni[,i] = fit.uni$sigma; 
   vari.Aop[,i] = fit.Aop$sigma; 
   vari.Aop.beta[,i] = fit.Aop.beta$sigma; 
   vari.Lop[,i] = fit.Lop$sigma
   p.Uni[,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
   p.Aop[,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
   p.Aop.beta[,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
   p.Lop[,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
   time.Uni[i] = fit.uni$time[3]; time.Aop[i] = fit.Aop$time[3];
   time.Aop.beta[i] = fit.Aop.beta$time[3]; time.Lop[i] = fit.Lop$time[3]
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]
}
###############################################################################
################# 5*10^5  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 5*10^5; pi <- c(.5, .5)
beta <- matrix(c( 1, 1, 1, 1,1,1,1,1,1,1, 4, 4, 4, 4,4,4,4,4,4,4), 10, 2); sigma <- c(2,2); corr = 0.5
set.seed(20210320)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 100; subsize <- 2000 
time.Aop <- time.Aop1 <- time.Aop.beta <- time.Lop <- time.Uni <- time.full <- c()
for(i in 1:re)
{
   ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
   init = kmeans(unif.y, centers = 2, nstart = 25);
   ind.c = order(init$centers)
   ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
   f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
   beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
   vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
   fit.uni = Uni(dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Aop = A.ops( dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Lop = L.ops(dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   beta.Uni[,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
   beta.Aop[,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
   beta.Aop.beta[,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); 
   beta.Lop[,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
   vari.Uni[,i] = fit.uni$sigma; 
   vari.Aop[,i] = fit.Aop$sigma; 
   vari.Aop.beta[,i] = fit.Aop.beta$sigma; 
   vari.Lop[,i] = fit.Lop$sigma
   p.Uni[,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
   p.Aop[,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
   p.Aop.beta[,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
   p.Lop[,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
   time.Uni[i] = fit.uni$time[3]; time.Aop[i] = fit.Aop$time[3]; 
   time.Aop.beta[i] = fit.Aop.beta$time[3]; time.Lop[i] = fit.Lop$time[3]
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]
}
###############################################################################
################# 10^6  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^6; pi <- c(.5, .5)
beta <- matrix(c( 1, 1, 1, 1,1,1,1,1,1,1, 4, 4, 4, 4,4,4,4,4,4,4), 10, 2) ; sigma <- c(2,2); corr = 0.5
set.seed(20210320)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 100; subsize <- 2000 
time.Aop <- time.Aop1 <- time.Aop.beta <- time.Lop <- time.Uni <- time.full <- c()
for(i in 1:re)
{
   ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
   init = kmeans(unif.y, centers = 2, nstart = 25);
   ind.c = order(init$centers)
   ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
   f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
   beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
   vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
   fit.uni = Uni(dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Aop = A.ops( dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Lop = L.ops(dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   beta.Uni[,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
   beta.Aop[,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
   beta.Aop.beta[,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); 
   beta.Lop[,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
   vari.Uni[,i] = fit.uni$sigma; 
   vari.Aop[,i] = fit.Aop$sigma; 
   vari.Aop.beta[,i] = fit.Aop.beta$sigma; 
   vari.Lop[,i] = fit.Lop$sigma
   p.Uni[,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
   p.Aop[,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
   p.Aop.beta[,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
   p.Lop[,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
   time.Uni[i] = fit.uni$time[3]; time.Aop[i] = fit.Aop$time[3]; 
   time.Aop.beta[i] = fit.Aop.beta$time[3]; time.Lop[i] = fit.Lop$time[3]
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]
}
###############################################################################
################# 5*10^6  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 5*10^6; pi <- c(.5, .5)
beta <- matrix(c( 1, 1, 1, 1,1,1,1,1,1,1, 4, 4, 4, 4,4,4,4,4,4,4), 10, 2) ; sigma <- c(2,2); corr = 0.5
set.seed(20210320)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 100; subsize <- 2000 
time.Aop <- time.Aop1 <- time.Aop.beta <- time.Lop <- time.Uni <- time.full <- c()
for(i in 1:re)
{
   ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
   init = kmeans(unif.y, centers = 2, nstart = 25);
   ind.c = order(init$centers)
   ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
   f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
   beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
   vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
   fit.uni = Uni(dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Aop = A.ops( dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Lop = L.ops(dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   beta.Uni[,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
   beta.Aop[,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
   beta.Aop.beta[,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); 
   beta.Lop[,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
   vari.Uni[,i] = fit.uni$sigma; 
   vari.Aop[,i] = fit.Aop$sigma; 
   vari.Aop.beta[,i] = fit.Aop.beta$sigma; 
   vari.Lop[,i] = fit.Lop$sigma
   p.Uni[,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
   p.Aop[,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
   p.Aop.beta[,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
   p.Lop[,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
   time.Uni[i] = fit.uni$time[3]; time.Aop[i] = fit.Aop$time[3];
   time.Aop.beta[i] = fit.Aop.beta$time[3]; time.Lop[i] = fit.Lop$time[3]
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]
}
###############################################################################
################# 10^7  ################
###############################################################################
rm(list=ls())
source('Subsampling.R');source('Datageneration.R')
library(cluster)
n <- 10^7; pi <- c(.5, .5)
beta <- matrix(c( 1, 1, 1, 1,1,1,1,1,1,1, 4, 4, 4, 4,4,4,4,4,4,4), 10, 2); sigma <- c(2,2); corr = 0.5
set.seed(20210320)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')
re = 100; subsize <- 2000 
iter.Aop <- iter.Aop1 <- iter.Lop  <- iter.Aop.beta  <- iter.Uni <- matrix(0, 2, re)
time.Aop <- time.Aop1 <- time.Aop.beta <- time.Lop <- time.Uni <- time.full <- iter.full <- c()
for(i in 1:re)
{
   ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 
   init = kmeans(unif.y, centers = 2, nstart = 25);
   ind.c = order(init$centers)
   ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
   f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
   beta.ini <- cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
   vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)
   fit.uni = Uni(dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Aop = A.ops( dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   fit.Lop = L.ops(dat$Y, dat$X, r = subsize, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
   beta.Uni[,i] = c(fit.uni$beta[,order(colSums(fit.uni$beta))]); 
   beta.Aop[,i] = c(fit.Aop$beta[,order(colSums(fit.Aop$beta))]); 
   beta.Aop.beta[,i] = c(fit.Aop.beta$beta[,order(colSums(fit.Aop.beta$beta))]); 
   beta.Lop[,i] = c(fit.Lop$beta[,order(colSums(fit.Lop$beta))])  
   vari.Uni[,i] = fit.uni$sigma; 
   vari.Aop[,i] = fit.Aop$sigma; 
   vari.Aop.beta[,i] = fit.Aop.beta$sigma; 
   vari.Lop[,i] = fit.Lop$sigma
   p.Uni[,i] = fit.uni$p[order(colSums(fit.uni$beta))]; 
   p.Aop[,i] = fit.Aop$p[order(colSums(fit.Aop$beta))] ; 
   p.Aop.beta[,i] = fit.Aop.beta$p[order(colSums(fit.Aop.beta$beta))] 
   p.Lop[,i] = fit.Lop$p[order(colSums(fit.Lop$beta))]
   time.Uni[i] = fit.uni$time[3]; time.Aop[i] = fit.Aop$time[3]; 
   time.Aop.beta[i] = fit.Aop.beta$time[3]; time.Lop[i] = fit.Lop$time[3]
   cat('re =', i, '\n')
   Full = FULL(dat$Y, dat$X, beta.ini, vari.ini = vari.ini)
   time.full[i] = Full$time[3]
}
