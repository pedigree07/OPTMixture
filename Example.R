source('Subsampling.R')
source('Datageneration.R')
library(cluster)
###################################################################
######################### Data Generation #########################
###################################################################
n <- 10^5; pi <- c(.5, .5)
#pi <- c(.3, .4, .3)
beta <- matrix(c( 1, 1, 1, 1, 4, 4, 4, 4), 4, 2) ; sigma <- c(0.5,0.5); corr = 0.5
set.seed(20210408)
dat = Datageneration(n, beta, sigma, pi, corr, Dist = 'mzNormal')

##################### Draw a pilot sample of size 500 ###########
ind = sample(1:n, 500, T); unif.x = dat$X[ind,]; unif.y = dat$Y[ind] 

######################### Obtain the initial values  ###########
######################### using  k-means clustering  ###########
init = kmeans(unif.y, centers = 2, nstart = 25);
ind.c = order(init$centers)
ind.f = which(init$cluster == ind.c[1]); ind.s = which(init$cluster == ind.c[2]); 
f.lm = lm(unif.y[ind.f]~unif.x[ind.f,]-1); s.lm = lm(unif.y[ind.s]~unif.x[ind.s,]-1)
beta.ini = cbind(f.lm$coefficients, s.lm$coefficients); beta.ini
vari.ini = c(sum(f.lm$residuals^2)/f.lm$df.residual, sum(s.lm$residuals^2)/s.lm$df.residual)

############# Subsample-based Estimation with subsample size 1000 #################
# using optimal subsampling probabilities minimizing the trace of V
fit.Aop = A.ops( dat$Y, dat$X, r = 1000, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
# using optimal subsampling probabilities minimizing the trace of V(beta)
fit.Aop.beta = A.ops.beta( dat$Y, dat$X, r = 1000, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)
# using optimal subsampling probabilities minimizing the trace of V(pi)
fit.Lop = L.ops(dat$Y, dat$X, r = 1000, pilot.ind = ind, beta.ini = beta.ini, vari.ini = vari.ini)

############# Output #################
#fit.Aop$beta : coefficient estimates
#fit.Aop$sigma : sd estimates
#fit.Aop$p : component weight estimates
