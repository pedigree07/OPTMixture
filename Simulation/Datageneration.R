#######################################################################################
library(mvtnorm)
Datageneration = function(n, beta, sigma, pi, corr, Dist = 'mzNormal')
{
   d <- nrow(beta)-1;   K <- ncol(beta)
   sigmax <- matrix(corr, d, d) + diag(1-corr, d)
   if( Dist == 'mzNormal' ){X  <- rmvnorm(n, rep(0, d), sigmax)}
   X = cbind(1, X)
   sig <- diag(sigma, K)
   error <- rmvnorm(n, rep(0, K), sig)
   ymat <- X %*% beta + error # n by K matrix
   ind <- t(rmultinom(n, size = 1, prob = pi))
   Y <- rowSums(ymat * ind)
   list(Y = Y, X = X, ind = ind)
}   