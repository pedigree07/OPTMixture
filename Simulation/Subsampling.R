EM.alg = function(x, y, pi, num.c, beta.ini, vari.ini, iter.max = 100, eps = 1e-6)
{   
   J = num.c; n = dim(x)[1]; p = dim(x)[2]
   cur.p = rep(1/J, J); cur.beta <- new.beta <- beta.ini
   new.vari <- cur.vari <- vari.ini; #w = matrix(0,n,J)
   iter = 1
   while(iter <= iter.max)
   {   
      w = t(dnorm(y - x %*% cur.beta,0,sqrt(cur.vari))) * cur.p
      #ind = rowSums(w) != 0
      w =  t(w)/colSums(w)
      #w[!ind,] = 0
      new.p = colSums(w/pi)/(sum(1/pi))
      for(j in 1:J)
      {   
         eig.Aos <- eigen( t(x) %*% (x *(w[,j]/pi)) )
         iI.Aos <- eig.Aos$vectors %*% (t(eig.Aos$vectors)/eig.Aos$values)
         #new.beta[,j] = solve(t(x) %*% (x *(w[,j]/pi)), t(x) %*% (y*(w[,j]/pi)) )
         new.beta[,j] =iI.Aos %*% t(x) %*% (y*(w[,j]/pi))
         new.vari[j] = sum((y - c(x%*% new.beta[,j]))^2*(w[,j]/pi))/(sum(w[,j]/pi))
      }
      tol = sum((new.p -cur.p)^2) + sum((c(new.beta)-c(cur.beta))^2) + sum((new.vari-cur.vari)^2)
      cur.p = new.p; cur.beta = new.beta; cur.vari = new.vari
      if( sqrt(tol) < eps) break
      iter = iter + 1
   }
   return(list(p = cur.p, beta = cur.beta, vari = cur.vari, w = w, iter = iter))
}
#######################################################################################
M.v = function(x, SP, w, res, v, w1t, w2t, w2t.w, w3t, vari = vari, pi = p, n, p, J)
{   
   vv = matrix(0,J*(p+1),J*(p+1)); 
   for(j in 1:J)  
   {   
      vv[(p*(j-1)+1):(p*j),(p*(j-1)+1):(p*j)] = t(x) %*% (x * w2t[j,]/(sqrt(vari[j])* SP)) 
      vv[(p*(j-1)+1):(p*j),(p*J+j)] = colSums(x * w1t[j,]* w2t.w[j,]/(SP) )
   }   
   vv22 = rowSums( w1t * (w2t.w - 2/sqrt(vari)) * t(res/SP)/sqrt(vari) )
   M = -v%*%(t(v)/SP)
   M[1:(J*p),1:(J*p)] = M[1:(J*p),1:(J*p)] + vv[1:(J*p),1:(J*p)];
   M[((J*p)+1):(J*(p+1)),1:(J*p)] <-M[1:(J*p),((J*p)+1):(J*(p+1))] + vv[1:(J*p),(J*p+1):(J*(p+1))]
   M[(J*p+1):(J*(p+1)), 1:(J*p)] <- t(M[1:(J*p),((J*p)+1):(J*(p+1))])
   M[((J*p)+1):(J*(p+1)),((J*p)+1):(J*(p+1))] = M[((J*p)+1):(J*(p+1)),((J*p)+1):(J*(p+1))] + diag(vv22)
   return(list(M = M, v = v))
} 
#######################################################################################
M.v1 = function(v, SP)
{   
   M = -v%*%(t(v)/SP)
   return(list(M = M))
}
#######################################################################################
M.v2 = function(x, SP, w, res, w1t, w2t, w2t.w, w3t, vari = vari, pi = p, n, p, J)
{   
   vv = matrix(0,J*(p+1),J*(p+1)); v = matrix(0,J*p+2*J-1,n)
   for(j in 1:J)  
   {   
      v[(p*(j-1)+1):(p*j),] = t(x *w1t[j,])
      vv[(p*(j-1)+1):(p*j),(p*(j-1)+1):(p*j)] = t(x) %*% (x * w2t[j,]/(sqrt(vari[j])* SP)) 
      vv[(p*(j-1)+1):(p*j),(p*J+j)] = colSums(x * w1t[j,]* w2t.w[j,]/(SP) )
   }   
   vv22 = rowSums( w1t * (w2t.w - 2/sqrt(vari)) * t(res/SP)/sqrt(vari) )
   v[((J*p)+1):(J*(p+1)),] = w2t;v[(J*(p+1)+1):(J*p+2*J-1),] = w3t
   #v = rbind(v1,w2t,w3t); 
   M = -v%*%(t(v)/SP)
   M[1:(J*p),1:(J*p)] = M[1:(J*p),1:(J*p)] + vv[1:(J*p),1:(J*p)];
   M[((J*p)+1):(J*(p+1)),1:(J*p)] <-M[1:(J*p),((J*p)+1):(J*(p+1))] + vv[1:(J*p),(J*p+1):(J*(p+1))]
   M[(J*p+1):(J*(p+1)), 1:(J*p)] <- t(M[1:(J*p),((J*p)+1):(J*(p+1))])
   M[((J*p)+1):(J*(p+1)),((J*p)+1):(J*(p+1))] = M[((J*p)+1):(J*(p+1)),((J*p)+1):(J*(p+1))] + diag(vv22)
   return(list(M = M, v = v))
} 
#######################################################################################
ini.theta = function(Y, X, k, init)
{   
   beta.ini = matrix(0, dim(X)[2], k)
   vari.ini = c()
   ind.c = order(init$centers)
   for(i in 1:k)
   {   
      ind = which(init$cluster == ind.c[i])
      fit = lm(Y[ind]~X[ind,]-1)
      beta.ini[,i] <- fit$coefficients
      vari.ini[i] = sum(fit$residuals^2)/fit$df.residual
   }
   return(list(beta.ini = beta.ini, vari.ini= vari.ini))
}
#######################################################################################
Uni = function(y, x, r, pilot.ind, beta.ini, vari.ini)
{
   n = dim(x)[1]; p = dim(x)[2]; J = dim(beta.ini)[2];sigma.Uni = matrix(0,p,J)
   r0 = length(pilot.ind)
   time.unif = system.time({
      ind = sample(1:n, r, T)
      unif.x = x[c(ind,pilot.ind),]; unif.y = y[c(ind,pilot.ind)] 
      w.Unif = rep(1/n, dim(unif.x)[1])
      fit.u = EM.alg(unif.x, unif.y, w.Unif, num.c = J, beta.ini = beta.ini, vari.ini = vari.ini)
   })
   ind = order(colSums(fit.u$beta))
   w = fit.u$w[,ind]; vari = fit.u$vari[ind]; beta = fit.u$beta[,ind];pr = fit.u$p[ind]
   Epsil_Unif = (unif.y - unif.x%*%beta)
   w1t = t(w * Epsil_Unif) / vari; 
   w2t.w = (t(Epsil_Unif^2) - vari)/(vari^(3/2)); 
   w2t = t(w) * w2t.w
   w3t = t(w[,-J])/pr[-J] - w[,J]/pr[J]
   A.M = M.v2(unif.x, w.Unif, w, Epsil_Unif, w1t, w2t, w2t.w, w3t, vari = vari, pi = pr, r + r0, p, J)
   iM.Unif = solve(A.M$M) * ( n * (r + r0) )
   V_Unif = ( A.M$v %*% (t(A.M$v)/(w.Unif^2)) ) / (( n * (r + r0) )^2) 
   sigma.Uni = diag(iM.Unif %*% V_Unif %*% iM.Unif)
   return(list(beta = beta, sigma = vari, V = V_Unif,iter = fit.u$iter, x = unif.x, y = unif.y,
               pi = w.Unif, iM = iM.Unif, w = w, p = pr , sig.est = sigma.Uni, time = time.unif, ind = ind))
} 
#######################################################################################
A.ops = function(y, x, r, pilot.ind, beta.ini, vari.ini)
{
   n = dim(x)[1]; p = dim(x)[2]; J = dim(beta.ini)[2]; v = matrix(0,J*p+2*J-1,n)
   p.Aop = rep(0,n); w = matrix(0, n, J); r0 = length(pilot.ind)
   time.Aop = system.time({
      pilot.x = x[pilot.ind,]; pilot.y = y[pilot.ind] 
      fit.Aop = EM.alg(pilot.x, pilot.y, rep(1/n,r0), num.c = J, beta.ini = beta.ini, vari.ini = vari.ini)
      
      res = (y - x%*%fit.Aop$beta)
      w =  t(dnorm(res,0,sqrt(fit.Aop$vari))) * fit.Aop$p
      w =  t(w)/colSums(w)
      w1t = t(w * res) / fit.Aop$vari; 
      w2t.w = (t(res^2) - fit.Aop$vari)/(fit.Aop$vari^(3/2))
      for(j in 1:J)     
         v[(p*(j-1)+1):(p*j),] = t(x *w1t[j,])
      v[((J*p)+1):(J*(p+1)),] =  t(w) * w2t.w 
      v[(J*(p+1)+1):(J*p+2*J-1),] = t(w[,-J])/fit.Aop$p[-J] - w[,J]/fit.Aop$p[J]
      M.vc = M.v(x[pilot.ind,], 1/n, w[pilot.ind,], res[pilot.ind,], v[,pilot.ind], w1t[,pilot.ind], w2t = v[((J*p)+1):(J*(p+1)),pilot.ind], 
                 w2t.w[,pilot.ind], w3t = v[(J*(p+1)+1):(J*p+2*J-1),pilot.ind], vari = fit.Aop$vari, pi = fit.Aop$p, r0, p, J)
      MvcM = (solve(M.vc$M/n) %*% v)*(r0); p.Aop = colSums(MvcM*MvcM)
      p.Aop <- sqrt(p.Aop)/sum(sqrt(p.Aop))
      idx.Aop <- sample(1:n, r, T, p.Aop)
      x.Aop <- x[c(idx.Aop,pilot.ind),]
      y.Aop <- y[c(idx.Aop,pilot.ind)]
      w.Aop <- c(p.Aop[idx.Aop],rep(1/n,r0))
      fit.A = EM.alg(x.Aop, y.Aop, w.Aop, num.c = J, beta.ini = fit.Aop$beta, vari.ini = fit.Aop$vari)
   })
   ind = order(colSums(fit.A$beta))
   w = fit.A$w[,ind]; vari = fit.A$vari[ind]; beta = fit.A$beta[,ind];pr = fit.A$p[ind]
   Epsil_Aop = y.Aop - x.Aop%*%beta
   w1t = t(w * Epsil_Aop) / vari; 
   w2t.w = (t(Epsil_Aop^2) - vari)/(vari^(3/2))
   w2t = t(w) * w2t.w; 
   w3t = t(w[,-J])/pr[-J] - w[,J]/pr[J]
   A.M = M.v2(x.Aop, w.Aop, w, Epsil_Aop, w1t, w2t, w2t.w, w3t, vari = vari, pi = pr, r0+r, p, J)
   iM.Aop = solve(A.M$M) * ( n * (r+r0) )
   V_Aop = ( A.M$v %*% (t(A.M$v)/(w.Aop^2)) ) / (( n * (r+r0) )^2) 
   sigma.Aop = diag(iM.Aop %*% V_Aop %*% iM.Aop)
   return(list(beta = beta, sigma = vari, V = V_Aop, first.iter = fit.Aop$iter, iter = fit.A$iter, 
               pi =  w.Aop, x = x.Aop, y = y.Aop, iM = iM.Aop,  w = w, p = pr , sig.est = sigma.Aop, time = time.Aop, ind= idx.Aop))
}
#######################################################################################
A.ops.beta = function(y, x, r, pilot.ind, beta.ini, vari.ini)
{
   n = dim(x)[1]; p = dim(x)[2]; J = dim(beta.ini)[2]; v = matrix(0,J*p+2*J-1,n)
   p.Aop = rep(0,n); w = matrix(0, n, J); r0 = length(pilot.ind)
   time.Aop = system.time({
      pilot.x = x[pilot.ind,]; pilot.y = y[pilot.ind] 
      fit.Aop = EM.alg(pilot.x, pilot.y, rep(1/n,r0), num.c = J, beta.ini = beta.ini, vari.ini = vari.ini)
      
      res = (y - x%*%fit.Aop$beta)
      w =  t(dnorm(res,0,sqrt(fit.Aop$vari))) * fit.Aop$p
      w =  t(w)/colSums(w)
      w1t = t(w * res) / fit.Aop$vari; 
      w2t.w = (t(res^2) - fit.Aop$vari)/(fit.Aop$vari^(3/2))
      for(j in 1:J)     
         v[(p*(j-1)+1):(p*j),] = t(x *w1t[j,])
      v[((J*p)+1):(J*(p+1)),] =  t(w) * w2t.w; 
      v[(J*(p+1)+1):(J*p+2*J-1),] = t(w[,-J])/fit.Aop$p[-J] - w[,J]/fit.Aop$p[J]
      M.vc = M.v(x[pilot.ind,], 1/n, w[pilot.ind,], res[pilot.ind,], v[,pilot.ind], w1t[,pilot.ind], w2t = v[((J*p)+1):(J*(p+1)),pilot.ind], 
                 w2t.w[,pilot.ind], w3t = v[(J*(p+1)+1):(J*p+2*J-1),pilot.ind], vari = fit.Aop$vari, pi = fit.Aop$p, r0, p, J)
      MvcM = solve(M.vc$M/n)*r0; M.beta = MvcM[1:(J*p),]
      M.b = M.beta %*% v 
      p.Aop = colSums(M.b*M.b)
      p.Aop <- sqrt(p.Aop)/sum(sqrt(p.Aop))
      idx.Aop <- sample(1:n, r, T, p.Aop)
      x.Aop <- x[c(idx.Aop,pilot.ind),]
      y.Aop <- y[c(idx.Aop,pilot.ind)]
      w.Aop <- c(p.Aop[idx.Aop],rep(1/n,r0))
      fit.A = EM.alg(x.Aop, y.Aop, w.Aop, num.c = J, beta.ini = fit.Aop$beta, vari.ini = fit.Aop$vari)
   })
   ind = order(colSums(fit.A$beta))
   w = fit.A$w[,ind]; vari = fit.A$vari[ind]; beta = fit.A$beta[,ind];pr = fit.A$p[ind]
   Epsil_Aop = y.Aop - x.Aop%*%beta
   w1t = t(w * Epsil_Aop) / vari; 
   w2t.w = (t(Epsil_Aop^2) - vari)/(vari^(3/2)); 
   w2t = t(w) * w2t.w
   w3t = t(w[,-J])/pr[-J] - w[,J]/pr[J]
   A.M = M.v2(x.Aop, w.Aop, w, Epsil_Aop, w1t, w2t,w2t.w, w3t, vari = vari, pi = pr, r0+r, p, J)
   iM.Aop = solve(A.M$M) * ( n * (r+r0) )
   V_Aop = ( A.M$v %*% (t(A.M$v)/(w.Aop^2)) ) / (( n * (r+r0) )^2) 
   sigma.Aop = diag(iM.Aop %*% V_Aop %*% iM.Aop)
   return(list(beta = beta, sigma = vari, V = V_Aop, first.iter = fit.Aop$iter, iter = fit.A$iter,
               pi =  w.Aop, x = x.Aop, y.Aop = y.Aop, iM = iM.Aop, w = w, p = pr , sig.est = sigma.Aop, time = time.Aop, ind = idx.Aop))
}
#######################################################################################
L.ops = function(y, x, r, pilot.ind, beta.ini, vari.ini)
{
   n = dim(x)[1]; p = dim(x)[2]; J = dim(beta.ini)[2];
   p.Lop = rep(0,n); w = matrix(0, n, J); r0 = length(pilot.ind)
   time.Lop = system.time({
      pilot.x = x[pilot.ind,]; pilot.y = y[pilot.ind] 
      fit.Lop = EM.alg(pilot.x, pilot.y, rep(1/n,r0), num.c = J, beta.ini = beta.ini, vari.ini = vari.ini)
      
      res = (y - x%*%fit.Lop$beta)
      w = t(dnorm(res,0,sqrt(fit.Lop$vari))) * fit.Lop$p
      w =  t(w)/colSums(w)
      w1t = t(w*res)/(fit.Lop$vari)
      w2t.w = (t(res)^2 - fit.Lop$vari)/(fit.Lop$vari^(3/2))
      w2t = t(w) * w2t.w
      w3t = t(w[,-J])/fit.Lop$p[-J] - w[,J]/fit.Lop$p[J]
      p.Lop =  rowSums(x*x) * colSums(w1t^2) + colSums(w2t^2) + colSums(w3t^2)
      p.Lop <- sqrt(p.Lop)/sum(sqrt(p.Lop))
      idx.Lop <- sample(1:n, r, T, p.Lop)
      x.Lop <- x[c(idx.Lop,pilot.ind),]; y.Lop <- y[c(idx.Lop,pilot.ind)]; 
      w.Lop = c(p.Lop[idx.Lop],rep(1/n,r0))
      fit.L = EM.alg(x.Lop, y.Lop, w.Lop, num.c = J, beta.ini = fit.Lop$beta, vari.ini = fit.Lop$vari)
   })
   ind = order(colSums(fit.L$beta))
   w = fit.L$w[,ind]; vari = fit.L$vari[ind]; beta = fit.L$beta[,ind];pr = fit.L$p[ind]
   Epsil_Lop = y.Lop - x.Lop%*%beta
   w1t = t(w * Epsil_Lop) / vari; 
   w2t.w = (t(Epsil_Lop^2) - vari)/(vari^(3/2)); 
   w2t = t(w) * w2t.w
   w3t = t(w[,-J])/pr[-J] - w[,J]/pr[J]
   A.M = M.v2(x.Lop, w.Lop, w, Epsil_Lop, w1t, w2t, w2t.w, w3t, vari = vari, pi = pr, r0+r, p, J)
   iM.Lop = solve(A.M$M) * ( n * (r+r0) )
   V_Lop = ( A.M$v %*% (t(A.M$v)/(w.Lop^2)) ) / (( n * (r+r0) )^2) 
   sigma.Lop = diag(iM.Lop %*% V_Lop %*% iM.Lop)
   return(list(beta = beta, sigma = vari, V = V_Lop, first.iter = fit.Lop$iter, iter = fit.L$iter,
               pi = w.Lop, x = x.Lop, y = y.Lop, iM = iM.Lop, w = w, p = pr , sig.est = sigma.Lop, time = time.Lop, ind = idx.Lop))
}
#######################################################################################
FULL = function(y, x, beta.ini, vari.ini)
{
   n = dim(x)[1]; p = dim(x)[2]; J = dim(beta.ini)[2]
   time.T = system.time({
      fit = EM.alg(x, y, rep(1/n,n), num.c = J, beta.ini = beta.ini, vari.ini = vari.ini)
   })
   return(list(beta = fit$beta, sigma = fit$vari, w = fit$w, p = fit$p, time = time.T, iter = fit$iter))
}
