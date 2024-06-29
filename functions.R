################################################
####### the functions for the PriorCGGM ########
################################################

#calculate the density function values at each sample point.
## @ data: n * p matrix, Y.
## @ covariate: n * (q+1) matrix, X. 
## @ Gamma: p * (q+1) vector, the coefficient matrix.
## @ Theta: p * p matrix, the precision matrix.
f.den.vec <- function(data, covariate, Gamma, Theta){
  n = as.integer(dim(data)[1])
  p = as.integer(dim(data)[2])
  fden = NULL
  for(i in 1:n){
    fden[i] = (2*pi)^(-p/2) * abs(det(Theta))^(1/2) * exp(-1/2*t(data[i,] - as.numeric(Gamma%*%covariate[i,])) %*% Theta %*% (data[i,] - as.numeric(Gamma%*%covariate[i,])))
  }
  return(fden)
}

#Symmetrize the precision matrices using the symmetrization strategy of Cai et al. (2016).
Symmetrize <- function(X){
  p = dim(X)[1]
  for(i in 1:p){
    for(j in i:p){
      if(X[i,j] < X[j, i]){
        X[j, i] = X[i, j]
      }else{
        X[i, j] = X[j, i]
      }
    }
  }
  return(X)
}

#Calculating the derivative of the MCP
mcp_d <- function(x,lambda,a=3){
  if(lambda!=0){
    rho <- lambda*( 1 > abs(x)/( lambda*a ) )*( 1 - abs(x)/( lambda*a ))
  } else{
    rho=0
  }
  return(rho)
}
#Calculating the derivative of the SCAD
SCAD_d <- function(x,lambda,a=3.7){
  if(lambda!=0){
    rho <- lambda*(   (abs(x) <= lambda) +
                        (abs(x) > lambda) * ( lambda*a > abs(x) )*( lambda*a - abs(x) )/((a-1)*lambda)   )
  } else{
    rho=0
  }
  return(rho)
}
#Calculating the derivative of the Lasso
lasso_d <- function(x,lambda,a=3){
  if(lambda!=0){
    rho <- a/a * lambda * (abs(x) > 0)
  } else{
    rho=0
  }
  return(rho)
}

#Define the soft threshold operator
S_soft <- function(z,lambda){
  return((abs(z) - lambda)*(abs(z) - lambda > 0)*sign(z))
}
#Define the MCP threshold operator.
MCP_soft <- function(z,lambda,a=3){
  return( S_soft(z,lambda)/(1-1/a) * (abs(z) - a*lambda <= 0) + z * (abs(z) - a*lambda > 0) )
}
#Define the SCAD threshold operator.
SCAD_soft <- function(z,lambda,a=3.7){
  return( S_soft(z,lambda) * (abs(z) - lambda/2 <= 0) +
            S_soft(z,lambda)/(1-1/a) * (abs(z) - lambda/2 > 0 & abs(z) - a*lambda <= 0) +
            z * (abs(z) - a*lambda > 0) )
}
# Define the lasso threshold operator.
lasso_soft <- function(z,lambda,a=3){
  return( S_soft(z,lambda) * a/a )
}


# Generating the initial values using  nonparametric clustering method and conditional GGM.
## @ K: the given number of groups (K is larger than K0)
initialize_fuc = function(data, covariate, K, n.start = 100){
  n <- dim(data)[1]
  p <- dim(data)[2]
  q <- dim(covariate)[2] - 1
  data.whole <- cbind(data, covariate)
  kmeans.clust <- kmeans(data.whole, K, nstart = n.start)
  memb <- kmeans.clust$cluster
  prob <- kmeans.clust$size/n
  Theta <- array(0, dim = c(p, p, K))
  Gamma <- array(0, dim = c(p, q+1, K))
  S <- array(0, dim = c(p, p, K))
  for(k in 1:K)
  {
    group.data <- data[memb == k, , drop = FALSE]
    group.cov <- covariate[memb == k, , drop = FALSE]
    group.list <- datacggm(Y=group.data, X=group.cov[,-1])
    group.fit <- cglasso(. ~ ., data=group.list, lambda=0.1, rho=0.1)
    Gamma[,,k] <- t(group.fit$B[,,1, 1])
    Theta[,,k] <- group.fit$Tht[,,1, 1]
  }
  L.mat = matrix(0,n,K)
  for(jj in 1:n) L.mat[jj, memb[jj]]=1
  
  int.res <- list()
  int.res$prob <-  prob
  int.res$Gamma <-  Gamma
  int.res$Theta <- Theta
  int.res$S <- S
  int.res$memb <- memb
  int.res$L.mat <- L.mat
  return(int.res)
}

### calculate BIC value
BIC <- function(data, covariate, Gamma_hat, Theta_hat, L.mat){
  n = nrow(data)
  p = ncol(data)
  K = dim(Theta_hat)[3]
  pi_vec = apply(L.mat, 2, sum)/n
  fit.error_mat = matrix(0, n, K)
  for (k in 1:K) {
    fit.error_mat[, k] = pi_vec[k] * f.den.vec(data, covariate, Gamma_hat[, , k], Theta_hat[, , k])
  }
  fit0 = apply(fit.error_mat, 1, sum)
  fit.error = sum(log(fit0))
  fit.error = -2 * fit.error
  df = log(n) * length(which(Theta_hat != 0)) + log(n) * length(which(Gamma_hat != 0))
  bic = fit.error + df
  P = list()
  P$fit.error = fit.error
  P$df = df
  P$bic = bic
  return(P)
}



#In the first step: updating the precision matrices using the ADMM algorithm after updating the mean vectors in the EM algorithm.
## @ S: p * p * K array, the estimated pseudo sample covariance matrices of given K subgroups in the iteration of the EM algorithm.
## @ nK: K * 1 vector, the estimated sample sizes of K subgroups in the iteration of the EM algorithm.
## @ lambda1: a float value, the tuning parameter controlling the sparse of the precision matrix.
## @ prior_Sr: prior information matrix
## @ Gamma_kk_diff: L2-norm differences in the mean vector between different subgroups.
## @ a: a float value, regularization parameter in MCP, the default setting is 3.
## @ rho: a float value, the penalty parameter in ADMM algorithm of updating precision matrix Theta, the default setting is 1.
## @ tol: a float value, algorithm termination threshold.
## @ maxiter: int, Maximum number of cycles of the ADMM algorithm.
## @ maxiter.AMA: int, Maximum number of cycles of the AMA algorithm.
## @ penalty: The type of the penalty, which can be selected from c("MCP", "SCAD", "lasso").
Update_Theta_prior = function(S, nK, lambda1, prior_Sr, a = 3, rho=1,
                              maxiter=10, maxiter.AMA=5, tol=1e-2, rho.increment=1, penalty="MCP"){
  p = dim(S)[1]
  K = dim(S)[3]
  nkk = rep(1,K)
  # initialize Theta:
  Theta = array(0, dim = c(p, p, K))
  for (k in 1:K) { Theta[,,k] =  diag(1,p)}
  # initialize Xi:
  Xi = array(0, dim = c(p, p, K))
  # initialize Phi:
  Phi = array(0, dim = c(p, p, K))
  
  iter=0
  diff_value = 10
  diff_value_Xi = 10
  while( iter<maxiter && !(diff_value < tol || diff_value_Xi < tol) )
  {
    Theta.prev = Theta
    Xi.prev = Xi
    Theta=Theta.prev
    # update Theta:
    for(k.ind in 1:K){
      Sa=S[,,k.ind] - rho*Xi[,,k.ind]/nkk[k.ind] + rho*Phi[,,k.ind]/nkk[k.ind] 
      edecomp = eigen(Sa)
      D = edecomp$values
      if(is.complex(D)){Sa=Symmetrize(Sa);edecomp = eigen(Sa);D = edecomp$values}
      V = edecomp$vectors
      D2 = nkk[k.ind]/(2*rho) * ( -D + sqrt(D^2 + 4*(rho/nkk[k.ind]) ) )
      Theta[,,k.ind] = V %*% diag(D2) %*% t(V)
    }
    # update Xi:
    B = array(0, dim = c(p, p, K))
    for(k in 1:K){ B[,,k] = Theta[,,k] + Phi[,,k] }
    
    Xi_out_list = AMA_XI_prior(B,lambda1,prior_Sr,maxiter=maxiter.AMA, penalty=penalty)
    Xi = Xi_out_list$Xi
    Xi[abs(Xi) < 1e-3] <- 0
    
    
    # update the dual variable Phi:
    Phi = Phi + Theta - Xi
    iter = iter+1
    diff_value = sum(abs(Theta - Theta.prev)) / sum(abs(Theta.prev))
    diff_value_Xi = sum(abs(Xi - Xi.prev)) / (sum(abs(Xi.prev))+0.001)
    rho = rho*rho.increment
  }
  print(paste("ADMMiter=", iter))
  diff = sum(abs(Theta-Xi))
  Theta_out = list(Theta=Theta,Xi=Xi,diff=diff,iters=iter)
  return(Theta_out)
}


#In the first step: solving S-AMA algorithm, which is the key step of updating the precision matrices.
## @ B: The output of the previous step in ADMM algorithm).
## @ lambda1: a float value, the tuning parameter controlling the sparse of the precision matrix.
## @ prior_Sr: the prior information matrix
## @ a: a float value, regularization parameter in MCP, the default setting is 3.
## @ kappa: a float value, the penalty parameter in S-AMA algorithm, the default setting is 1.
## @ tol: a float value, algorithm termination threshold.
## @ maxiter: int, Maximum number of cycles of the algorithm.
## @ penalty: The type of the penalty, which can be selected from c("MCP", "SCAD", "lasso").
AMA_XI_prior = function(B, lambda1, prior_Sr, a = 3,
                        kappa=1, maxiter=5, tol=1e-2, penalty="MCP"){
  p = dim(B)[1]
  K = dim(B)[3]
  # initialize Xi:
  Xi = array(0, dim = c(p, p, K))
  for (k in 1:K) { Xi[,,k] =  diag(1,p)}
  diag_index = Xi + prior_Sr
  
  
  Z = array(0, dim = c(p, p, K))
  
  # iterations
  iter=0
  diff_value = 10
  while( iter<maxiter && diff_value > tol )
  {
    Xi.prev = Xi
    for (i in 1:p) {
      for (j in 1:p) {
        Z[i,j,] = B[i,j,]
      }
    }
    
    if(penalty=="MCP"){
      Xi = S_soft(Z,lambda1)/(1-1/a) * (abs(Z) <= a*lambda1 & diag_index==0) + Z * (abs(Z) > a*lambda1 | diag_index==1)
    }
    if(penalty=="SCAD"){
      a = 3.7
      Xi = S_soft(Z,lambda1) * (abs(Z) - 2*lambda1 <= 0  & diag_index==0) +
        S_soft(Z,lambda1*a/(a-1))/(1-1/(a-1)) * (abs(Z) - 2*lambda1 > 0 & abs(Z) - 2*a*lambda1 <= 0 & diag_index==0) +
        Z * (abs(Z) - a*lambda1 > 0 | diag_index==1)
    }
    if(penalty=="lasso"){
      Xi = S_soft(Z,lambda1) * (diag_index==0) + Z * (diag_index==1)
    }
    
    iter = iter+1
    diff_value = sum(abs(Xi - Xi.prev)^2) / sum(abs(Xi.prev)^2)
  }
  Xi_out = list(Xi=Xi,diff=diff_value,iters=iter)
  return(Xi_out)
}

## In the first step: the main function to get the estimation
## @ lambda2:  a float value, the tuning parameter controlling the sparse of the coefficient matrix.
## @ eps: a float value, algorithm termination threshold for the algorithm.
## @ niter: int, Maximum number of cycles of the EM algorithm.
## @ maxiter: int, Maximum number of cycles of the ADMM algorithm.
## @ maxiter.AMA: int, Maximum number of cycles of the AMA algorithm.
## @ initialize: the inital value.
## @ asymmetric: whether need symmetrize the precision matrix, default is TRUE.
## @ local_appro: whether use local_appro to estimate coefficient matrix, default is TRUE.
FCGGM_prior <- function(data, covariate, K, lambda1, lambda2, prior_Sr, a = 3, rho = 1,
                        eps = 5e-2, niter = 20, maxiter=10, maxiter.AMA=5, initialize=init,
                        asymmetric=TRUE, local_appro=TRUE, penalty = "MCP"){
  
  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])
  q <- as.integer(dim(covariate)[2])-1
  K_c <- as.matrix(combn(K,2))
  
  #initialization
  Theta = initialize$Theta
  Gamma = initialize$Gamma
  prob = initialize$prob
  L.mat = initialize$L.mat
  memb = apply(L.mat,1,which.max)
  
  
  # EM algorithm
  f.mat = matrix(0,n,K)
  L.mat.old = matrix(0,n,K)
  Gamma.old = array(0, dim = c(p, q+1, K))
  Theta.old = array(0, dim = c(p, p, K))
  Cyy = array(0, dim = c(p, p, K))
  Cyx = array(0, dim = c(p, q+1, K))
  Cxx = array(0, dim = c(q+1, q+1, K))
  
  t = 0
  diff_gamma = 10
  diff_theta = 10
  L_list = list()
  divide_ind = likel.before = likel.after = NULL
  
  while( !(diff_theta < eps || diff_gamma < eps) && t < niter ){
    prob.old = prob
    Gamma.old = Gamma
    Theta.old = Theta
    L.mat.old = L.mat
    
    K_c <- as.matrix(combn(K,2))
    f.mat = matrix(0,n,K)
    Cyy = array(0, dim = c(p, p, K))
    Cyx = array(0, dim = c(p, q+1, K))
    Cxx = array(0, dim = c(q+1, q+1, K))
    
    for(k.ind in 1:K) {
      f.mat[,k.ind]=f.den.vec(data, covariate, Gamma.old[,,k.ind], Theta.old[,,k.ind])
    }
    
    # update L and pi
    for(k.ind in 1:K) {
      for(i in 1:n) {
        L.mat[i,k.ind] = prob.old[k.ind] * f.mat[i,k.ind] / prob.old %*% f.mat[i,]
      }
      prob[k.ind] = mean(L.mat[,k.ind])
      if(sum(L.mat[,k.ind])!=0){
        Cyy[,,k.ind] = (t(data) %*% diag(L.mat[,k.ind]) %*% data)/sum(L.mat[,k.ind])
        Cyx[,,k.ind] = (t(data) %*% diag(L.mat[,k.ind]) %*% covariate)/sum(L.mat[,k.ind])
        Cxx[,,k.ind] = (t(covariate) %*% diag(L.mat[,k.ind]) %*% covariate)/sum(L.mat[,k.ind])
      }
      else{
        Cyy[,,k.ind] = matrix(0, nrow = p, ncol = p)
        Cyx[,,k.ind] = matrix(0, nrow = p, ncol = q+1)
        Cxx[,,k.ind] = matrix(0, nrow= q+1, ncol = q+1)
      }
      
    }
    nK = apply(L.mat,2,sum)
    
    Theta_kk_diff <- rep(0,dim(K_c)[2])
    for (l in 1:dim(K_c)[2]) {
      Theta_kk_diff[l] <- sum((Theta.old[,,K_c[1,l]] - Theta.old[,,K_c[2,l]])^2)
    }
    
    Gamma_kk_diff <- rep(0,dim(K_c)[2])
    for (l in 1:dim(K_c)[2]) {
      Gamma_kk_diff[l] <- sum((Gamma.old[,,K_c[1,l]] - Gamma.old[,,K_c[2,l]])^2)
    }
    
    # update coefficient matrices Gamma
    for(k.ind in 1:K){
      for(i in 1:p){
        for(j in 1:(q+1)){
          e_j = rep(0, q+1); e_j[j] = 1
          e_i = rep(0, p); e_i[i]  = 1
          tmp = t(e_j) %*% Cxx[,,k.ind] %*% t(Gamma[,,k.ind]) %*% Theta.old[,,k.ind] %*% e_i - (Theta.old[i,i,k.ind]*Cxx[j,j,k.ind]*Gamma.old[i,j,k.ind]) - t(e_j) %*% t(Cyx[,,k.ind]) %*% Theta.old[,,k.ind] %*% e_i
          hij = -nK[k.ind] * tmp 
          if(penalty=="MCP"){
            mcp_lambda1 = mcp_d(abs(Gamma[i,j,k.ind]), lambda2, a)
          }
          if(penalty=="SCAD"){
            v_k = sum(SCAD_d(tau_k, lambda3, a=3.7) / (tau_k+0.00001) * Gamma[i,j,-k.ind])
            mcp_lambda1 = SCAD_d(abs(Gamma[i,j,k.ind]), lambda2, a=3.7)
            v_k_hat = sum(SCAD_d(tau_k, lambda3, a=3.7) / (tau_k+0.00001))
          }
          if(penalty=="lasso"){
            v_k = sum(lasso_d(tau_k, lambda3, a) / (tau_k+0.00001) * Gamma[i,j,-k.ind])
            mcp_lambda1 = lasso_d(abs(Gamma[i,j,k.ind]), lambda2, a)
            v_k_hat = sum(lasso_d(tau_k, lambda3, a) / (tau_k+0.00001))
          }
          
          if(local_appro){
            Gamma[i,j,k.ind] = (hij) / (nK[k.ind]*Theta.old[i,i,k.ind]*Cxx[j,j,k.ind] + n*mcp_lambda1/(abs(Gamma[i,j,k.ind])+0.00001)+0.00001)
          } else {
            if(n*mcp_lambda1 >= abs(hij + n*v_k)){
              Gamma[i,j,k.ind] = 0
            }else{
              Gamma[i,j,k.ind] = (hij + n*v_k - n*mcp_lambda1*sign(Gamma[i,j,k.ind])) / (nK[k.ind]*Theta.old[i,i,k.ind]*Cxx[j,j,k.ind] + n*v_k_hat)
            }
          }
          
        }
      }
    }
    Gamma[abs(Gamma) < 1e-3] <- 0
    
    # update precision matrices Theta
    Gamma_kk_diff <- rep(0,dim(K_c)[2])
    for (l in 1:dim(K_c)[2]) {
      Gamma_kk_diff[l] <- sum((Gamma[,,K_c[1,l]] - Gamma[,,K_c[2,l]])^2)
    }
    # define the pseudo sample covariance matrices S
    S = array(0, dim = c(p, p, K))
    for (k.ind in 1:K) {
      S[,,k.ind] = Cyy[,,k.ind] - Cyx[,,k.ind] %*% t(Gamma[,,k.ind]) - Gamma[,,k.ind] %*% t(Cyx[,,k.ind]) + Gamma[,,k.ind] %*% Cxx[,,k.ind] %*% t(Gamma[,,k.ind])
    }
    Theta_out = Update_Theta_prior(S,nK,lambda1, prior_Sr, tol=eps, maxiter=maxiter, maxiter.AMA=maxiter.AMA, penalty=penalty)
    Theta = Theta_out$Theta
    Xi = Theta_out$Xi
    
    Theta[abs(Theta) < 1e-3] <- 0
    
    diff_gamma = norm(Gamma.old-Gamma,type="2")/(norm(Gamma,type="2")+0.001)
    diff_theta = norm(Theta.old-Theta,type="2")/(norm(Theta,type="2")+0.001)
    
    L_list[[t+1]] = L.mat
    t = t + 1
    
  }
  print(paste("EMiter=", t))
  K_0 = K
  Theta_final <- Theta; Gamma_final <- Gamma; Xi_final <- Xi; prob0 <- prob; L.mat0 <- L.mat
  if(asymmetric){
    for(k in 1:K_0) {
      Theta_final[,,k] = Symmetrize(Theta_final[,,k])
      Xi_final[,,k] = Symmetrize(Xi_final[,,k])
    }
  }
  
  Theta_final[abs(Theta_final) < 1e-2] <- 0
  Gamma_final[abs(Gamma_final) < 1e-2] <- 0
  member = apply(L.mat0,1,function(a){which(a == max(a))[1]})
  
  FGGM_res <- list();FGGM_res$Gamma <-  Gamma_final;FGGM_res$Theta <- Theta_final
  FGGM_res$Xi <- Xi_final; FGGM_res$niter <- t; 
  FGGM_res$prob0 <- prob0; FGGM_res$L.mat0 <- L.mat0;
  FGGM_res$Theta0 <- Theta; FGGM_res$member <- member;
  FGGM_res$Xi0 <- Xi;FGGM_res$Gamma0 <- Gamma;
  FGGM_res$Gamma_kk_diff <- Gamma_kk_diff; FGGM_res$Theta_kk_diff <- Theta_kk_diff
  return(FGGM_res)
}

############# obtain membership matrix z ##############
clu_mat <- function(Y_prior){
  clu_matrix <- matrix(0, nrow = length(Y_prior), ncol = max(Y_prior))
  for (i in 1:length(Y_prior)) {
    clu_matrix[i, Y_prior[i]] <- 1
  }
  return(clu_matrix)
}


#In the second step: updating the precision matrices using the ADMM algorithm after updating the coefficient matrices in the EM algorithm.
## @ S: p * p * K array, the estimated pseudo sample covariance matrices of given K subgroups in the iteration of the EM algorithm.
## @ nK: K * 1 vector, the estimated sample sizes of K subgroups in the iteration of the EM algorithm.
## @ lambda1: a float value, the tuning parameter controlling the sparse of the precision matrix.
## @ lambda3: a float value, the tuning parameter controlling the number of subgroup.
## @ Gamma_kk_diff: L2-norm differences in the coefficient matrices between different subgroups.
## @ K_c: All combinations of natural numbers within K.
## @ a: a float value, regularization parameter in MCP, the default setting is 3.
## @ rho: a float value, the penalty parameter in ADMM algorithm of updating precision matrix Theta, the default setting is 1.
## @ tol: a float value, algorithm termination threshold.
## @ maxiter: int, Maximum number of cycles of the ADMM algorithm.
## @ maxiter.AMA: int, Maximum number of cycles of the AMA algorithm.
## @ penalty: The type of the penalty, which can be selected from c("MCP", "SCAD", "lasso").
## @ theta.fusion: Whether or not the fusion penalty term contains elements of the precision matrices.
Update_Theta = function(S, nK, lambda1, lambda3, Gamma_kk_diff, K_c, a = 3, rho=1,
                        maxiter=10, maxiter.AMA=5, tol=1e-2, rho.increment=1, penalty="MCP", theta.fusion=TRUE){
  p = dim(S)[1]
  K = dim(S)[3]
  nkk = rep(1,K)
  # initialize Theta:
  Theta = array(0, dim = c(p, p, K))
  for (k in 1:K) { Theta[,,k] =  diag(1,p)}
  # initialize Xi:
  Xi = array(0, dim = c(p, p, K))
  # initialize Phi:
  Phi = array(0, dim = c(p, p, K))
  
  iter=0
  diff_value = 10
  diff_value_Xi = 10
  while( iter<maxiter && !(diff_value < tol || diff_value_Xi < tol) )
  {
    Theta.prev = Theta
    Xi.prev = Xi
    Theta=Theta.prev
    # update Theta:
    for(k.ind in 1:K){
      Sa=S[,,k.ind] - rho*Xi[,,k.ind]/nkk[k.ind] + rho*Phi[,,k.ind]/nkk[k.ind] 
      edecomp = eigen(Sa)
      D = edecomp$values
      if(is.complex(D)){Sa=Symmetrize(Sa);edecomp = eigen(Sa);D = edecomp$values}
      V = edecomp$vectors
      D2 = nkk[k.ind]/(2*rho) * ( -D + sqrt(D^2 + 4*(rho/nkk[k.ind]) ) )
      Theta[,,k.ind] = V %*% diag(D2) %*% t(V)
    }
    # update Xi:
    B = array(0, dim = c(p, p, K))
    for(k in 1:K){ B[,,k] = Theta[,,k] + Phi[,,k] }
    
    Xi_out_list = AMA_XI(B,K_c,lambda1,lambda3,Gamma_kk_diff,maxiter=maxiter.AMA, penalty=penalty, theta.fusion=theta.fusion)
    Xi = Xi_out_list$Xi
    Xi[abs(Xi) < 1e-3] <- 0
    V = Xi_out_list$V
    V_kk = round(apply(V^2,3,sum),4)
    
    # update the dual variable Phi:
    Phi = Phi + Theta - Xi
    iter = iter+1
    diff_value = sum(abs(Theta - Theta.prev)) / sum(abs(Theta.prev))
    diff_value_Xi = sum(abs(Xi - Xi.prev)) / (sum(abs(Xi.prev))+0.001)
    rho = rho*rho.increment
  }
  print(paste("ADMMiter=", iter))
  diff = sum(abs(Theta-Xi))
  Theta_out = list(Theta=Theta,Xi=Xi,V_kk=V_kk,diff=diff,iters=iter)
  return(Theta_out)
}

#In the second step: Solving S-AMA algorithm, which is the key step of updating the precision matrices.
## @ B: The output of the previous step in ADMM algorithm).
## @ K_c: All combinations of natural numbers within K.
## @ lambda1: a float value, the tuning parameter controlling the sparse of the precision matrix.
## @ lambda3: a float value, the tuning parameter controlling the number of subgroup.
## @ Gamma_kk_diff: L2-norm differences in the mean vector between different subgroups.
## @ a: a float value, regularization parameter in MCP, the default setting is 3.
## @ kappa: a float value, the penalty parameter in S-AMA algorithm, the default setting is 1.
## @ tol: a float value, algorithm termination threshold.
## @ maxiter: int, Maximum number of cycles of the algorithm.
## @ penalty: The type of the penalty, which can be selected from c("MCP", "SCAD", "lasso").
## @ theta.fusion: Whether or not the fusion penalty term contains elements of the precision matrices.
AMA_XI = function(B, K_c, lambda1, lambda3, Gamma_kk_diff, a = 3,
                  kappa=1, maxiter=5, tol=1e-2, penalty="MCP", theta.fusion=TRUE){
  p = dim(B)[1]
  K = dim(B)[3]
  num_Kc = dim(K_c)[2]
  # initialize Xi:
  Xi = array(0, dim = c(p, p, K))
  for (k in 1:K) { Xi[,,k] =  diag(1,p)}
  diag_index = Xi
  # initialize V:
  V = array(0, dim = c(p, p, num_Kc))
  # initialize Delta:
  Delta = array(0, dim = c(p, p, num_Kc))
  
  Z = array(0, dim = c(p, p, K))
  Omega = array(0, dim = c(p, p, num_Kc))
  e_k12 = matrix(0,K,num_Kc)
  for (i in 1:K) {e_k12[i,which(K_c[1,] == i)] = 1;e_k12[i,which(K_c[2,] == i)] = -1}
  
  # iterations
  iter=0
  diff_value = 10
  while( iter<maxiter && diff_value > tol )
  {
    V.prev = V
    Xi.prev = Xi
    Delta.prev = Delta
    for (i in 1:p) {
      for (j in 1:p) {
        Z[i,j,] = B[i,j,] + apply(t(Delta[i,j,] * t(e_k12)),1,sum)
      }
    }
    
    if(penalty=="MCP"){
      Xi = S_soft(Z,lambda1)/(1-1/a) * (abs(Z) <= a*lambda1 & diag_index==0) + Z * (abs(Z) > a*lambda1 | diag_index==1)
    }
    if(penalty=="SCAD"){
      a = 3.7
      Xi = S_soft(Z,lambda1) * (abs(Z) - 2*lambda1 <= 0  & diag_index==0) +
        S_soft(Z,lambda1*a/(a-1))/(1-1/(a-1)) * (abs(Z) - 2*lambda1 > 0 & abs(Z) - 2*a*lambda1 <= 0 & diag_index==0) +
        Z * (abs(Z) - a*lambda1 > 0 | diag_index==1)
    }
    if(penalty=="lasso"){
      Xi = S_soft(Z,lambda1) * (diag_index==0) + Z * (diag_index==1)
    }
    
    if(!theta.fusion){lambda3 = 0}
    for (l in 1:num_Kc) {
      Omega[,,l] = Xi[,,K_c[1,l]] - Xi[,,K_c[2,l]] - Delta[,,l]/kappa
      Omegal2 = sum(Omega[,,l]^2)
      if(penalty=="MCP"){
        V[,,l] = Omega[,,l] * MCP_soft(sqrt(Omegal2 + Gamma_kk_diff[l]),lambda3/kappa) / (sqrt(Omegal2+Gamma_kk_diff[l])+1e-8)
      }
      if(penalty=="SCAD"){
        V[,,l] = Omega[,,l] * SCAD_soft(sqrt(Omegal2 + Gamma_kk_diff[l]),lambda3/kappa) / (sqrt(Omegal2+Gamma_kk_diff[l])+1e-8)
      }
      if(penalty=="lasso"){
        V[,,l] = Omega[,,l] * lasso_soft(sqrt(Omegal2 + Gamma_kk_diff[l]),lambda3/kappa) / (sqrt(Omegal2+Gamma_kk_diff[l])+1e-8)
      }
      Delta[,,l] = Delta[,,l] + kappa * ( V[,,l] - Xi[,,K_c[1,l]] + Xi[,,K_c[2,l]] )* (sum(V[,,l]^2) > 0 )
    }
    
    iter = iter+1
    diff_value = sum(abs(Xi - Xi.prev)^2) / sum(abs(Xi.prev)^2)
  }
  Xi_out = list(Xi=Xi,V=V,Delta=Delta,diff=diff_value,iters=iter)
  return(Xi_out)
}

#Truncation of small differences in parameter estimates between subgroups.
cut_diff_ama <- function(V_kk,K_c,K,cutoff=0.01){
  V_kk_num = which(V_kk < cutoff)
  K_group_final = list()
  if(length(V_kk_num) > 0){
    K_group = list()
    for (j in 1:length(V_kk_num)) {
      K_group[[j]] = K_c[,V_kk_num[j]]
    }
    outnum = setdiff(1:K,Reduce(union,K_group))
    if(length(outnum) > 0){
      for (j in 1:length(outnum)) {
        K_group[[length(V_kk_num)+j]] = outnum[j]
      }
      K_group[[length(V_kk_num)+j+1]] = K
    } else{
      K_group[[length(V_kk_num)+1]] = K
    }
    
    kk = 1
    repeat{
      repeat{
        K_group_old=K_group
        k_del = NULL
        for (kkk in setdiff(1:length(K_group),1) ) {
          if(length(Reduce(intersect,list(K_group[[1]],K_group_old[[kkk]]))) > 0){
            K_group[[1]] = sort(unique(c(K_group[[1]],K_group_old[[kkk]])))
            k_del = c(k_del,kkk)
          }
        }
        if(length(k_del) > 0){
          for (j in sort(k_del,decreasing = TRUE)) {
            K_group[[j]] = NULL
          }
        }
        if(length(K_group_old) == length(K_group)){break}
      }
      K_group_final[[kk]] = K_group[[1]]
      if(kk==1 && length(K_group) == 1){
        # print("Warning: Only one cluster!")s
        break}
      if(length(K_group) == 2){K_group_final[[kk+1]] = K_group[[2]];break}
      if(kk>1 && length(K_group) == 1){break}
      K_group[[1]] = NULL
      kk = kk+1
    }
  }else {
    for (k in 1:K) {
      K_group_final[[k]] = k
    }
  }
  return(K_group_final)
}

## In the second step: the main function to get the estimation
## @ z: the membership matrix obtained from the first step
## @ eta: the weight parameter in the second step
FCGGM_main <- function(data, covariate, K, lambda1, lambda2, lambda3, z, 
                        eta, a = 3, rho = 1, eps = 0.05, niter = 20, maxiter = 10, 
                        maxiter.AMA = 5, initialize = init, 
                        average = FALSE, asymmetric = TRUE, local_appro = TRUE, 
                        penalty = "MCP", theta.fusion = TRUE){
  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])
  q <- as.integer(dim(covariate)[2]) - 1
  K_c <- as.matrix(combn(K, 2))
  Theta = initialize$Theta
  Gamma = initialize$Gamma
  prob = initialize$prob
  L.mat = initialize$L.mat
  memb = apply(L.mat, 1, which.max)
  eta.mat = z
  f.mat = matrix(0, n, K)
  L.mat.old = matrix(0, n, K)
  Gamma.old = array(0, dim = c(p, q + 1, K))
  Theta.old = array(0, dim = c(p, p, K))
  Cyy = array(0, dim = c(p, p, K))
  Cyx = array(0, dim = c(p, q + 1, K))
  Cxx = array(0, dim = c(q + 1, q + 1, K))
  t = 0
  diff_gamma = 10
  diff_theta = 10
  L_list = list()
  divide_ind = likel.before = likel.after = NULL
  while (!(diff_theta < eps || diff_gamma < eps) && t < niter){
    prob.old = prob
    Gamma.old = Gamma
    Theta.old = Theta
    L.mat.old = L.mat
    K_c <- as.matrix(combn(K, 2))
    f.mat = eta.mat = matrix(0, n, K)
    Cyy = array(0, dim = c(p, p, K))
    Cyx = array(0, dim = c(p, q + 1, K))
    Cxx = array(0, dim = c(q + 1, q + 1, K))
    for (k.ind in 1:K){
      f.mat[, k.ind] = f.den.vec(data, covariate, Gamma.old[, , k.ind], Theta.old[, , k.ind])
    }
    for (k.ind in 1:K){
      for (i in 1:n){
        L.mat[i, k.ind] = prob.old[k.ind] * f.mat[i, k.ind]/prob.old %*% f.mat[i, ]
      }
      prob[k.ind] = mean(L.mat[, k.ind])
      eta.mat[, k.ind] = (1 - eta) * L.mat[, k.ind] + eta * z[, k.ind]
      if (sum(eta.mat[, k.ind]) != 0) {
        Cyy[, , k.ind] = (t(data) %*% diag(eta.mat[, k.ind]) %*% data)/sum(eta.mat[, k.ind])
        Cyx[, , k.ind] = (t(data) %*% diag(eta.mat[, k.ind]) %*% covariate)/sum(eta.mat[, k.ind])
        Cxx[, , k.ind] = (t(covariate) %*% diag(eta.mat[, k.ind]) %*% covariate)/sum(eta.mat[, k.ind])
      }
      else {
        Cyy[, , k.ind] = matrix(0, nrow = p, ncol = p)
        Cyx[, , k.ind] = matrix(0, nrow = p, ncol = q + 1)
        Cxx[, , k.ind] = matrix(0, nrow = q + 1, ncol = q + 1)
      }
    }
    nK = apply(eta.mat, 2, sum)
    Theta_kk_diff <- rep(0, dim(K_c)[2])
    for (l in 1:dim(K_c)[2]){
      Theta_kk_diff[l] <- sum((Theta.old[, , K_c[1, l]] - Theta.old[, , K_c[2, l]])^2)
      if (!theta.fusion){
        Theta_kk_diff[l] = 0
      }
    }
    Gamma_kk_diff <- rep(0, dim(K_c)[2])
    for (l in 1:dim(K_c)[2]){
      Gamma_kk_diff[l] <- sum((Gamma.old[, , K_c[1, l]] - Gamma.old[, , K_c[2, l]])^2)
    }
    for (k.ind in 1:K){
      for (i in 1:p){
        for (j in 1:(q + 1)){
          e_j = rep(0, q + 1)
          e_j[j] = 1
          e_i = rep(0, p)
          e_i[i] = 1
          tmp = t(e_j) %*% Cxx[, , k.ind] %*% t(Gamma[, , k.ind]) %*% Theta.old[, , k.ind] %*% e_i - (Theta.old[i, i, k.ind] * Cxx[j, j, k.ind] * Gamma.old[i, j, k.ind]) - t(e_j) %*% t(Cyx[, , k.ind]) %*% Theta.old[, , k.ind] %*% e_i
          hij = -nK[k.ind] * tmp
          tau_k = sqrt(Gamma_kk_diff[ceiling(which(K_c == k.ind)/2)] + Theta_kk_diff[ceiling(which(K_c == k.ind)/2)])
          if (penalty == "MCP"){
            v_k = sum(mcp_d(tau_k, lambda3, a)/(tau_k + 1e-05) * Gamma[i, j, -k.ind])
            mcp_lambda1 = mcp_d(abs(Gamma[i, j, k.ind]), lambda2, a)
            v_k_hat = sum(mcp_d(tau_k, lambda3, a)/(tau_k + 1e-05))
          }
          if (penalty == "SCAD"){
            v_k = sum(SCAD_d(tau_k, lambda3, a = 3.7)/(tau_k + 1e-05) * Gamma[i, j, -k.ind])
            mcp_lambda1 = SCAD_d(abs(Gamma[i, j, k.ind]), lambda2, a = 3.7)
            v_k_hat = sum(SCAD_d(tau_k, lambda3, a = 3.7)/(tau_k + 1e-05))
          }
          if (penalty == "lasso"){
            v_k = sum(lasso_d(tau_k, lambda3, a)/(tau_k + 1e-05) * Gamma[i, j, -k.ind])
            mcp_lambda1 = lasso_d(abs(Gamma[i, j, k.ind]), lambda2, a)
            v_k_hat = sum(lasso_d(tau_k, lambda3, a)/(tau_k + 1e-05))
          }
          if (local_appro){
            Gamma[i, j, k.ind] = (hij + n * v_k)/(nK[k.ind] * Theta.old[i, i, k.ind] * Cxx[j, j, k.ind] + n * v_k_hat + n * mcp_lambda1/(abs(Gamma[i,j, k.ind]) + 1e-05) + 1e-05)
          }
          else{
            if(n * mcp_lambda1 >= abs(hij + n * v_k)){
              Gamma[i, j, k.ind] = 0
            }
            else{
              Gamma[i, j, k.ind] = (hij + n * v_k - n * mcp_lambda1 * sign(Gamma[i, j, k.ind]))/(nK[k.ind] * Theta.old[i, i, k.ind] * Cxx[j, j, k.ind] + n * v_k_hat)
            }
          }
        }
      }
    }
    Gamma[abs(Gamma) < 0.001] <- 0
    Gamma_kk_diff <- rep(0, dim(K_c)[2])
    for (l in 1:dim(K_c)[2]){
      Gamma_kk_diff[l] <- sum((Gamma[, , K_c[1, l]] - Gamma[, , K_c[2, l]])^2)
    }
    S = array(0, dim = c(p, p, K))
    for (k.ind in 1:K){
      S[, , k.ind] = Cyy[, , k.ind] - Cyx[, , k.ind] %*% 
        t(Gamma[, , k.ind]) - Gamma[, , k.ind] %*% t(Cyx[, , k.ind]) + Gamma[, , k.ind] %*% Cxx[, , k.ind] %*% t(Gamma[, , k.ind])
    }
    Theta_out = Update_Theta(S, nK, lambda1, lambda3, Gamma_kk_diff, 
                             K_c, tol = eps, maxiter = maxiter, maxiter.AMA = maxiter.AMA, 
                             penalty = penalty, theta.fusion = theta.fusion)
    Theta = Theta_out$Theta
    Xi = Theta_out$Xi
    V_kk = Theta_out$V_kk
    Theta[abs(Theta) < 0.001] <- 0
    diff_gamma = norm(Gamma.old - Gamma, type = "2")/(norm(Gamma, type = "2") + 0.001)
    diff_theta = norm(Theta.old - Theta, type = "2")/(norm(Theta, type = "2") + 0.001)
    likel.before[t + 1] <- BIC(data, covariate, Gamma, Theta, eta.mat)$bic
    likel.after[t + 1] <- loglikel.after(data, covariate, eta.mat, Theta, Gamma)
    if (diff_theta > eps && diff_gamma > eps){
      group_final = comb.group(L.mat, K, likel.before[t + 1], likel.after[t + 1])
      K = length(group_final)
      Gamma_final = array(0, dim = c(p, q + 1, K))
      Theta_final = array(0, dim = c(p, p, K))
      Xi_final = array(0, dim = c(p, p, K))
      prob0 = rep(0, K)
      L.mat0 = matrix(0, n, K)
      z0 = matrix(0, n, K)
      for (l in 1:K){
        gg = group_final[[l]]
        prob0[l] = sum(prob[gg])
        if (length(gg) > 1){
          L.mat0[, l] = apply(L.mat[, gg], 1, sum)
          z0[, l] = apply(z[, gg], 1, sum)
          Gamma_final[, , l] = Gamma[, , gg[which.max(nK[gg])]]
          Theta_final[, , l] = Theta[, , gg[which.max(nK[gg])]]
          Xi_final[, , l] = Xi[, , gg[which.max(nK[gg])]]
        }
        else{
          Gamma_final[, , l] = Gamma[, , gg]
          Theta_final[, , l] = Theta[, , gg]
          Xi_final[, , l] = Xi[, , gg]
          L.mat0[, l] = L.mat[, gg]
          z0[, l] = z[, gg]
        }
      }
      if (asymmetric){
        for (k in 1:K){
          Theta_final[, , k] = Symmetrize(Theta_final[, , k])
          Xi_final[, , k] = Symmetrize(Xi_final[, , k])
        }
      }
      prob = prob0
      L.mat = L.mat0
      Theta = Theta_final
      Gamma = Gamma_final
      z = z0
    }
    L_list[[t + 1]] = L.mat
    t = t + 1
  }
  print(paste("EMiter=", t))
  group_final = cut_diff_ama(V_kk, K_c, K, cutoff = 0.01)
  K_0 = length(group_final)
  Gamma_final = array(0, dim = c(p, q + 1, K_0))
  Theta_final = array(0, dim = c(p, p, K_0))
  Xi_final = array(0, dim = c(p, p, K_0))
  prob0 = rep(0, K_0)
  L.mat0 = matrix(0, n, K_0)
  for (l in 1:K_0){
    gg = group_final[[l]]
    prob0[l] = sum(prob[gg])
    if (length(gg) > 1){
      L.mat0[, l] = apply(L.mat[, gg], 1, sum)
      if (!average){
        Gamma_final[, , l] = Gamma[, , gg[which.max(nK[gg])]]
        Theta_final[, , l] = Theta[, , gg[which.max(nK[gg])]]
        Xi_final[, , l] = Xi[, , gg[which.max(nK[gg])]]
      }
      else {
        Gamma_final[, , l] = Gamma[, , gg[1]]/length(gg)
        Theta_final[, , l] = Theta[, , gg[1]]/length(gg)
        Xi_final[, , l] = Xi[, , gg[1]]/length(gg)
        for (gi in gg[-1]) {
          Gamma_final[, , l] = Gamma_final[, , l] + Gamma[, , gi]/length(gg)
          Theta_final[, , l] = Theta_final[, , l] + Theta[, , gi]/length(gg)
          Xi_final[, , l] = Xi_final[, , l] + Xi[, , gi]/length(gg)
        }
      }
    }
    else{
      Gamma_final[, , l] = Gamma[, , gg]
      Theta_final[, , l] = Theta[, , gg]
      Xi_final[, , l] = Xi[, , gg]
      L.mat0[, l] = L.mat[, gg]
    }
  }
  if (asymmetric) {
    for (k in 1:K_0) {
      Theta_final[, , k] = Symmetrize(Theta_final[, , k])
      Xi_final[, , k] = Symmetrize(Xi_final[, , k])
    }
  }
  Theta_final[abs(Theta_final) < 0.01] <- 0
  Gamma_final[abs(Gamma_final) < 0.01] <- 0
  member = apply(L.mat0, 1, function(a) {which(a == max(a))[1]})
  FGGM_res <- list()
  FGGM_res$Gamma <- Gamma_final
  FGGM_res$Theta <- Theta_final
  FGGM_res$Xi <- Xi_final
  FGGM_res$niter <- t
  FGGM_res$diff_Xi <- V_kk
  FGGM_res$prob0 <- prob0
  FGGM_res$L.mat0 <- L.mat0
  FGGM_res$Theta0 <- Theta
  FGGM_res$member <- member
  FGGM_res$Xi0 <- Xi
  FGGM_res$Gamma0 <- Gamma
  FGGM_res$group <- group_final
  FGGM_res$Gamma_kk_diff <- Gamma_kk_diff
  FGGM_res$Theta_kk_diff <- Theta_kk_diff
  FGGM_res$likel.before <- likel.before
  FGGM_res$likel.after <- likel.after
  return(FGGM_res)
}

#determine whether or not need to adjustment
#like1: BIC before adjustment
#like2: BIC after adjustment
#Only re-assign the smallest group if BIC getting smaller and K > 2.
comb.group <- function(L.mat, K, like1, like2){
  group_list <- list()
  mem.now <- apply(L.mat,1,function(a){which(a == max(a))[1]})
  if(length(table(mem.now)) != K){
    which.miss <- setdiff(1:K, unique(mem.now))
    for(gr in 1:length(setdiff(1:K, which.miss))){
      group_list[[gr]] <- setdiff(1:K, which.miss)[gr]
    }
  }else{
    min.size <- min(as.vector(table(mem.now)))
    if(like2 < like1 & K > 2 & is.na(like1)==F & is.na(like2)==F & min.size >= n_all/5 ){
      ind.now <- which.min(as.numeric(table(mem.now)))
      if(length(which(mem.now==ind.now))!=1){
        second.group <- apply(L.mat[which(mem.now==ind.now),], 1, function(a) order(a, decreasing = T))[2,]
      }else{
        second.group <- order(L.mat[which(mem.now==ind.now),], decreasing = T)[2]
      }
      mem.now[which(mem.now==ind.now)] <- second.group
      group_comb <- c(ind.now, as.numeric(names(table(second.group))[which.max(table(second.group))]))
      for(gr in 1:length(setdiff(1:K, group_comb))){
        group_list[[gr]] <- setdiff(1:K, group_comb)[gr]
      }
      group_list[[gr+1]] <- group_comb
    }
    else{
      for (k in 1:K) {
        group_list[[k]] = k
      }
    }
  }
  return(group_list)
}

#BIC before adjustment
loglikel.before <- function(data, covariate, L.mat, Theta, Gamma){
  K.est <- dim(Theta)[3]
  n <- dim(data)[1]
  pi_vec = apply(L.mat, 2, sum)/n
  f.mat.est <- matrix(NA, nrow = n, ncol = K.est)
  for (k.ind in 1:K.est) {
    f.mat.est[, k.ind] = f.den.vec(data, covariate, Gamma[, , k.ind], Theta[, , k.ind])
  }
  log.each <- NULL
  for (i in 1:n) {
    log.each[i] <- log(t(as.vector(pi_vec)) %*% as.vector(f.mat.est[i,]))
  }
  bic.value <- -sum(log.each) + sum(Theta != 0) + sum(Gamma != 0)
  return(bic.value)
}

#BIC after adjustment
loglikel.after <- function(data, covariate, L.mat, Theta, Gamma){
  n <- dim(data)[1]
  p <- dim(data)[2]
  q <- dim(covariate)[2] - 1
  nK = apply(L.mat, 2, sum)
  prob = apply(L.mat, 2, sum)/n
  mem.now <- apply(L.mat, 1, function(a) {which(a == max(a))[1]})
  if (dim(Theta)[3] > 2 & length(table(mem.now)) == dim(Theta)[3]){
    group_list <- list()
    ind.now <- which.min(as.numeric(table(mem.now)))
    if (length(which(mem.now == ind.now)) != 1){
      second.group <- apply(L.mat[which(mem.now == ind.now), ], 1, function(a) order(a, decreasing = T))[2,]
    }
    else{
      second.group <- order(L.mat[which(mem.now == ind.now),], decreasing = T)[2]
    }
    mem.now[which(mem.now == ind.now)] <- second.group
    group_comb <- c(ind.now, as.numeric(names(table(second.group))[which.max(table(second.group))]))
    for (gr in 1:length(setdiff(1:dim(Theta)[3], group_comb))) {
      group_list[[gr]] <- setdiff(1:dim(Theta)[3], group_comb)[gr]
    }
    group_list[[gr + 1]] <- group_comb
    K.t = length(group_list)
    Gamma_final = array(0, dim = c(p, q + 1, K.t))
    Theta_final = array(0, dim = c(p, p, K.t))
    prob0 = rep(0, K.t)
    L.mat0 = matrix(0, n, K.t)
    for (l in 1:K.t) {
      gg = group_list[[l]]
      prob0[l] = sum(prob[gg])
      if (length(gg) > 1) {
        L.mat0[, l] = apply(L.mat[, gg], 1, sum)
        Gamma_final[, , l] = Gamma[, , gg[which.max(nK[gg])]]
        Theta_final[, , l] = Theta[, , gg[which.max(nK[gg])]]
      }
      else {
        Gamma_final[, , l] = Gamma[, , gg]
        Theta_final[, , l] = Theta[, , gg]
        L.mat0[, l] = L.mat[, gg]
      }
    }
    for (k in 1:K.t) {
      Theta_final[, , k] = Symmetrize(Theta_final[, , k])
    }
    bic.value <- BIC(data, covariate, Gamma_final, Theta_final, L.mat0)$bic
  }
  else {
    bic.value <- Inf
  }
  return(bic.value)
}
