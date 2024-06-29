load("metabric_brca_pampath.RData")
source("functions.R")

p <- dim(data)[2]
q <- dim(covariate)[2] - 1
n_all <- dim(data)[1]

K = 6

library(mvtnorm)
library(cglasso)
library(flexclust)
library(mixtools)
library(combinat)
set.seed(1234567)

## prior information matrix
prior_Sr <- array(0, dim = c(p, p, K))
for (k in 1:K) {
  prior_Sr[,,k] <- prior_omega
}

## initialization
init <- initialize_fuc(data, covariate, K=K, n.start = 100)

## the first step
lambda1.prior <- 0.2
lambda2.prior <- 0.1

res_prior <- FCGGM_prior(data, covariate, K=K, lambda1=lambda1.prior, lambda2=lambda2.prior, prior_Sr=prior_Sr,
                         a = 3, rho = 1,
                         eps = 0.05, niter = 40, maxiter=40, maxiter.AMA=10, initialize=init,
                         asymmetric=TRUE, local_appro=TRUE, penalty = "MCP")

## the second step
lambda1 <- 0.15
lambda2 <- 0.05
lambda3 <- 3.6

z <- clu_mat(res_prior$member)
res <- FCGGM_main(data, covariate, K=K, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, z = z, eta = 0.2, 
                  a = 3, rho = 1, eps = 0.05, niter = 40, maxiter=40, maxiter.AMA=10, 
                  initialize = init, 
                  average=FALSE, asymmetric=TRUE, local_appro=TRUE, 
                  penalty = "MCP", theta.fusion = TRUE)

## the sample size for each group
table(res$member)
