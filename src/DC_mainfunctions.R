library(igraph)
library(ellipse)
library(CovTools)
library(spherepc)
library(abind)
library(fdapace)
library(colorspace)
library(plotrix)
library(expm)
library(vegan)
library(readxl)
library(fdadensity)
library(ggplot2)
library(fields)
library(scatterplot3d)
library(jsonlite)
library(dslabs)
library(frechet)


## Obtain ambient distance from data frame 
ambient.dist = function(df){
  
  # df : (k x p) data frame on R^p
  
  final.dist = as.matrix(dist(df,  method = "euclidean"))
  
  return(final.dist)
  
}

## Obtain ambient distance from ambient distance matrix
intrinsic.dist = function(distmat, L = NULL){
  
  # distmat : k x k symmetric distance matrix
  # L : number of adjacent nodes connected
  
  
  if(!is.matrix(distmat)){
    
    stop('distmat must be a matrix')
  
  }
  
  if(is.null(L)){
    
    stop("Select the number of adjacent nodes L")
  
  }
  
  # k: number of observations
  k = nrow(distmat)
  
  for(i in 1:k){
    
    distmat[i,][distmat[i,]>sort(distmat[i,])[L]] = 0
    
  }
  
  G <- graph_from_adjacency_matrix(distmat, mode='undirected', weighted = TRUE)
  final.dist = shortest.paths(G, v = V(G), to = V(G),algorithm = "dijkstra")
  
  return(final.dist)
  
}

# Main function: Estimate the intrinsic variance, curvature and test statistics
intrinsic.curv.est = function(df = NULL, distmat=NULL, L = NULL, pop = FALSE) {
  # df : (k x p) data frame on R^p
  # L : number of adjacent nodes connected
  
  if(is.null(df) & is.null(distmat)){
    
    stop("Inputs should be either data frame or distance matrix")
  }
  
  if(is.null(L)){
    
    stop("Select the number of adjacent nodes L")
    
  }
  
  if(!is.null(df)){
    ## Obtain ambient distance from data frame 
    amb.dist.mat = ambient.dist(df)
    
    ## Obtain intrinsic distance from ambient distance matrix
    intrin.dist.mat = intrinsic.dist(amb.dist.mat, L = L)
  } else{
    
    intrin.dist.mat = intrinsic.dist(distmat, L = L)
    
  }
  
  # k: number of observations
  k = nrow(intrin.dist.mat)
  
  # Fmean: Index of Frechet mean
  Fmean = which.min(colMeans(intrin.dist.mat^2))
  
  if(pop){
    # Vf: Frechet variance estimate
    # Vm: Metric variance estimate
    Vf = sum(intrin.dist.mat[Fmean,]^2)/k
    Vm = sum(intrin.dist.mat^2)/(2*k*k)
  } else{
    # Vf: Frechet variance estimate
    # Vm: Metric variance estimate
    Vf = sum(intrin.dist.mat[Fmean,]^2)/k
    Vm = sum(intrin.dist.mat^2)/(2*k*(k-1))
  }
  
  
  # sigma.Vf: asymptotic variance of Frechet variance
  # sigma.Vm: asymptotic variance of Metric variance
  # sigma.Vmf: asymptotic variance of cross product of Vf and Vm
  sigma.Vf = mean(intrin.dist.mat[Fmean,]^4) - Vf^2
  
  asym.dist = 0
  for(i in 1:k){
    asym.dist = asym.dist +(sum(intrin.dist.mat[i,]^2)/(k-1))^2
  }
  sigma.Vm = asym.dist/k - (2*Vm)^2
  
  coasym.dist = 0
  for(i in 1:k){
    coasym.dist = coasym.dist + sum(intrin.dist.mat[i,]^2)/(k-1)*(intrin.dist.mat[Fmean,i]^2)
  }
  sigma.Vmf = coasym.dist/k - 2*Vf*Vm
  
  cov = matrix(data = c(sigma.Vm,sigma.Vmf,sigma.Vmf,sigma.Vf),nrow=2)
  cov.normalized = matrix(data = c(sigma.Vm,sigma.Vmf,sigma.Vmf,sigma.Vf)/k,nrow=2)
  
  rho = Vf/Vm - 1
  
  a = c(-Vf/(Vm^2),1/Vm)
  
  sd = sqrt(a %*% cov %*% a)
  
  test.stat = sqrt(k)*rho/sd
  
  pval = 2*min(1-pnorm(test.stat), pnorm(test.stat))
  
  return(list(Fmean = Fmean, Vf = Vf, Vm = Vm, cov = cov, cov.normalized = cov.normalized, test.stat= test.stat, pval = pval, rho= rho, sd = sd))
  
}


# Main function: Estimate the variances, curvature and test statistics when we have SPD matrices as inputs
cov.curv.est = function(cov.array, metric = NULL) {
  
  if (!metric %in% c("AIRM","Cholesky","Euclidean","LERM","Procrustes.SS","RootEuclidean")){
    stop("metric choice not supported.")
  }
  
  # k: number of observations
  k = dim(cov.array)[3]
  
  # Fmean: Frechet mean
  Fmean = CovMean(cov.array, method=metric)
  
  samples = abind(cov.array,Fmean)
  
  dist.mat = CovDist(samples, method=metric)
  
  # Vf: Frechet variance estimate
  # Vm: Metric variance estimate
  Vf = sum(dist.mat[k+1,1:k]^2)/k
  Vm = sum(dist.mat[1:k,1:k]^2)/(2*k*k)
  
  # sigma.Vf: asymptotic variance of Frechet variance
  # sigma.Vm: asymptotic variance of Metric variance
  # sigma.Vmf: asymptotic variance of cross product of Vf and Vm
  sigma.Vf = mean(dist.mat[k+1,1:k]^4) - Vf^2
  
  asym.dist = 0
  for(i in 1:k){
    asym.dist = asym.dist +(sum(dist.mat[i,1:k]^2)/(k-1))^2
  }
  sigma.Vm = asym.dist/k - (2*Vm)^2
  
  coasym.dist = 0
  for(i in 1:k){
    coasym.dist = coasym.dist + sum(dist.mat[i,1:k]^2)/(k-1)*(dist.mat[k+1,i]^2)
  }
  sigma.Vmf = coasym.dist/k - 2*Vf*Vm
  
  cov = matrix(data = c(sigma.Vm,sigma.Vmf,sigma.Vmf,sigma.Vf),nrow=2)
  cov.normalized = matrix(data = c(sigma.Vm,sigma.Vmf,sigma.Vmf,sigma.Vf)/k,nrow=2)
  
  rho = Vf/Vm - 1
  
  a = c(-Vf/(Vm^2),1/Vm)
  
  sd = sqrt(a %*% cov %*% a)
  
  test.stat = sqrt(k)*rho/sd
  
  pval = 2*min(1-pnorm(test.stat), pnorm(test.stat))
  
  return(list(Fmean = Fmean, Vf = Vf, Vm = Vm, cov = cov, cov.normalized = cov.normalized, test.stat= test.stat, pval = pval, rho= rho, sd = sd))
  
}


# Main function: Estimate the variances, curvature and test statistics when we have spherical data as inputs
sph.curv.est = function(df) {
  
  # k: number of observations
  k = nrow(df)
  
  # p: dimensions
  p = ncol(df)
  
  # Fmean: Frechet mean
  df.s = matrix(0,nrow=k,ncol=p-1)
  for(i in 1:k){
    df.s[i,] = Trans.sph(df[i,])
  }
  
  Fmean = Trans.Euclid(IntrinsicMean(df.s))
  
  samples = rbind(df,Fmean)
  
  dist.mat = matrix(0,nrow=k+1,ncol=k+1)
  for(row in 1:(k+1)){
    for(col in 1:(k+1)){
      
      if (row == col){
        
        dist.mat[row,col] = 0
        
      } else{
        
        dist.mat[row,col] = acos(samples[row,] %*% samples[col,])
        
      }
      
    }
  }
  
  # Vf: Frechet variance estimate
  # Vm: Metric variance estimate
  Vf = sum(dist.mat[k+1,1:k]^2)/k
  Vm = sum(dist.mat[1:k,1:k]^2)/(2*k*(k-1))
  
  # sigma.Vf: asymptotic variance of Frechet variance
  # sigma.Vm: asymptotic variance of Metric variance
  # sigma.Vmf: asymptotic variance of cross product of Vf and Vm
  sigma.Vf = mean(dist.mat[k+1,1:k]^4) - Vf^2
  
  asym.dist = 0
  for(i in 1:k){
    asym.dist = asym.dist +(sum(dist.mat[i,1:k]^2)/(k-1))^2
  }
  sigma.Vm = asym.dist/k - (2*Vm)^2
  
  coasym.dist = 0
  for(i in 1:k){
    coasym.dist = coasym.dist + sum(dist.mat[i,1:k]^2)/(k-1)*(dist.mat[k+1,i]^2)
  }
  sigma.Vmf = coasym.dist/k - 2*Vf*Vm
  
  cov = matrix(data = c(sigma.Vm,sigma.Vmf,sigma.Vmf,sigma.Vf),nrow=2)
  cov.normalized = matrix(data = c(sigma.Vm,sigma.Vmf,sigma.Vmf,sigma.Vf)/k,nrow=2)
  
  rho = Vf/Vm - 1
  
  a = c(-Vf/(Vm^2),1/Vm)
  
  sd = sqrt(a %*% cov %*% a)
  
  test.stat = sqrt(k)*rho/sd
  
  pval = 2*min(1-pnorm(test.stat), pnorm(test.stat))
  
  return(list(Fmean = Fmean, Vf = Vf, Vm = Vm, cov = cov, cov.normalized = cov.normalized, test.stat= test.stat, pval = pval, rho= rho, sd = sd))
  
}

# Main function: Estimate the variances, curvature and test statistics when we have functional data as inputs
functional.curv.est = function(qSup, qin){
  
  # qSup: A numeric vector holding the grid on [0,1] quantile functions take value on
  # qin : A matrix holding the quantile functions of the response.
  
  # k: number of observations
  k = ncol(qin)
  
  Fmean = rowMeans(qin)
  
  samples = cbind(qin,Fmean)
  
  dist.mat = matrix(0,nrow=k+1,ncol=k+1)
  for(row in 1:(k+1)){
    for(col in 1:(k+1)){
      dist.mat[row,col] = sqrt(trapzRcpp(qSup,(samples[,row]-samples[,col])^2))
    }
  }
  
  # Vf: Frechet variance estimate
  # Vm: Metric variance estimate
  Vf = sum(dist.mat[k+1,1:k]^2)/k
  Vm = sum(dist.mat[1:k,1:k]^2)/(2*k*(k-1))
  
  # sigma.Vf: asymptotic variance of Frechet variance
  # sigma.Vm: asymptotic variance of Metric variance
  # sigma.Vmf: asymptotic variance of cross product of Vf and Vm
  sigma.Vf = mean(dist.mat[k+1,1:k]^4) - Vf^2
  
  asym.dist = 0
  for(i in 1:k){
    asym.dist = asym.dist +(sum(dist.mat[i,1:k]^2)/(k-1))^2
  }
  sigma.Vm = asym.dist/k - (2*Vm)^2
  
  coasym.dist = 0
  for(i in 1:k){
    coasym.dist = coasym.dist + sum(dist.mat[i,1:k]^2)/(k-1)*(dist.mat[k+1,i]^2)
  }
  sigma.Vmf = coasym.dist/k - 2*Vf*Vm
  
  cov = matrix(data = c(sigma.Vm,sigma.Vmf,sigma.Vmf,sigma.Vf),nrow=2)
  cov.normalized = matrix(data = c(sigma.Vm,sigma.Vmf,sigma.Vmf,sigma.Vf)/k,nrow=2)
  
  rho = Vf/Vm - 1
  
  return(list(Fmean = Fmean, Vf = Vf, Vm = Vm, cov = cov, cov.normalized = cov.normalized, rho= rho))
  
}

fisher_rao_distance <- function(Sigma1, Sigma2) {
  # Compute the generalized eigenvalues (eigenvalues of Sigma1^(-1) Sigma2)
  eigen_vals <- eigen(solve(Sigma1) %*% Sigma2)$values
  
  # Compute the Fisher-Rao distance
  dist <- sqrt(sum(log(eigen_vals)^2))
  
  return(dist)
}

# Function to compute 2-Wasserstein distance between two zero-mean multivariate Gaussian distributions
wasserstein_distance <- function(Sigma1, Sigma2) {
  # Calculate the square root of Sigma1
  sqrt_Sigma1 <- sqrtm(Sigma1)
  
  # Compute the matrix product: sqrt(Sigma1) * Sigma2 * sqrt(Sigma1)
  A <- sqrt_Sigma1 %*% Sigma2 %*% sqrt_Sigma1
  
  # Calculate the square root of A
  sqrt_A <- sqrtm(A)
  
  # Calculate the trace of each component
  trace_Sigma1 <- sum(diag(Sigma1))
  trace_Sigma2 <- sum(diag(Sigma2))
  trace_sqrt_A <- sum(diag(sqrt_A))
  
  # Compute the squared Wasserstein distance
  W2_squared <- trace_Sigma1 + trace_Sigma2 - 2 * trace_sqrt_A
  
  # Return the Wasserstein distance
  W2 <- sqrt(W2_squared)
  return(W2)
}

wasserstein_geodesic <- function(Sigma1, Sigma2, t) {
  
  sqrt_Sigma2 = sqrtm(Sigma2)
  
  K = sqrt_Sigma2 %*% Sigma1 %*% sqrt_Sigma2
  
  sqrt_K = sqrtm(K)
  
  A =  sqrt_Sigma2 %*% solve(sqrt_K) %*% sqrt_Sigma2
  
  W_geodesic = ( (1-t)*diag(2) + t*A ) %*% Sigma1 %*% ( (1-t)*diag(2) + t*A )
  
  return(W_geodesic)
  
}

wasserstein_barycenter = function(Sigma_list,  lambda, tol = 1e-6, max_iter = 1000) {
  
  sqrt_Sigma_init = matrix(c(0,0,0,0),ncol=2)
  k = length(Sigma_list)
  for (i in 1:k) {
    sqrt_Sigma_init = sqrt_Sigma_init + lambda[i] * sqrtm(Sigma_list[[i]])
  }
  sqrt_Sigma_old = sqrt_Sigma_init
  
  for (iter in 1:max_iter) {
    # Compute the weighted sum of square roots
    sqrt_terms <- lapply(1:k, function(i) {
      sqrtm(sqrt_Sigma_old %*% Sigma_list[[i]] %*% sqrt_Sigma_old)
    })
    
    # Update Sigma
    Sigma_new <- Reduce('+', lapply(1:k, function(i) lambda[i] * sqrt_terms[[i]]))
    sqrt_Sigma_new = sqrtm(Sigma_new)
    
    # Check for convergence
    if (norm(sqrt_Sigma_new - sqrt_Sigma_old, type = "F") < tol) {
      cat("Converged in", iter, "iterations.\n")
      return(Sigma_new)
    }
    
    sqrt_Sigma_old <- sqrt_Sigma_new
  }
  
}

kerFctn <- function(kernel_type){
  if (kernel_type=='gauss'){
    ker <- function(x){
      dnorm(x) #exp(-x^2 / 2) / sqrt(2*pi)
    }
  } else if(kernel_type=='rect'){
    ker <- function(x){
      as.numeric((x<=1) & (x>=-1))
    }
  } else if(kernel_type=='epan'){
    ker <- function(x){
      n <- 1
      (2*n+1) / (4*n) * (1-x^(2*n)) * (abs(x)<=1)
    }
  } else if(kernel_type=='gausvar'){
    ker <- function(x) {
      dnorm(x)*(1.25-0.25*x^2)
    }
  } else if(kernel_type=='quar'){
    ker <- function(x) {
      (15/16)*(1-x^2)^2 * (abs(x)<=1)
    }
  } else {
    stop('Unavailable kernel')
  }
  return(ker)
}

Kern <- kerFctn('gauss')

K <- function(x, h) {
  k <- 1
  for (i in 1:length(h)) {
    k <- k * Kern(x[i] / h[i])
  }
  return(as.numeric(k))
}

find_weight = function(a1, a2, t, isomap_rep, h) {
  
  score = (1-t)* a1 + t *a2
  weight = c()
  
  if (length(a1) == 1){
    
    k = length(isomap_rep)
    
    for(i in 1:k){
      weight = c(weight, K(isomap_rep[i]-score,h))
    }
    
  } else if (length(a1) > 1){
    
    k = nrow(isomap_rep)
    
    for(i in 1:k){
      weight = c(weight, K(isomap_rep[i,]-score,h))
    }
    
  }
  
  weight = weight / sum(weight)
  
  
  return(weight)
  
}


find_weighted_barycenter_on_sphere <- function(points, weights, max_iter = 100, tol = 1e-6) {
  # points: matrix of points on the unit sphere (rows are points, columns are coordinates)
  # weights: vector of weights corresponding to each point
  # max_iter: maximum number of iterations
  # tol: tolerance to stop the iteration when changes are small
  
  # Function to normalize a vector to lie on the unit sphere
  normalize <- function(v) {
    return(v / sqrt(sum(v^2)))
  }
  
  # Ensure that weights sum to 1
  weights <- weights / sum(weights)
  
  # Initial guess: use the weighted centroid normalized to the sphere
  weighted_mean <- colSums(points * weights)
  current_mean <- normalize(weighted_mean)
  
  for (i in 1:max_iter) {
    # Project each point onto the tangent space at the current mean
    tangent_vectors <- t(apply(points, 1, function(p) {
      # Find the tangent vector by subtracting the component in the direction of current_mean
      projection <- p - sum(p * current_mean) * current_mean
      return(weights[which(points == p, arr.ind = TRUE)[1]] * projection)
    }))
    
    # Calculate the weighted average tangent vector
    avg_tangent <- colSums(tangent_vectors)
    
    # If the average tangent vector is near zero, convergence is reached
    if (sqrt(sum(avg_tangent^2)) < tol) {
      break
    }
    
    # Move along the tangent vector and re-project to the sphere
    new_mean <- normalize(current_mean + avg_tangent)
    
    # Check for convergence
    if (sqrt(sum((new_mean - current_mean)^2)) < tol) {
      break
    }
    
    # Update current mean
    current_mean <- new_mean
  }
  
  return(current_mean)
}
