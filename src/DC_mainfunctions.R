library(igraph)
library(ellipse)
library(CovTools)
library(spherepc)
library(abind)
library(fdapace)
library(colorspace)
library(plotrix)

## Obtain ambient distance from data frame 
ambient.dist = function(df){
  
  # df : (k x p) data frame on R^p
  
  k = nrow(df)
  final.dist = matrix(NA,nrow = k, ncol= k)
  for(row in 1:k){
    for(col in 1:k){
      
      final.dist[row,col] = sqrt(sum((df[row,] - df[col,])^2))
      
    }
  }
  
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
intrinsic.curv.est = function(df, L = NULL) {
  # df : (k x p) data frame on R^p
  # L : number of adjacent nodes connected
  
  if(is.null(L)){
    
    stop("Select the number of adjacent nodes L")
    
  }
  
  ## Obtain ambient distance from data frame 
  amb.dist.mat = ambient.dist(df)
  
  ## Obtain intrinsic distance from ambient distance matrix
  intrin.dist.mat = intrinsic.dist(amb.dist.mat, L = L)
  
  # k: number of observations
  k = nrow(intrin.dist.mat)
  
  # Fmean: Index of Frechet mean
  Fmean = which.min(colMeans(intrin.dist.mat^2))
  
  # Vf: Frechet variance estimate
  # Vm: Metric variance estimate
  Vf = sum(intrin.dist.mat[Fmean,]^2)/k
  Vm = sum(intrin.dist.mat^2)/(2*k*(k-1))
  
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