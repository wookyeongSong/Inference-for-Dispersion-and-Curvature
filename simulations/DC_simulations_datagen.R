
dat.gen.pc.pos = function(k) {
  
  # k : the number of samples from point cloud on spherical manifold
  set.seed(999)
  theta = pi * runif(k) / 2
  psi = 2 * pi * runif(k)
  
  z <- cos(theta) + rnorm(1, mean = 0, sd = 1/10)
  x <- sin(theta)*cos(psi)  + rnorm(1, mean = 0, sd = 1/10)
  y <- sin(theta)*sin(psi)  + rnorm(1, mean = 0, sd = 1/10)
  
  df.sph = as.matrix(data.frame(x,y,z))
  
  return(df.sph)
  
}

## Generate simulation data from point cloud on negatively curved manifold
dat.gen.pc.neg = function(k) {
  
  # k : the number of samples from point cloud on spherical manifold
  set.seed(1129)
  x = rnorm(k, mean = 0, sd = 1)
  
  y = sapply(1:k, function(i){
    
    set.seed(i)
    theta = runif(1,min=0, max = pi)
    sqrt(1+x[i]^2)*cos(theta) + rnorm(1, mean = 0, sd = 1/10)
    
  })
  
  z = sapply(1:k, function(i){
    
    set.seed(i)
    theta = runif(1,min=0, max = pi)
    sqrt(1+x[i]^2)*sin(theta) + rnorm(1, mean = 0, sd = 1/10)
    
  })

  df.hyper = as.matrix(data.frame(x,y,z))
  
  return(df.hyper)
  
}

## Generate simulation data from point cloud on negatively curved manifold
dat.gen.pc.flat = function(k) {
  
  # k : the number of samples from point cloud on spherical manifold
  set.seed(1129)
  x <- runif(k) +  rnorm(k,mean = 0 , sd = 1/10)    # uniform on [0, 2pi]
  y <- runif(k) +  rnorm(k,mean = 0 , sd = 1/10)
  z = rep(0,k) +  rnorm(k,mean = 0 , sd = 1/10)
  df.plane = as.matrix(data.frame(x,y,z))
  
  return(df.plane)
  
}

## Generate simulation data from space of symmetric positive definite matrices
dat.gen.spd = function(k){
  
  # k : the number of samples from point cloud on spherical manifold
  cov.array = array(0,dim = c(5,5,k))
  for (i in 1:k){
    set.seed(1231+i)
    eig.vec = rnorm(25)
    eig.val = 100*rbeta(5,3,5)
    v1 = eig.vec[1:5]/sqrt(sum(eig.vec[1:5]^2))
    v2 = eig.vec[6:10]/sqrt(sum(eig.vec[6:10]^2))
    v3 = eig.vec[11:15]/sqrt(sum(eig.vec[11:15]^2))
    v4 = eig.vec[16:20]/sqrt(sum(eig.vec[16:20]^2))
    v5 = eig.vec[21:25]/sqrt(sum(eig.vec[21:25]^2))
    V = matrix(c(v1,v2,v3,v4,v5),nrow=5)
    cov.array[,,i] = V %*% diag(eig.val) %*% t(V)
  }
  
  return(cov.array)

}


## Generate simulation data from open upper hemisphere with geodesic metric
dat.gen.sph = function(k){
  
  # k : the number of samples from point cloud on spherical manifold
  set.seed(999)
  z <- runif(k,min=0,max =1)          # uniform on [0, 1]
  theta <- runif(k, min = 0 , max = 2*pi)    # uniform on [0, 2pi]
  x <- cos(theta)*sqrt(1-z^2)  # based on angle
  y <- sin(theta)*sqrt(1-z^2) 
  
  df = as.matrix(data.frame(x,y,z))
  
  return(df)
  
}

## Generate simulation data from point cloud on positively curved manifold
dat.gen.power.anal = function(k, curv) {
  
  # k : the number of samples from point cloud on spherical manifold
  # No error
  
  if ( curv < 0 | curv >2){
    stop("error")
  }
  
  if (curv == 0){
    
    z = rep(0, k)
    x = runif(k, min = -pi/(4*sqrt(2)), max = pi/(4*sqrt(2)))
    y = runif(k, min = -pi/(4*sqrt(2)), max = pi/(4*sqrt(2)))
    df = as.matrix(data.frame(x,y,z))
    L = 50
    
  } else {
    
    radi = 1/sqrt(curv)
    z <- runif(k,min=1/sqrt(curv)*cos(pi * sqrt(curv) / 4),max =1/sqrt(curv))          # uniform on [0, 1]
    theta <- runif(k, min = 0 , max = 2*pi)    # uniform on [0, 2pi]
    x <- cos(theta)*sqrt(1/curv-z^2)  # based on angle
    y <- sin(theta)*sqrt(1/curv-z^2) 
    
    df = as.matrix(data.frame(x,y,z))
    L = 20
  }
  
  return(list(df= df, L = L))
  
}


dat.gen.high.dim = function(k, p, noise0){
  
  # k : the number of samples from point cloud on spherical manifold
  theta = pi * runif(k) / 2
  psi = 2 * pi * runif(k)
  
  noise = noise0
  z <- cos(theta) + rnorm(k, mean = 0, sd = noise )
  x <- sin(theta)*cos(psi)  + rnorm(k, mean = 0, sd =noise)
  y <- sin(theta)*sin(psi)  + rnorm(k, mean = 0, sd = noise)
  
  
  extra = matrix(rnorm((p-3)*k, mean = 0, sd = noise), nrow = k, ncol = p-3)
  df.sph = as.matrix(data.frame(x,y,z))
  
  df = cbind(df.sph,extra)
  
  return(df)
  
}

dat.gen.high.dim.fixed = function(k, p, noise0){
  
  # k : the number of samples from point cloud on spherical manifold
  theta = pi * runif(k) / 2
  psi = 2 * pi * runif(k)
  
  noise = noise0/sqrt(p)
  z <- cos(theta) + rnorm(k, mean = 0, sd = noise )
  x <- sin(theta)*cos(psi)  + rnorm(k, mean = 0, sd =noise)
  y <- sin(theta)*sin(psi)  + rnorm(k, mean = 0, sd = noise)
  
  
  extra = matrix(rnorm((p-3)*k, mean = 0, sd = noise), nrow = k, ncol = p-3)
  df.sph = as.matrix(data.frame(x,y,z))
  
  df = cbind(df.sph,extra)
  
  return(df)
  
}





