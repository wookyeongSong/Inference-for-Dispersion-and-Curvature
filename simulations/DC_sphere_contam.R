# -------------------------------------------------------------------------
# Simulation Study: Noise-Contaminated Spherical Data
#
# Manuscript reference : Section 4.2
# Figure reproduced    : Figure 2
# -------------------------------------------------------------------------


## Import main functions
source("src/DC_mainfunctions.R")


## Generate simulation data
k=1000

set.seed(9999)
theta = pi * runif(k) / 2
psi = 2 * pi * runif(k)
z <- cos(theta) 
x <- sin(theta)*cos(psi) 
y <- sin(theta)*sin(psi)

## Intrinsic curvature test for Noise-Contaminated Spherical Data with different error sigma
sd.vec = c(1/32, 1/16, 1/8, 1/4, 1/2)
rho.vec = c()
lower.vec = c()
upper.vec = c()

sph.list = list()
for(i in 1:length(sd.vec)){
  
  sd = sd.vec[i]
  
  set.seed(i)
  
  z.new <- z + rnorm(k, mean = 0, sd =sd)
  x.new <- x + rnorm(k, mean = 0, sd =sd) 
  y.new <- y + rnorm(k, mean = 0, sd =sd)
  
  df.sph = as.matrix(data.frame(x.new,y.new,z.new))
  sph.list[[i]] = df.sph
  
  fit = intrinsic.curv.est(df = df.sph, L = 5)
  
  rho.vec = c(rho.vec, fit$rho)
  lower.vec = c(lower.vec, fit$rho - qnorm(0.975)*fit$sd/sqrt(k))
  upper.vec = c(upper.vec, fit$rho + qnorm(0.975)*fit$sd/sqrt(k))
}
