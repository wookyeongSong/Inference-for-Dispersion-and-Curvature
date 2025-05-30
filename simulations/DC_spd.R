# -------------------------------------------------------------------------
# Additional simulation results: Spaces of symmetric positive definite matrices
#
# Manuscript reference : Section S.1.1
# Figure reproduced    : Figure S.1
# -------------------------------------------------------------------------


## Import main functions
source("src/DC_mainfunctions.R")
source("simulations/DC_simulations_datagen.R")


## Generate Simulation Data
k = 100 # Number of samples 

## Generate simulation data from space of symmetric positive definite matrices
cov.array = dat.gen.spd(k)

# (a) Frobenius metric
cov.curv.fit.Frob = cov.curv.est(cov.array, metric = "Euclidean")

# (b) log-Euclidean metric
cov.curv.fit.LE = cov.curv.est(cov.array, metric = "LERM")

# (c) power Frobenius metric with power 1/2
cov.curv.fit.Root = cov.curv.est(cov.array, metric = "RootEuclidean")

# (d) Cholesky metric
cov.curv.fit.Chol = cov.curv.est(cov.array, metric = "Cholesky")

# (e) Affine-Invariant Riemannian (AIR) metric
cov.curv.fit.AIRM = cov.curv.est(cov.array, metric = "AIRM")

# (f) Procrustes size-and-shape (PSS) metric (a.k.a Bures-Wasserstein (BW) metric)
cov.curv.fit.PSS = cov.curv.est(cov.array, metric = "Procrustes.SS")

## Ratio between metric variance over Frechet variance 
Vm.spd = c(cov.curv.fit.Frob$Vm, cov.curv.fit.LE$Vm, cov.curv.fit.Root$Vm, cov.curv.fit.Chol$Vm, cov.curv.fit.AIRM$Vm, cov.curv.fit.PSS$Vm)
Vf.spd = c(cov.curv.fit.Frob$Vf, cov.curv.fit.LE$Vf, cov.curv.fit.Root$Vf, cov.curv.fit.Chol$Vf, cov.curv.fit.AIRM$Vf, cov.curv.fit.PSS$Vf)


## Plot confidence region for the joint distribution of metric and Frechet variance with PSS metric and AIR metric
cov.est.AIRM = cov.curv.fit.AIRM$cov.normalized
Vm.est.AIRM = cov.curv.fit.AIRM$Vm
Vf.est.AIRM = cov.curv.fit.AIRM$Vf

cov.est.PSS = cov.curv.fit.PSS$cov.normalized
Vm.est.PSS = cov.curv.fit.PSS$Vm
Vf.est.PSS = cov.curv.fit.PSS$Vf
