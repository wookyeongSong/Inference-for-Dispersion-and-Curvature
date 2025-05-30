# -------------------------------------------------------------------------
# Additional simulation results: Spherical data with geodesic distance
#
# Manuscript reference : Section S.1.2
# Figure reproduced    : Figure S.2
# -------------------------------------------------------------------------


## Import main functions
source("src/DC_mainfunctions.R")
source("simulations/DC_simulations_datagen.R")

k = 50 # Number of samples 

## Generate simulation data from open upper hemisphere with geodesic metric
df.sph = dat.gen.sph(k)

# Main function: Estimate the variances, curvature and test statistics when we have spherical data as inputs
sph.curv.fit = sph.curv.est(df.sph)

## Confidence region for the joint distribution of metric and Frechet variance when we have spherical data with geodesic distance
cov.est.sph = sph.curv.fit$cov.normalized
Vm.est.sph = sph.curv.fit$Vm
Vf.est.sph = sph.curv.fit$Vf

