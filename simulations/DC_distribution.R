# -------------------------------------------------------------------------
# Simulation Study: Distributional data with 2-Wasserstein distance
#
# Manuscript reference : Section 5.1, Section S.5, Section S.6
# Figures reproduced    : Figure 3, Figure 4, Figure S.6, Figure S.8
# -------------------------------------------------------------------------


## Import main functions
source("src/DC_mainfunctions.R")


## Generate simulation data
base_cov_mat = matrix(c(4,0,0,1), nrow = 2)

set.seed(999)
k = 100

theta = pi/2 * rbeta(k, 2, 2)
rot = lapply(1:k, function(i){
  matrix(c(cos(theta[i]), sin(theta[i]), -sin(theta[i]), cos(theta[i])), nrow = 2)
})

cov_mat = lapply(1:k, function(i){
  rot[[i]] %*% base_cov_mat %*% t(rot[[i]])
})

base_left = matrix(c(4,0,0,1), nrow = 2)
base_right = matrix(c(1,0,0,4), nrow = 2)

cov_mat[[k+1]] = base_left
cov_mat[[k+2]] = base_right

min_idx = k+1
max_idx = k+2


#################################################################################################
### Section 5.1: Intrinsic curvature test for Distributional data with 2-Wasserstein distance ###
#################################################################################################

## Obtain 2-Wasserstein pairwise distance
dist.mat.Wass = matrix(0,nrow=(k+2),ncol=(k+2)) 
for(row in 1:(k+2)){
  for(col in 1:(k+2)){
    
    if(row == col){
      dist.mat.Wass[row,col] = 0
      next
    }
    
    dist.mat.Wass[row,col] =  wasserstein_distance(cov_mat[[row]], cov_mat[[col]])
    
  }
}

## Intrinsic curvature test
gauss_wass_func.Wass = intrinsic.curv.est(distmat = dist.mat.Wass, L = 6, pop=TRUE)
cov.est.wass = gauss_wass_func.Wass$cov.normalized
Vm.est.wass = gauss_wass_func.Wass$Vm
Vf.est.wass = gauss_wass_func.Wass$Vf

## ISOMAP results
isomap_result <- isomap(dist.mat.Wass, k = 6, ndim = 1)
isomap_rep = as.vector(isomap_result$points)


#############################################################################################
### Section S.6: Sensitivity Analysis of Input distance: Wasserstein, Fisher-Rao distance ###
#############################################################################################

## Obtain Fisher-Rao pairwise distance
dist.mat.FR = matrix(0,nrow=(k+2),ncol=(k+2))
for(row in 1:(k+2)){
  for(col in 1:(k+2)){
    
    if(row == col){
      dist.mat.FR[row,col] = 0
      next
    }
    
    dist.mat.FR[row,col] =  fisher_rao_distance(cov_mat[[row]], cov_mat[[col]])
    
  }
}

gauss_wass_func.FR = intrinsic.curv.est(distmat = dist.mat.FR, L = 6, pop=TRUE)

## Plot confidence region for the joint distribution of metric and Frechet variance when we have spherical data with geodesic distance
cov.est.wass.FR = gauss_wass_func.FR$cov.normalized
Vm.est.wass.FR = gauss_wass_func.FR$Vm
Vf.est.wass.FR = gauss_wass_func.FR$Vf



