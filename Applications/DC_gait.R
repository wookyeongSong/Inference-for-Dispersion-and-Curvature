# -------------------------------------------------------------------------
# Data Application: Gait synchronization analysis
#
# Data available at https://github.com/deepcharles/gait-data
# Reference: Truong et al. 2019, 'A data set for the study of human locomotion with inertial measurements units', Image Processing On Line 9, 381–390."
# Manuscript reference : Section 6.1, Section S.6
# Figure reproduced    : Figure 6, Figure S.7
# Table reproduced     : Table 2
# -------------------------------------------------------------------------


## Import main functions
source("src/DC_mainfunctions.R")

## Path to data
path = "set_directory_path_to_data"


############################
## Data Preprocessing 
filenames.csv = c()
filenames.json = c()

### NA list
na.subjects = c(67, 68, 128, 155, 158, 201)

### Generate the filenames
for (i in 1:230) {
  
  if (i %in% na.subjects) {next}
  
  filenames.csv = c(filenames.csv, paste0(i, "-1.csv"))
  filenames.json = c(filenames.json, paste0(i, "-1.json"))
  
}

k = 224
gait.list = list(signal = NULL, metadata =NULL)
gait.list$signal = list()
gait.list$metadata = list()
for(i in 1:k) {
  
  gait.list$signal[[i]] = read.csv(paste0(path,filenames.csv[i],sep=""))
  gait.list$metadata[[i]] = fromJSON(paste0(path,filenames.json[i],sep=""))
  
}

path.group = c()
for(i in 1:k){
  
  path.group = c(path.group, gait.list$metadata[[i]]$PathologyGroup)
  
}

healthy_idx = which(path.group == "Healthy")
ortho_idx = which(path.group == "Orthopedic")


## Healthy Group Data
healthy.gait.list = list()
for(i in 1:length(healthy_idx)){
  
  healthy.gait.list[[i]] = cov((gait.list$signal[[healthy_idx[i]]][,c(5,8)]))
  
}

## Orthopedic Group Data
ortho.gait.list = list()
for(i in 1:length(ortho_idx)){
  
  ortho.gait.list[[i]] = cov((gait.list$signal[[ortho_idx[i]]][,c(5,8)]))
  
}



################################################################################
### Section 6.1: Curvature Inference with Bures-Wasserstein ambient distance ###
################################################################################

## Generate Bures-Wasserstein distance matrix of healthy and orthopedic group
healthy.k = length(healthy_idx)
ortho.k = length(ortho_idx)

healthy.dist.mat = matrix(0,nrow=healthy.k,ncol=healthy.k)
for(row in 1:healthy.k){
  for(col in 1:healthy.k){
    
    if(row == col){
      healthy.dist.mat[row,col] = 0
      next
    }
    
    healthy.dist.mat[row,col] =  wasserstein_distance(healthy.gait.list[[row]], healthy.gait.list[[col]])
    
  }
}

ortho.dist.mat = matrix(0,nrow=ortho.k,ncol=ortho.k)
for(row in 1:ortho.k){
  for(col in 1:ortho.k){
    
    if(row == col){
      ortho.dist.mat[row,col] = 0
      next
    }
    
    ortho.dist.mat[row,col] =  wasserstein_distance(ortho.gait.list[[row]], ortho.gait.list[[col]])
    
  }
}


## Intrinsic curvature of healthy and orthopedic group
healthy.intrin.curv = intrinsic.curv.est(distmat = healthy.dist.mat, L = 4)
ortho.intrin.curv = intrinsic.curv.est(distmat = ortho.dist.mat, L = 4)

healthy.isomap.rep <- isomap(healthy.dist.mat, k = 4, ndim = 3)
ortho.isomap.rep <- isomap(ortho.dist.mat, k = 4, ndim = 3)

healthy.isomap.rep = as.matrix(healthy.isomap.rep$points)
ortho.isomap.rep = as.matrix(ortho.isomap.rep$points)

healthy.idx.min = as.vector(which.min(healthy.isomap.rep[,1]))
healthy.idx.max = as.vector(which.max(healthy.isomap.rep[,1]))
ortho.idx.max = as.vector(which.min(ortho.isomap.rep[,1]))
ortho.idx.min = as.vector(which.max(ortho.isomap.rep[,1]))

k=healthy.k
## Confidence intervals for the intrinsic curvatures rho_{I}
# alpha = 0.01 (0.072, 0.132)
C.I.lower.healthy.01 = healthy.intrin.curv$rho - qnorm(0.995)*healthy.intrin.curv$sd/sqrt(k)
C.I.upper.healthy.01 = healthy.intrin.curv$rho + qnorm(0.995)*healthy.intrin.curv$sd/sqrt(k)

# alpha = 0.05 (0.080, 0.125)
C.I.lower.healthy.05 = healthy.intrin.curv$rho - qnorm(0.975)*healthy.intrin.curv$sd/sqrt(k)
C.I.upper.healthy.05 = healthy.intrin.curv$rho + qnorm(0.975)*healthy.intrin.curv$sd/sqrt(k)

# alpha = 0.1 (0.083, 0.122)
C.I.lower.healthy.1 = healthy.intrin.curv$rho - qnorm(0.95)*healthy.intrin.curv$sd/sqrt(k)
C.I.upper.healthy.1 = healthy.intrin.curv$rho + qnorm(0.95)*healthy.intrin.curv$sd/sqrt(k)


## Table 2
ci.tbl.healthy <- data.frame(
  Level = c("99%", "95%", "90%"),
  Lower = c(C.I.lower.healthy.01,
            C.I.lower.healthy.05,
            C.I.lower.healthy.1),
  Upper = c(C.I.upper.healthy.01,
            C.I.upper.healthy.05,
            C.I.upper.healthy.1)
)

ci.tbl.healthy <- transform(ci.tbl.healthy,
                        Lower = round(Lower, 3),
                        Upper = round(Upper, 3))

print("Healthy Group: Confidence interval for intrinsic curvature")
print(ci.tbl.healthy, row.names = FALSE)


k=ortho.k
## Confidence intervals for the intrinsic curvatures rho_{I}
# alpha = 0.01 (0.072, 0.132)
C.I.lower.ortho.01 = ortho.intrin.curv$rho - qnorm(0.995)*ortho.intrin.curv$sd/sqrt(k)
C.I.upper.ortho.01 = ortho.intrin.curv$rho + qnorm(0.995)*ortho.intrin.curv$sd/sqrt(k)

# alpha = 0.05 (0.080, 0.125)
C.I.lower.ortho.05 = ortho.intrin.curv$rho - qnorm(0.975)*ortho.intrin.curv$sd/sqrt(k)
C.I.upper.ortho.05 = ortho.intrin.curv$rho + qnorm(0.975)*ortho.intrin.curv$sd/sqrt(k)

# alpha = 0.1 (0.083, 0.122)
C.I.lower.ortho.1 = ortho.intrin.curv$rho - qnorm(0.95)*ortho.intrin.curv$sd/sqrt(k)
C.I.upper.ortho.1 = ortho.intrin.curv$rho + qnorm(0.95)*ortho.intrin.curv$sd/sqrt(k)


## Table 2
ci.tbl.ortho <- data.frame(
  Level = c("99%", "95%", "90%"),
  Lower = c(C.I.lower.ortho.01,
            C.I.lower.ortho.05,
            C.I.lower.ortho.1),
  Upper = c(C.I.upper.ortho.01,
            C.I.upper.ortho.05,
            C.I.upper.ortho.1)
)

ci.tbl.ortho <- transform(ci.tbl.ortho,
                            Lower = round(Lower, 3),
                            Upper = round(Upper, 3))

print("Orthopedic Group: Confidence interval for intrinsic curvature")
print(ci.tbl.ortho, row.names = FALSE)


## Plot confidence region for the joint distribution of intrinsic metric variance and Frechet variance
cov.est.healthy = healthy.intrin.curv$cov.normalized
Vm.est.healthy = healthy.intrin.curv$Vm
Vf.est.healthy = healthy.intrin.curv$Vf

cov.est.ortho = ortho.intrin.curv$cov.normalized
Vm.est.ortho = ortho.intrin.curv$Vm
Vf.est.ortho = ortho.intrin.curv$Vf



################################################################################################################
### Section S.6: Sensitivity Analysis of Input distance: Bures-Wasserstein, Cholesky, and Frobenius distance ###
################################################################################################################

## Generate Bures-Wasserstein distance matrix of healthy and orthopedic group
healthy.k = length(healthy_idx)
ortho.k = length(ortho_idx)

healthy.dist.mat.BW = matrix(0,nrow=healthy.k,ncol=healthy.k)
for(row in 1:healthy.k){
  for(col in 1:healthy.k){
    
    if(row == col){
      healthy.dist.mat.BW[row,col] = 0
      next
    }
    
    healthy.dist.mat.BW[row,col] =  wasserstein_distance(healthy.gait.list[[row]], healthy.gait.list[[col]])
    
  }
}

ortho.dist.mat.BW = matrix(0,nrow=ortho.k,ncol=ortho.k)
for(row in 1:ortho.k){
  for(col in 1:ortho.k){
    
    if(row == col){
      ortho.dist.mat.BW[row,col] = 0
      next
    }
    
    ortho.dist.mat.BW[row,col] =  wasserstein_distance(ortho.gait.list[[row]], ortho.gait.list[[col]])
    
  }
}


## Generate Frobenius distance matrix of healthy and orthopedic group
healthy.dist.mat.Fb = matrix(0,nrow=healthy.k,ncol=healthy.k)
for(row in 1:healthy.k){
  for(col in 1:healthy.k){
    
    if(row == col){
      healthy.dist.mat.Fb[row,col] = 0
      next
    }
    
    healthy.dist.mat.Fb[row,col] = sqrt(sum((healthy.gait.list[[row]]-healthy.gait.list[[col]])^2)) 
    
  }
}

ortho.dist.mat.Fb = matrix(0,nrow=ortho.k,ncol=ortho.k)
for(row in 1:ortho.k){
  for(col in 1:ortho.k){
    
    if(row == col){
      ortho.dist.mat.Fb[row,col] = 0
      next
    }
    
    ortho.dist.mat.Fb[row,col] =  sqrt(sum((ortho.gait.list[[row]]-ortho.gait.list[[col]])^2)) 
    
  }
}


## Generate Cholesky distance matrix of healthy and orthopedic group
healthy.dist.mat.Chol = matrix(0,nrow=healthy.k,ncol=healthy.k)
for(row in 1:healthy.k){
  for(col in 1:healthy.k){
    
    if(row == col){
      healthy.dist.mat.Chol[row,col] = 0
      next
    }
    
    healthy.dist.mat.Chol[row,col] = sqrt(sum((chol(healthy.gait.list[[row]])-chol(healthy.gait.list[[col]]))^2)) 
    
  }
}


ortho.dist.mat.Chol = matrix(0,nrow=ortho.k,ncol=ortho.k)
for(row in 1:ortho.k){
  for(col in 1:ortho.k){
    
    if(row == col){
      ortho.dist.mat.Chol[row,col] = 0
      next
    }
    
    ortho.dist.mat.Chol[row,col] =  sqrt(sum((chol(ortho.gait.list[[row]])-chol(ortho.gait.list[[col]]))^2)) 
    
  }
}


## Implement intrinsic curvature test with Bures-Wasserstein metric

healthy.intrin.curv.BW = intrinsic.curv.est(distmat = healthy.dist.mat.BW, L = 4)
ortho.intrin.curv.BW = intrinsic.curv.est(distmat = ortho.dist.mat.BW, L = 4)

cov.est.healthy.BW = healthy.intrin.curv.BW$cov.normalized
Vm.est.healthy.BW = healthy.intrin.curv.BW$Vm
Vf.est.healthy.BW = healthy.intrin.curv.BW$Vf

cov.est.ortho.BW = ortho.intrin.curv.BW$cov.normalized
Vm.est.ortho.BW = ortho.intrin.curv.BW$Vm
Vf.est.ortho.BW = ortho.intrin.curv.BW$Vf


## Implement intrinsic curvature test with Frobenius metric

healthy.intrin.curv.Fb = intrinsic.curv.est(distmat = healthy.dist.mat.Fb, L = 4)
ortho.intrin.curv.Fb = intrinsic.curv.est(distmat = ortho.dist.mat.Fb, L = 4)

cov.est.healthy.Fb = healthy.intrin.curv.Fb$cov.normalized
Vm.est.healthy.Fb = healthy.intrin.curv.Fb$Vm
Vf.est.healthy.Fb = healthy.intrin.curv.Fb$Vf

cov.est.ortho.Fb = ortho.intrin.curv.Fb$cov.normalized
Vm.est.ortho.Fb = ortho.intrin.curv.Fb$Vm
Vf.est.ortho.Fb = ortho.intrin.curv.Fb$Vf


## Implement intrinsic curvature test with Cholesky metric

healthy.intrin.curv.Chol = intrinsic.curv.est(distmat = healthy.dist.mat.Chol, L = 4)
ortho.intrin.curv.Chol = intrinsic.curv.est(distmat = ortho.dist.mat.Chol, L = 4)

cov.est.healthy.Chol = healthy.intrin.curv.Chol$cov.normalized
Vm.est.healthy.Chol = healthy.intrin.curv.Chol$Vm
Vf.est.healthy.Chol = healthy.intrin.curv.Chol$Vf

cov.est.ortho.Chol = ortho.intrin.curv.Chol$cov.normalized
Vm.est.ortho.Chol = ortho.intrin.curv.Chol$Vm
Vf.est.ortho.Chol = ortho.intrin.curv.Chol$Vf

