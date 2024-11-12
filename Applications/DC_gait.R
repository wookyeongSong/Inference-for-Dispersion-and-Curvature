source("~/src/DC_mainfunctions.R")

library("jsonlite")

path = "path of data"
## Data available at https://github.com/deepcharles/gait-data, 
## Reference: Truong et al. 2019, 'A data set for the study of human locomotion with inertial measurements units', Image Processing On Line 9, 381–390."

##########################
### Data Preprocessing ###
##########################
# Initialize an empty vector to store the filenames
filenames.csv = c()
filenames.json = c()

# NA list
na.subjects = c(67, 68, 128, 155, 158, 201)

# Generate the filenames
for (i in 1:230) {
  
  if (i %in% na.subjects) {
    
    next
    
  }
  
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

### Healthy
### Orthopedic
path.group = c()
for(i in 1:k){
  
  path.group = c(path.group, gait.list$metadata[[i]]$PathologyGroup)
  
}

healthy_idx = which(path.group == "Healthy")
ortho_idx = which(path.group == "Orthopedic")

healthy.gait.list = list()
for(i in 1:length(healthy_idx)){
  
  healthy.gait.list[[i]] = cov((gait.list$signal[[healthy_idx[i]]][,c(5,8)]))
  
}

ortho.gait.list = list()
for(i in 1:length(ortho_idx)){
  
  ortho.gait.list[[i]] = cov((gait.list$signal[[ortho_idx[i]]][,c(5,8)]))
  
}


## Generate distance matrix of healthy and orthopedic group
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



## Plot confidence region for the joint distribution of intrinsic metric variance and Frechet variance
cov.est.healthy = healthy.intrin.curv$cov.normalized
Vm.est.healthy = healthy.intrin.curv$Vm
Vf.est.healthy = healthy.intrin.curv$Vf

cov.est.ortho = ortho.intrin.curv$cov.normalized
Vm.est.ortho = ortho.intrin.curv$Vm
Vf.est.ortho = ortho.intrin.curv$Vf

t=2
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((t+1))

par(mar=c(5,5,3,2))
plot(ellipse(cov.est.healthy,center = c(Vm.est.healthy,Vf.est.healthy),level=0.99), xlim=c(400, 1900), ylim = c(400,1900), type='l',col=col[1] ,lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Healthy Group")
points(Vm.est.healthy,Vf.est.healthy,pch=16)
lines(ellipse(cov.est.healthy,center = c(Vm.est.healthy,Vf.est.healthy),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.healthy,center = c(Vm.est.healthy,Vf.est.healthy),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

par(mar=c(5,5,3,2))
plot(ellipse(cov.est.ortho,center = c(Vm.est.ortho,Vf.est.ortho),level=0.99), xlim=c(200, 1700), ylim = c(200,1700), type='l',col=col[1] ,lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Orthopedic Group")
points(Vm.est.ortho,Vf.est.ortho,pch=16)
lines(ellipse(cov.est.ortho,center = c(Vm.est.ortho,Vf.est.ortho),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.ortho,center = c(Vm.est.ortho,Vf.est.ortho),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

## Plot Intrinsic Geodesic
t=5
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((t))

par(mfrow=c(1,1),mar = c(5,5,3,2))
plot(ellipse(healthy.gait.list[[1]],center = c(0,0)),xlim= c(-400,400), ylim = c(-400,400), type='l',col='grey80',cex.lab=2,cex.axis=2,cex.main=2.5,lwd=2,xlab="", ylab = "",main="Intrinsic Geodesic (Healthy Group)")
for(i in 2:healthy.k){
  lines(ellipse(healthy.gait.list[[i]],center = c(0,0)),col='grey80')
}
lines(ellipse(healthy.gait.list[[healthy.idx.min]]), col = col[1], lwd =4)
lines(ellipse(wasserstein_barycenter(healthy.gait.list, find_weight(healthy.isomap.rep[healthy.idx.min,1],healthy.isomap.rep[healthy.idx.max,1],0.25,healthy.isomap.rep[,1], h = 4))), col = col[2], lwd =4)
lines(ellipse(wasserstein_barycenter(healthy.gait.list, find_weight(healthy.isomap.rep[healthy.idx.min,1],healthy.isomap.rep[healthy.idx.max,1],0.5,healthy.isomap.rep[,1], h = 4))), col = col[3], lwd =4)
lines(ellipse(wasserstein_barycenter(healthy.gait.list, find_weight(healthy.isomap.rep[healthy.idx.min,1],healthy.isomap.rep[healthy.idx.max,1],0.75,healthy.isomap.rep[,1], h = 4))),col = col[4], lwd =4)
lines(ellipse(healthy.gait.list[[healthy.idx.max]]), col = col[5], lwd =4)
legend("topleft", legend = c("t = 0","t = 0.25", "t = 0.5", "t = 0.75", "t = 1"),lty=1,lwd=4, col = col,cex=2)

par(mfrow=c(1,1),mar = c(5,5,3,2))
plot(ellipse(ortho.gait.list[[1]],center = c(0,0)),xlim= c(-400,400), ylim = c(-400,400), type='l',col='grey80',cex.lab=2,cex.axis=2,cex.main=2.5,lwd=2,xlab="", ylab = "",main="Intrinsic Geodesic (Orthopedic Group)")
for(i in 2:ortho.k){
  lines(ellipse(ortho.gait.list[[i]],center = c(0,0)),col='grey80')
}
lines(ellipse(ortho.gait.list[[ortho.idx.min]]), col = col[1], lwd =4)
lines(ellipse(wasserstein_barycenter(ortho.gait.list, find_weight(ortho.isomap.rep[ortho.idx.min,1],ortho.isomap.rep[ortho.idx.max,1],0.25,ortho.isomap.rep[,1], h = 4))), col = col[2], lwd =4)
lines(ellipse(wasserstein_barycenter(ortho.gait.list, find_weight(ortho.isomap.rep[ortho.idx.min,1],ortho.isomap.rep[ortho.idx.max,1],0.5,ortho.isomap.rep[,1], h = 4))), col = col[3], lwd =4)
lines(ellipse(wasserstein_barycenter(ortho.gait.list, find_weight(ortho.isomap.rep[ortho.idx.min,1],ortho.isomap.rep[ortho.idx.max,1],0.75,ortho.isomap.rep[,1], h = 4))),col = col[4], lwd =4)
lines(ellipse(ortho.gait.list[[ortho.idx.max]]), col = col[5], lwd =4)
legend("topleft", legend = c("t = 0","t = 0.25", "t = 0.5", "t = 0.75", "t = 1"),lty=1,lwd=4, col = col,cex=2)




