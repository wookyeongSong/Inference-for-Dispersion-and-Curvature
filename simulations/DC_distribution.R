source("~/src/DC_mainfunctions.R")

#######################
### Data Generation ###
#######################
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


################################################################################
### Sensitivity Analysis of Input distance: Wasserstein, Fisher-Rao distance ###
################################################################################


## Calculate the distance
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

gauss_wass_func.Wass = intrinsic.curv.est(distmat = dist.mat.Wass, L = 6, pop=TRUE)

## Plot confidence region for the joint distribution of metric and Frechet variance when we have spherical data with geodesic distance
cov.est.wass = gauss_wass_func.Wass$cov.normalized
Vm.est.wass = gauss_wass_func.Wass$Vm
Vf.est.wass = gauss_wass_func.Wass$Vf

isomap_result <- isomap(dist.mat.Wass, k = 6, ndim = 1)
isomap_rep = as.vector(isomap_result$points)

## Calculate the distance
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



t=2
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((t+1))
par(mfrow = c(1,2), mar=c(5,5,3,2),oma = c(0, 1, 5, 0))

plot(ellipse(cov.est.wass,center = c(Vm.est.wass,Vf.est.wass),level=0.99),type='l',col=col[1],lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2,cex.lab=2,cex.axis=2,main="Input distance: Wasserstein Metric")
points(Vm.est.wass,Vf.est.wass,pch=16)
lines(ellipse(cov.est.wass,center = c(Vm.est.wass,Vf.est.wass),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.wass,center = c(Vm.est.wass,Vf.est.wass),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

plot(ellipse(cov.est.wass.FR,center = c(Vm.est.wass.FR,Vf.est.wass.FR),level=0.99),type='l',col=col[1],lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2,cex.lab=2,cex.axis=2,main="Input distance: Fisher-Rao Metric")
points(Vm.est.wass.FR,Vf.est.wass.FR,pch=16)
lines(ellipse(cov.est.wass.FR,center = c(Vm.est.wass.FR,Vf.est.wass.FR),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.wass.FR,center = c(Vm.est.wass.FR,Vf.est.wass.FR),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

mtext("Distributional Data Simulation", side = 3, font = 2, outer = TRUE, line = 1, cex = 2.5)


#############################################################
### Curvature Inference with Wasserstein ambient distance ###
#############################################################

t=2
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((t+1))
par(mar = c(5, 5, 4, 1))  # Adjust margins for the main plot
plot(ellipse(cov.est.wass,center = c(Vm.est.wass,Vf.est.wass),level=0.99),type='l',col=col[1],lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Proposed intrinsic curvature test")
points(Vm.est.wass,Vf.est.wass,pch=16)
lines(ellipse(cov.est.wass,center = c(Vm.est.wass,Vf.est.wass),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.wass,center = c(Vm.est.wass,Vf.est.wass),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

## Plot ISOMAP Representation
t=5
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((t))
t1 = 0.25
t2 = 0.5
t3 = 0.75
par(mfrow=c(1,1),mar = c(5,5,3,2))
plot(c(theta)/(pi/2), isomap_rep[1:100], xlim = c(0,1), ylim = c(-1.2, 1.2), xlab = expression(theta), ylab = "ISOMAP Representation",cex.lab=2,cex.axis=2,cex.main=2.5)
points(0,isomap_rep[101], col = col[1], pch = 19, cex = 3)
points(find_weight(isomap_rep[min_idx],isomap_rep[max_idx],t1,isomap_rep, h =0.01) %*% c(theta/(pi/2),0,1), (1-t1)*isomap_rep[101] + t1*isomap_rep[102], col = col[2], pch = 19, cex = 3 )
points(find_weight(isomap_rep[min_idx],isomap_rep[max_idx],t2,isomap_rep, h =0.01) %*% c(theta/(pi/2),0,1), (1-t2)*isomap_rep[101] + t2*isomap_rep[102], col = col[3], pch = 19, cex = 3 )
points(find_weight(isomap_rep[min_idx],isomap_rep[max_idx],t3,isomap_rep, h =0.01) %*% c(theta/(pi/2),0,1), (1-t3)*isomap_rep[101] + t3*isomap_rep[102], col = col[4], pch = 19, cex = 3 )
points(1,isomap_rep[102], col = col[5], pch = 19, cex = 3)
legend("topleft", legend = c("t = 0","t = 0.25", "t = 0.5", "t = 0.75", "t = 1"),pch = 19, col = col,cex=1.5)

## Plot Wasserstein Geodesic
par(mfrow=c(1,1),mar = c(5,5,3,2))
plot(ellipse(cov_mat[[1]],center = c(0,0)),xlim= c(-6,6), ylim = c(-6,6), type='l',col='grey80',cex.lab=2,cex.axis=2,cex.main=2.5,lwd=2,xlab="X1", ylab = "X2",main="Wasserstein Geodesic")
for(i in 2:k){
  lines(ellipse(cov_mat[[i]],center = c(0,0)),col='grey80')
}
lines(ellipse(cov_mat[[min_idx]]), col = col[1], lwd =5)
lines(ellipse(wasserstein_geodesic(cov_mat[[min_idx]], cov_mat[[max_idx]], 0.25)), col = col[2], lwd =2)
lines(ellipse(wasserstein_geodesic(cov_mat[[min_idx]], cov_mat[[max_idx]], 0.5)), col = col[3], lwd =5)
lines(ellipse(wasserstein_geodesic(cov_mat[[min_idx]], cov_mat[[max_idx]], 0.75)), col = col[4], lwd =2)
lines(ellipse(cov_mat[[max_idx]]), col = col[5], lwd =5)
legend("topleft", legend = c("t = 0","t = 0.25", "t = 0.5", "t = 0.75", "t = 1"),lty=1,lwd=c(5,2,5,2,5), col = col,cex=1.5)

## Plot Intrinsic Geodesic
par(mfrow=c(1,1),mar = c(5,5,3,2))
plot(ellipse(cov_mat[[1]],center = c(0,0)),xlim= c(-6,6), ylim = c(-6,6), type='l',col='grey80',cex.lab=2,cex.axis=2,cex.main=2.5,lwd=2,xlab="X1", ylab = "X2",main="Intrinsic Geodesic")
for(i in 2:k){
  lines(ellipse(cov_mat[[i]],center = c(0,0)),col='grey80')
}
lines(ellipse(cov_mat[[min_idx]]), col = col[1], lwd =5)
lines(ellipse(wasserstein_barycenter(cov_mat, find_weight(isomap_rep[min_idx],isomap_rep[max_idx],0.25,isomap_rep, h =0.05))), col = col[2], lwd =2)
lines(ellipse(wasserstein_barycenter(cov_mat, find_weight(isomap_rep[min_idx],isomap_rep[max_idx],0.5,isomap_rep, h =0.05))), col = col[3], lwd =5)
lines(ellipse(wasserstein_barycenter(cov_mat, find_weight(isomap_rep[min_idx],isomap_rep[max_idx],0.75,isomap_rep, h =0.05))), col = col[4], lwd =2)
lines(ellipse(cov_mat[[max_idx]]), col = col[5], lwd =5)
legend("topleft", legend = c("t = 0","t = 0.25", "t = 0.5", "t = 0.75", "t = 1"),lty=1,lwd=c(5,2,5,2,5), col = col,cex=1.5)

