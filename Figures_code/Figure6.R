# -------------------------------------------------------------------------
# Figure 6
#
# Data Application: Gait synchronization analysis
# Data available at https://github.com/deepcharles/gait-data
# Reference: Truong et al. 2019, 'A data set for the study of human locomotion with inertial measurements units', Image Processing On Line 9, 381–390."
# Manuscript reference : Section 6.1
# Data application code reference : Applications/DC_gait.R
# -------------------------------------------------------------------------


## Import main functions
source("Applications/DC_gait.R")


## Generate Figure 6
pdf("Figure6.pdf", width = 25, height = 6)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings
par(mfrow = c(1, 4))                 # 1 row, 4 columns

rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((3))


### Figure 6(a) ─ Intrinsic curvature test for gait data from healthy group
par(mar=c(5,6,3,2))
plot(ellipse(cov.est.healthy,center = c(Vm.est.healthy,Vf.est.healthy),level=0.99), xlim=c(400, 1900), ylim = c(400,1900), type='l',col=col[1] ,lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Healthy Group")
points(Vm.est.healthy,Vf.est.healthy,pch=16)
lines(ellipse(cov.est.healthy,center = c(Vm.est.healthy,Vf.est.healthy),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.healthy,center = c(Vm.est.healthy,Vf.est.healthy),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)


### Figure 6(b) ─ Intrinsic curvature test for gait data from orthopedic group
par(mar=c(5,6,3,2))
plot(ellipse(cov.est.ortho,center = c(Vm.est.ortho,Vf.est.ortho),level=0.99), xlim=c(200, 1700), ylim = c(200,1700), type='l',col=col[1] ,lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Orthopedic Group")
points(Vm.est.ortho,Vf.est.ortho,pch=16)
lines(ellipse(cov.est.ortho,center = c(Vm.est.ortho,Vf.est.ortho),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.ortho,center = c(Vm.est.ortho,Vf.est.ortho),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)


### Figure 6(c) ─ Intrinsic Geodesic for gait data from healthy group
col = rgPal((5))

par(mar = c(5,6,3,2))
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


### Figure 6(d) ─ Intrinsic Geodesic for gait data from orthopedic group
par(mar = c(5,6,3,2))
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

par(old_par)                          # restore original settings
dev.off()

