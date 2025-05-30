# -------------------------------------------------------------------------
# Figure 4
#
# Manuscript reference : Section 5.1
# Simulation code reference : simulations/DC_distribution.R
# -------------------------------------------------------------------------


## Import main functions
source("simulations/DC_distribution.R")

## Generate Figure 4
pdf("Figure4.pdf", width = 20, height = 7)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings
par(mfrow = c(1, 3))                 # 1 row, 3 columns

### Figure 4(a) ─ ISOMAP Representation
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((5))
t1 = 0.25
t2 = 0.5
t3 = 0.75
par(mar = c(5,5,3,2))
plot(c(theta)/(pi/2), isomap_rep[1:100], xlim = c(0,1), ylim = c(-1.2, 1.2), xlab = expression(theta), ylab = "ISOMAP Representation",cex.lab=2,cex.axis=2,cex.main=2.5)
points(0,isomap_rep[101], col = col[1], pch = 19, cex = 3)
points(find_weight(isomap_rep[min_idx],isomap_rep[max_idx],t1,isomap_rep, h =0.01) %*% c(theta/(pi/2),0,1), (1-t1)*isomap_rep[101] + t1*isomap_rep[102], col = col[2], pch = 19, cex = 3 )
points(find_weight(isomap_rep[min_idx],isomap_rep[max_idx],t2,isomap_rep, h =0.01) %*% c(theta/(pi/2),0,1), (1-t2)*isomap_rep[101] + t2*isomap_rep[102], col = col[3], pch = 19, cex = 3 )
points(find_weight(isomap_rep[min_idx],isomap_rep[max_idx],t3,isomap_rep, h =0.01) %*% c(theta/(pi/2),0,1), (1-t3)*isomap_rep[101] + t3*isomap_rep[102], col = col[4], pch = 19, cex = 3 )
points(1,isomap_rep[102], col = col[5], pch = 19, cex = 3)
legend("topleft", legend = c("t = 0","t = 0.25", "t = 0.5", "t = 0.75", "t = 1"),pch = 19, col = col,cex=1.5)

### Figure 4(b) ─ Estimated Intrinsic Geodesic
par(mar = c(5,5,3,2))
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

### Figure 4(c) ─ Estimated Wasserstein Geodesic
par(mar = c(5,5,3,2))
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

par(old_par)                          # restore original settings
dev.off()