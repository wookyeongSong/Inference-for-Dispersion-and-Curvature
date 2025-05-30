# -------------------------------------------------------------------------
# Figure 3
#
# Manuscript reference : Section 5.1
# Simulation code reference : simulations/DC_distribution.R
# -------------------------------------------------------------------------

## Import main functions
source("simulations/DC_distribution.R")


## Generate Figure 3
pdf("Figure3.pdf", width = 20, height = 7)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings
layout(matrix(1:3, ncol = 3), widths = c(2,0.5,3))
colors <- colorRampPalette(c('red', 'blue'))(k)
idx_theta = order(theta)


### Figure 3a: Density 95% contour plots
par(mar = c(5, 5, 4, 1))  # Adjust margins for the main plot
plot(ellipse(cov_mat[[idx_theta[1]]],center = c(0,0)),xlim= c(-6,6), ylim = c(-6,6), type='l',col=colors[1],cex.lab=2,cex.axis=2,cex.main=3,lwd=1.5,xlab="", ylab = "",main="Density 95% contour plots")
for(i in 2:length(idx_theta)){
  obs = idx_theta[i]
  lines(ellipse(cov_mat[[obs]],center = c(0,0)),col=colors[i])
}

plot(1:10, (1:10)*10, type="n", bty="n", xaxt = "n", yaxt = "n", ann = FALSE) # Set up margins for the color bar and plot an empty plot area

image.plot(legend.only = TRUE, 
           zlim = c(1, 120),
           col = colorRampPalette(c('red', 'blue'))(100), 
           axis.args = list(at = seq(1, 120, length.out = 5), labels = c("0", "0.25", "0.5", "0.75", "1"), cex.axis = 2),
           legend.mar = 32,  
           legend.width = 8,  
           legend.args = list(text = expression(theta), side = 3, line = 1, cex = 2)) # Add a color bar legend using image.plot()




## Figure 3b: Intrinsic curvature test - Plot confidence region for the joint distribution of intrinsic metric and Frechet variance
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((3))
par(mar = c(5, 5, 4, 1))  # Adjust margins for the main plot
plot(ellipse(cov.est.wass,center = c(Vm.est.wass,Vf.est.wass),level=0.99),type='l',col=col[1],lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Proposed intrinsic curvature test")
points(Vm.est.wass,Vf.est.wass,pch=16)
lines(ellipse(cov.est.wass,center = c(Vm.est.wass,Vf.est.wass),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.wass,center = c(Vm.est.wass,Vf.est.wass),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

par(old_par)                          # restore original settings
dev.off()