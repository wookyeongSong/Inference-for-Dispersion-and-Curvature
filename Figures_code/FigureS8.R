# -------------------------------------------------------------------------
# Sensitivity Analysis for the choice of input distance using distributional simulations with Wasserstein distance
#
# Manuscript reference : Section S.6
# Figure reproduced    : Figure S.8
# Simulation code reference : simulations/DC_distribution.R
# -------------------------------------------------------------------------


## Import main functions
source("simulations/DC_distribution.R")

## Generate Figure S.8
pdf("FigureS8.pdf", width = 14, height = 8)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings

rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal(3)
par(mfrow = c(1,2), mar=c(5,5,3,2),oma = c(0, 1, 5, 0))


## Left panel of Figure S.8: Intrinsic curvature test using 2-Wasserstein metric
plot(ellipse(cov.est.wass,center = c(Vm.est.wass,Vf.est.wass),level=0.99),type='l',col=col[1],lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2,cex.lab=2,cex.axis=2,main="Input distance: Wasserstein Metric")
points(Vm.est.wass,Vf.est.wass,pch=16)
lines(ellipse(cov.est.wass,center = c(Vm.est.wass,Vf.est.wass),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.wass,center = c(Vm.est.wass,Vf.est.wass),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)


## Right panel of Figure S.8: Intrinsic curvature test using Fisher-Rao metric
plot(ellipse(cov.est.wass.FR,center = c(Vm.est.wass.FR,Vf.est.wass.FR),level=0.99),type='l',col=col[1],lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2,cex.lab=2,cex.axis=2,main="Input distance: Fisher-Rao Metric")
points(Vm.est.wass.FR,Vf.est.wass.FR,pch=16)
lines(ellipse(cov.est.wass.FR,center = c(Vm.est.wass.FR,Vf.est.wass.FR),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.wass.FR,center = c(Vm.est.wass.FR,Vf.est.wass.FR),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

mtext("Distributional Data Simulation", side = 3, font = 2, outer = TRUE, line = 1, cex = 2.5)

par(old_par)                          # restore original settings
dev.off()

