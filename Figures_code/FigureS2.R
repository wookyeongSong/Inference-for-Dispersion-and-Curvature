# -------------------------------------------------------------------------
# Figure S.2
#
# Manuscript reference : Section S.1.2
# Simulation code reference : simulations/DC_sphere.R
# -------------------------------------------------------------------------


## Import main functions
source("simulations/DC_sphere.R")


## Generate Figure S.2
pdf("FigureS2.pdf", width = 6, height = 6)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings

par(mfrow=c(1,1),mar = c(5,5,3,2))
plot(ellipse(cov.est.sph,center = c(Vm.est.sph,Vf.est.sph),level=0.99),type='l',xlim=c(0.7,1.5),ylim=c(0.7,1.5),col=col[1],lwd=2,xlab=expression(V[M]), ylab = expression(V[F]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Spherical Data (n = 50)")
points(Vm.est.sph,Vf.est.sph,pch=16)
lines(ellipse(cov.est.sph,center = c(Vm.est.sph,Vf.est.sph),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.sph,center = c(Vm.est.sph,Vf.est.sph),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=1)

par(old_par)                          # restore original settings
dev.off()