source("~/src/DC_mainfunctions.R")
source("~/DC_simulations_datagen.R")

############################################################
#### Simulations 1: Spherical data with geodesic metric ####
############################################################

## Number of samples 
k = 50

## Generate simulation data from open upper hemisphere with geodesic metric
df.sph = dat.gen.sph(k)

# Main function: Estimate the variances, curvature and test statistics when we have spherical data as inputs
sph.curv.fit = sph.curv.est(df.sph)

## Plot confidence region for the joint distribution of metric and Frechet variance when we have spherical data with geodesic distance
cov.est.sph = sph.curv.fit$cov.normalized
Vm.est.sph = sph.curv.fit$Vm
Vf.est.sph = sph.curv.fit$Vf

par(mfrow=c(1,1),mar = c(5,5,3,2))
plot(ellipse(cov.est.sph,center = c(Vm.est.sph,Vf.est.sph),level=0.99),type='l',xlim=c(0.7,1.5),ylim=c(0.7,1.5),col=col[1],lwd=2,xlab=expression(V[M]), ylab = expression(V[F]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Spherical Data (n = 50)")
points(Vm.est.sph,Vf.est.sph,pch=16)
lines(ellipse(cov.est.sph,center = c(Vm.est.sph,Vf.est.sph),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.sph,center = c(Vm.est.sph,Vf.est.sph),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=1)


