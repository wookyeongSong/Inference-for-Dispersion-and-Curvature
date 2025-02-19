source("~/Desktop/WK_UCD/Research/varianceproj/src/DC_mainfunctions.R")
source("~/Desktop/WK_UCD/Research/varianceproj/DC_simulations_datagen.R")

#########################################
#### Simulations 1: Point Cloud Data ####
#########################################

#########################################
# B1: Positively curved manifold
## Number of samples 
k = 1000

## Generate simulation data from point cloud on positively curved manifold
df.pc.pos = dat.gen.pc.pos(k)

## Estimate intrinsic curvature and test statistics (main function)
curv.fit.pc.pos = intrinsic.curv.est(df.pc.pos, L = 4)

## Confidence intervals for the intrinsic curvatures rho_{I}
# alpha = 0.01 (0.072, 0.132)
C.I.lower.pc.pos.01 = curv.fit.pc.pos$rho - qnorm(0.995)*curv.fit.pc.pos$sd/sqrt(k)
C.I.upper.pc.pos.01 = curv.fit.pc.pos$rho + qnorm(0.995)*curv.fit.pc.pos$sd/sqrt(k)

# alpha = 0.05 (0.080, 0.125)
C.I.lower.pc.pos.05 = curv.fit.pc.pos$rho - qnorm(0.975)*curv.fit.pc.pos$sd/sqrt(k)
C.I.upper.pc.pos.05 = curv.fit.pc.pos$rho + qnorm(0.975)*curv.fit.pc.pos$sd/sqrt(k)

# alpha = 0.1 (0.083, 0.122)
C.I.lower.pc.pos.1 = curv.fit.pc.pos$rho - qnorm(0.95)*curv.fit.pc.pos$sd/sqrt(k)
C.I.upper.pc.pos.1 = curv.fit.pc.pos$rho + qnorm(0.95)*curv.fit.pc.pos$sd/sqrt(k)


## Plot confidence region for the joint distribution of intrinsic metric and Frechet variance
cov.est.pc.pos = curv.fit.pc.pos$cov.normalized
Vm.est.pc.pos = curv.fit.pc.pos$Vm
Vf.est.pc.pos = curv.fit.pc.pos$Vf

rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((3))

par(mar=c(5,5,3,2))
plot(ellipse(cov.est.pc.pos,center = c(Vm.est.pc.pos,Vf.est.pc.pos),level=0.99),type='l',col=col[1],xlim = c(1.6,2.6), ylim =c(1.6,2.6) ,lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Point cloud data on spherical manifold")
points(Vm.est.pc.pos,Vf.est.pc.pos,pch=16)
lines(ellipse(cov.est.pc.pos,center = c(Vm.est.pc.pos,Vf.est.pc.pos),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.pc.pos,center = c(Vm.est.pc.pos,Vf.est.pc.pos),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

#########################################
# B2: Negatively curved manifold

## Number of samples 
k = 1000

## Generate simulation data from point cloud on positively curved manifold
df.pc.neg = dat.gen.pc.neg(k)

## Estimate intrinsic curvature and test statistics (main function)
curv.fit.pc.neg = intrinsic.curv.est(df.pc.neg, L = 10)

## Confidence intervals for the intrinsic curvatures rho_{I}
# alpha = 0.01 (-0.098, -0.072)
C.I.lower.pc.neg.01 = curv.fit.pc.neg$rho - qnorm(0.995)*curv.fit.pc.neg$sd/sqrt(k)
C.I.upper.pc.neg.01 = curv.fit.pc.neg$rho + qnorm(0.995)*curv.fit.pc.neg$sd/sqrt(k)

# alpha = 0.05 (-0.095, -0.075)
C.I.lower.pc.neg.05 = curv.fit.pc.neg$rho - qnorm(0.975)*curv.fit.pc.neg$sd/sqrt(k)
C.I.upper.pc.neg.05 = curv.fit.pc.neg$rho + qnorm(0.975)*curv.fit.pc.neg$sd/sqrt(k)

# alpha = 0.1 (-0.093, -0.077)
C.I.lower.pc.neg.1 = curv.fit.pc.neg$rho - qnorm(0.95)*curv.fit.pc.neg$sd/sqrt(k)
C.I.upper.pc.neg.1 = curv.fit.pc.neg$rho + qnorm(0.95)*curv.fit.pc.neg$sd/sqrt(k)

## Plot confidence region for the joint distribution of intrinsic metric and Frechet variance
cov.est.pc.neg = curv.fit.pc.neg$cov.normalized
Vm.est.pc.neg = curv.fit.pc.neg$Vm
Vf.est.pc.neg = curv.fit.pc.neg$Vf

rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((3))

par(mar=c(5,5,3,2))
plot(ellipse(cov.est.pc.neg,center = c(Vm.est.pc.neg,Vf.est.pc.neg),level=0.99),type='l',col=col[1],xlim = c(5,8), ylim =c(5,8),lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Point cloud data on hyperbolic manifold")
points(Vm.est.pc.neg,Vf.est.pc.neg,pch=16)
lines(ellipse(cov.est.pc.neg,center = c(Vm.est.pc.neg,Vf.est.pc.neg),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.pc.neg,center = c(Vm.est.pc.neg,Vf.est.pc.neg),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

#########################################
# B3: Flat manifold

## Number of samples 
k = 1000

## Generate simulation data from point cloud on positively curved manifold
df.pc.flat = dat.gen.pc.flat(k)

## Estimate intrinsic curvature and test statistics (main function)
curv.fit.pc.flat = intrinsic.curv.est(df.pc.flat, L = 4)

## Confidence intervals for the intrinsic curvatures rho_{I}
# alpha = 0.01 (-0.023, 0.051)
C.I.lower.pc.flat.01 = curv.fit.pc.flat$rho - qnorm(0.995)*curv.fit.pc.flat$sd/sqrt(k)
C.I.upper.pc.flat.01 = curv.fit.pc.flat$rho + qnorm(0.995)*curv.fit.pc.flat$sd/sqrt(k)

# alpha = 0.05 (-0.014, 0.042)
C.I.lower.pc.flat.05 = curv.fit.pc.flat$rho - qnorm(0.975)*curv.fit.pc.flat$sd/sqrt(k)
C.I.upper.pc.flat.05 = curv.fit.pc.flat$rho + qnorm(0.975)*curv.fit.pc.flat$sd/sqrt(k)

# alpha = 0.1 (-0.009, 0.038)
C.I.lower.pc.flat.1 = curv.fit.pc.flat$rho - qnorm(0.95)*curv.fit.pc.flat$sd/sqrt(k)
C.I.upper.pc.flat.1 = curv.fit.pc.flat$rho + qnorm(0.95)*curv.fit.pc.flat$sd/sqrt(k)

## Plot confidence region for the joint distribution of intrinsic metric and Frechet variance
cov.est.pc.flat = curv.fit.pc.flat$cov.normalized
Vm.est.pc.flat = curv.fit.pc.flat$Vm
Vf.est.pc.flat = curv.fit.pc.flat$Vf

rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((3))

par(mar=c(5,5,3,2))
plot(ellipse(cov.est.pc.flat,center = c(Vm.est.pc.flat,Vf.est.pc.flat),level=0.99),type='l',col=col[1],xlim = c(0.4,0.65), ylim =c(0.4,0.65),lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Point cloud data on flat manifold")
points(Vm.est.pc.flat,Vf.est.pc.flat,pch=16)
lines(ellipse(cov.est.pc.flat,center = c(Vm.est.pc.flat,Vf.est.pc.flat),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.pc.flat,center = c(Vm.est.pc.flat,Vf.est.pc.flat),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)