source("~/src/DC_mainfunctions.R")
source("~/DC_simulations_datagen.R")

##################################################################################
#### Simulations 2: Symmetric Positive Definite Matrices with Various Metrics ####
##################################################################################

## Number of samples 
k = 100

## Generate simulation data from space of symmetric positive definite matrices
cov.array = dat.gen.spd(k)

# (a) Frobenius metric
cov.curv.fit.Frob = cov.curv.est(cov.array, metric = "Euclidean")

# (b) log-Euclidean metric
cov.curv.fit.LE = cov.curv.est(cov.array, metric = "LERM")

# (c) power Frobenius metric with power 1/2
cov.curv.fit.Root = cov.curv.est(cov.array, metric = "RootEuclidean")

# (d) Cholesky metric
cov.curv.fit.Chol = cov.curv.est(cov.array, metric = "Cholesky")

# (e) Affine-Invariant Riemannian (AIR) metric
cov.curv.fit.AIRM = cov.curv.est(cov.array, metric = "AIRM")

# (f) Procrustes size-and-shape (PSS) metric
cov.curv.fit.PSS = cov.curv.est(cov.array, metric = "Procrustes.SS")

## Ratio between metric variance over Frechet variance 
Vm.spd = c(cov.curv.fit.Frob$Vm, cov.curv.fit.LE$Vm, cov.curv.fit.Root$Vm, cov.curv.fit.Chol$Vm, cov.curv.fit.AIRM$Vm, cov.curv.fit.PSS$Vm)
Vf.spd = c(cov.curv.fit.Frob$Vf, cov.curv.fit.LE$Vf, cov.curv.fit.Root$Vf, cov.curv.fit.Chol$Vf, cov.curv.fit.AIRM$Vf, cov.curv.fit.PSS$Vf)

par(mar=c(5,5,3,2))
plot(Vf.spd/Vm.spd,col='black',xlab = "Metrics",ylab = "Ratio of Variances",xaxt='n',ylim=c(0.7,1.3),cex=2,cex.axis=2,cex.lab=2,pch=19)
axis(side = 1, at = 1:6, labels = c("Frob", "log-E", "p-Frob", "Chol","AIR", "PSS"), cex.axis=2)
legend("topleft",legend = c(expression(V[F]~"/"~V[M])),fill=c("black"),cex=2)

## Plot confidence region for the joint distribution of metric and Frechet variance with PSS metric and AIRM metric
cov.est.AIRM = cov.curv.fit.AIRM$cov.normalized
Vm.est.AIRM = cov.curv.fit.AIRM$Vm
Vf.est.AIRM = cov.curv.fit.AIRM$Vf

cov.est.PSS = cov.curv.fit.PSS$cov.normalized
Vm.est.PSS = cov.curv.fit.PSS$Vm
Vf.est.PSS = cov.curv.fit.PSS$Vf

rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((3))

par(mar=c(5,5,3,2))
plot(ellipse(cov.est.AIRM,center = c(Vm.est.AIRM,Vf.est.AIRM),level=0.99),type='l',col=col[1],main="Space of SPD matrix with AIR metric",cex.lab=2,cex.axis=2,cex.main=2.5,lwd=2,xlab=expression(V[M]), ylab = expression(V[F]),xlim=c(15,45),ylim=c(15,45))
points(Vm.est.AIRM,Vf.est.AIRM,pch=16)
lines(ellipse(cov.est.AIRM,center = c(Vm.est.AIRM,Vf.est.AIRM),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.AIRM,center = c(Vm.est.AIRM,Vf.est.AIRM),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

par(mar=c(5,5,3,2))
plot(ellipse(cov.est.PSS,center = c(Vm.est.PSS,Vf.est.PSS),level=0.99),type='l',col=col[1],main="Space of SPD matrix with PSS metric",cex.lab=2,cex.axis=2,cex.main=2.5,lwd=2,xlab=expression(V[M]), ylab = expression(V[F]),xlim=c(4500,9500),ylim=c(4500,9500))
points(Vm.est.PSS,Vf.est.PSS,pch=16)
lines(ellipse(cov.est.PSS,center = c(Vm.est.PSS,Vf.est.PSS),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.PSS,center = c(Vm.est.PSS,Vf.est.PSS),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

