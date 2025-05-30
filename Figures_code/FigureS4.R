# -------------------------------------------------------------------------
# Figure S.4
#
# Manuscript reference : Section S.1.4
# Simulation code reference : simulations/DC_highdim.R
# -------------------------------------------------------------------------


## Import main functions
source("simulations/DC_highdim.R")


## Generate Figure S.4
pdf("FigureS4.pdf", width = 12, height = 12)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings
par(mfrow = c(2, 2))                 # 2 row, 2 columns

### Figure S.4(a): Power analysis for increasing ambient space dimension with low noise sigma = 0.05
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal(3)
par(mar = c(5,5,3,2))
plot(dim_param, hdim_power_list[[1]][1:11],xlim = c(4, 100), ylim = c(0, 1.05), lty = 1, type = 'l',col=col[1],cex.lab=2,cex.axis=2,cex.main=2,lwd=3,xlab="Dimension p" , ylab = "Power", main=expression(paste("Fixed Noise ", sigma, " = 0.05 (low)")))
lines(dim_param, hdim_power_list[[1]][12:22], col=col[2], lty = 2, lwd = 3)
lines(dim_param, hdim_power_list[[1]][23:33], col=col[3], lty = 3, lwd = 3)
legend("bottomright", legend = c("n = 200", "n = 500", "n = 1000"),lty=c(1,2,3),lwd=3, col = col,cex=2)
abline(h = 0.05, lty = 2)
abline(h = 1, lty = 1)

### Figure S.4(b): Power analysis for increasing ambient space dimension with middle noise sigma = 0.1
par(mar = c(5,5,3,2))
plot(dim_param, hdim_power_list[[2]][1:11],xlim = c(4, 100), ylim = c(0, 1.05), lty = 1, type = 'l',col=col[1],cex.lab=2,cex.axis=2,cex.main=2,lwd=3,xlab="Dimension p" , ylab = "Power", main=expression(paste("Fixed Noise ", sigma, " = 0.1 (middle)")))
lines(dim_param, hdim_power_list[[2]][12:22], col=col[2], lty = 2, lwd = 3)
lines(dim_param, hdim_power_list[[2]][23:33], col=col[3], lty = 3, lwd = 3)
legend("bottomright", legend = c("n = 200", "n = 500", "n = 1000"),lty=c(1,2,3),lwd=3, col = col,cex=2)
abline(h = 0.05, lty = 2)
abline(h = 1, lty = 1)

### Figure S.4(c): Power analysis for increasing ambient space dimension with high noise sigma = 0.2
par(mar = c(5,5,3,2))
plot(dim_param, hdim_power_list[[3]][1:11],xlim = c(4, 100), ylim = c(0, 1.05), lty = 1, type = 'l',col=col[1],cex.lab=2,cex.axis=2,cex.main=2,lwd=3,xlab="Dimension p" , ylab = "Power", main=expression(paste("Fixed Noise ", sigma, " = 0.2 (high)")))
lines(dim_param, hdim_power_list[[3]][12:22], col=col[2], lty = 2, lwd = 3)
lines(dim_param, hdim_power_list[[3]][23:33], col=col[3], lty = 3, lwd = 3)
legend("bottomright", legend = c("n = 200", "n = 500", "n = 1000"),lty=c(1,2,3),lwd=3, col = col,cex=2)
abline(h = 0.05, lty = 2)
abline(h = 1, lty = 1)

### Figure S.4(d): Power analysis for increasing ambient space dimension with fixed signal-to-noise ratio
par(mar = c(5,5,3,2))
plot(dim_param, hdim_power_list2[[1]][1:11],xlim = c(4, 100), ylim = c(0, 1.05), lty = 1, type = 'l',col=col[1],cex.lab=2,cex.axis=2,cex.main=2,lwd=3,xlab="Dimension p" , ylab = "Power", main=expression(paste("Fixed Signal-to-Noise-Ratio ", sigma, " = 3/(10",sqrt(p),")")))
lines(dim_param, hdim_power_list2[[1]][12:22], col=col[2], lty = 2, lwd = 3)
lines(dim_param, hdim_power_list2[[1]][23:33], col=col[3], lty = 3, lwd = 3)
legend("bottomright", legend = c("n = 200", "n = 500", "n = 1000"),lty=c(1,2,3),lwd=3, col = col,cex=2)
abline(h = 0.05, lty = 2)
abline(h = 1, lty = 1)

par(old_par)                          # restore original settings
dev.off()