# -------------------------------------------------------------------------
# Figure S.3
#
# Manuscript reference : Section S.1.3
# Simulation code reference : simulations/DC_power_analysis.R
# -------------------------------------------------------------------------

## Import main functions
source("simulations/DC_power_analysis.R")

## Generate Figure S.3
pdf("FigureS3.pdf", width = 12, height = 10)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings

rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal(3)
par(mfrow=c(1,1),mar = c(5,5,3,2))
plot(kappa_pos_param, res_power[1:11], ylim = c(-0.05, 1.05), lty = 1, type = 'l',col=col[1],cex.lab=2,cex.axis=2,cex.main=2,lwd=3,xlab=expression(paste("Alexandrov curvature ", kappa)) , ylab = "Power", main="Intrinsic Curvature Test")
lines(kappa_pos_param, res_power[12:22], col=col[2], lty = 2, lwd = 3)
lines(kappa_pos_param, res_power[23:33], col=col[3], lty = 3, lwd = 3)
legend("right", legend = c("n = 200", "n = 500", "n = 1000"),lty=c(1,2,3),lwd=3, col = col,cex=2)
abline(h = 0.05, lty = 2)
abline(h = 1, lty = 1)

par(old_par)                          # restore original settings
dev.off()