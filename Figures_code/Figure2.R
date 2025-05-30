# -------------------------------------------------------------------------
# Figure 2
#
# Manuscript reference : Section 4.2
# Simulation code reference : simulations/DC_sphere_contam.R
# -------------------------------------------------------------------------


## Import main functions
source("simulations/DC_sphere_contam.R")

## Generate Figure 2
pdf("Figure2.pdf", width = 20, height = 6)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings
par(mfrow = c(1, 4))                 # 1 row, 4 columns
col = rainbow_hcl(3)                 # set color


### Figure 2(a) ─ Random samples with noise sigma = 1/32
par(mar = c(1, 1, 2, 1))             
scatterplot3d::scatterplot3d(
  sph.list[[1]], pch = 19, color = col[1],
  xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), zlim = c(-1.5, 3),cex.main=2.5,cex.lab=2,cex.axis=2,
  xlab = expression(~~~~sigma == 1/32), ylab = "", zlab = "",
  box = FALSE, x.ticklabs = "", y.ticklabs = "", z.ticklabs = ""
)

### Figure 2(b) ─ Random samples with noise sigma = 1/8
par(mar = c(1, 1, 2, 1))
scatterplot3d::scatterplot3d(
  sph.list[[3]], pch = 19, color = col[2],
  xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), zlim = c(-1.5, 3),cex.main=2.5,cex.lab=2,cex.axis=2,
  xlab = expression(~~~~sigma == 1/8), ylab = "", zlab = "",
  box = FALSE, x.ticklabs = "", y.ticklabs = "", z.ticklabs = ""
)

### Figure 2(c) ─ Random samples with noise sigma = 1/2
par(mar = c(1, 1, 2, 1))
scatterplot3d::scatterplot3d(
  sph.list[[5]], pch = 19, color = col[3],
  xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), zlim = c(-1.5, 3),cex.main=2.5,cex.lab=2,cex.axis=2,
  xlab = expression(~~~~sigma == 1/2), ylab = "", zlab = "",
  box = FALSE, x.ticklabs = "", y.ticklabs = "", z.ticklabs = ""
)

### Figure 2(d) ─ Confidence interval for the intrinsic curvature
par(mar = c(5, 5, 12, 7))             # bigger margins for axes
plotCI(
  x = 1:5, y = rho.vec,
  uiw = upper.vec - rho.vec, liw = upper.vec - rho.vec,
  ylab = expression(rho[I]), xlab = expression(sigma),
  ylim = c(-0.02, 0.08), axes = FALSE,
  scol = "black", col = "black", pch = 19, font.axis = 2,
  cex.axis = 2, cex.lab = 2, cex.main = 2.5,
  main = "Intrinsic metric curvature"
)
lines(1:5, rho.vec, lwd = 1)
abline(h = 0, lty = 2)
axis(2, cex.axis = 2)
axis(1, at = 1:5, labels = c("1/32", "1/16", "1/8", "1/4", "1/2"),
     font.axis = 1, cex.axis = 2)
box(bty = "l")

par(old_par)                          # restore original settings
dev.off()
