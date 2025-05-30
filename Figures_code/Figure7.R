# -------------------------------------------------------------------------
# Figure 7
#
# Data Application: Energy source data analysis
# Data available at Applications/energy_data
# Manuscript reference : Section 6.2
# Data application code reference : Applications/DC_energy.R
# -------------------------------------------------------------------------


## Import main functions
source("Applications/DC_energy.R")

## Generate Figure 7
pdf("Figure7.pdf", width = 26, height = 7)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings
par(mfrow = c(1, 3),                     # ← 1 row × 3 columns
    mar   = c(5, 6, 3, 2),               # inner margins for each cell
    oma   = c(0, 0, 2, 0))               # outer top space for a title (optional)


## Figure 7(a): Intrinsic curvature test using confidence region for the joint distribution of intrinsic metric and Frechet variance
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal(3)
plot(ellipse(cov.est.energy,center = c(Vm.est.energy,Vf.est.energy),level=0.99),type='l',col=col[1],lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="")
points(Vm.est.energy,Vf.est.energy,pch=16)
lines(ellipse(cov.est.energy,center = c(Vm.est.energy,Vf.est.energy),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.energy,center = c(Vm.est.energy,Vf.est.energy),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)


## Figure 7(b): Energy proportional data trajectory (Estimated geodesics)
time.stamp = seq(0,1,0.1)

geod.sph = data.frame()
for (tt in time.stamp){
  
  temp.geo = find_weight(min(isomap_rep),max(isomap_rep),tt,isomap_rep, h =0.01)
  geod.sph = rbind(geod.sph, find_weighted_barycenter_on_sphere(comp.data, temp.geo))
  
}
colnames(geod.sph) = c("CoalPet","NaturalGas","Other")

colors <- colorRampPalette(c("cadetblue", "coral3"))(nrow(comp.data))
s3d <- scatterplot3d(comp.data[, 1], comp.data[, 2], comp.data[, 3], color = colors, pch = 16,
                     main = "Proportion of energy source",
                     xlab = "sqrt( Coal + Pet )", ylab = "sqrt( Natural Gas )", zlab = "sqrt( Others )",
                     cex.main = 3, cex.axis = 2, cex.lab = 2, box = FALSE)

for (i in 1:(nrow(geod.sph) - 1)) {
  s3d$points3d(x = c(geod.sph[i, 1], geod.sph[i + 1, 1]),
               y = c(geod.sph[i, 2], geod.sph[i + 1, 2]),
               z = c(geod.sph[i, 3], geod.sph[i + 1, 3]),
               type = "l", col = "grey30", lty = 1,lwd = 3)  # Line color set to black
}



## Figure 7(c): ISOMAP representation interpolation
years <- seq(2012, 2021, length.out = 10)
par(fig = c(2 / 3, 1, 0, 1), new = TRUE, mar = c(5, 5, 3, 9))

plot(1:120, isomap_rep, col = colors, pch = 16,
     xlab = "Year", ylab = "ISOMAP Representation",
     cex.main = 2, cex.axis = 2, cex.lab = 2, main = "", xaxt = "n")
axis(1, at = seq(1, 120, length.out = 4), labels = as.character(years)[c(1,4,7,10)], cex.axis = 2)

par(new = TRUE, fig = c(0.94, 0.99, 0.1, 0.9), mar = c(0, 0, 0, 4))
plot(1:10, (1:10)*10, type="n", bty="n", xaxt = "n", yaxt = "n", ann = FALSE)

years <- seq(2012, 2021, length.out = 10)
image.plot(legend.only = TRUE, 
           zlim = c(1, 120),  # Scale according to the data size
           col = colorRampPalette(c("cadetblue", "coral3"))(100), 
           axis.args = list(at = seq(1, 120, length.out = 10), labels = as.character(years), cex.axis = 2),
           legend.mar = 6,  # Increased space on the right for the color bar
           legend.width = 8,  # Adjust the width of the legend bar
           legend.args = list(text = "Year", side = 3, line = 1, cex = 2))
par(fig = c(0, 1, 0, 1), new = FALSE)   # reset fig

par(old_par)                          # restore original settings
dev.off()
