# -------------------------------------------------------------------------
# Figure S.6
#
# Manuscript reference : Section S.5
# Simulation code reference : simulations/DC_distribution.R
# -------------------------------------------------------------------------

## Import main functions
source("simulations/DC_distribution.R")

## Generate Figure S.6
pdf("FigureS6.pdf", width = 8, height = 7)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings

### Figure 4(a) ─ ISOMAP Representation
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((5))
t1 = 0.25
t2 = 0.5
t3 = 0.75
par(mfrow = c(1, 1), mar = c(5,5,3,2))
plot(c(theta)/(pi/2), isomap_rep[1:100], xlim = c(0,1), ylim = c(-1.2, 1.2), xlab = expression(theta), ylab = "ISOMAP Representation",cex.lab=2,cex.axis=2,cex.main=2.5)
points(0,isomap_rep[101], col = col[1], pch = 19, cex = 3)
points(find_weight(isomap_rep[min_idx],isomap_rep[max_idx],t1,isomap_rep, h =0.01) %*% c(theta/(pi/2),0,1), (1-t1)*isomap_rep[101] + t1*isomap_rep[102], col = col[2], pch = 19, cex = 3 )
points(find_weight(isomap_rep[min_idx],isomap_rep[max_idx],t2,isomap_rep, h =0.01) %*% c(theta/(pi/2),0,1), (1-t2)*isomap_rep[101] + t2*isomap_rep[102], col = col[3], pch = 19, cex = 3 )
points(find_weight(isomap_rep[min_idx],isomap_rep[max_idx],t3,isomap_rep, h =0.01) %*% c(theta/(pi/2),0,1), (1-t3)*isomap_rep[101] + t3*isomap_rep[102], col = col[4], pch = 19, cex = 3 )
points(1,isomap_rep[102], col = col[5], pch = 19, cex = 3)
legend("topleft", legend = c("t = 0","t = 0.25", "t = 0.5", "t = 0.75", "t = 1"),pch = 19, col = col,cex=1.5)

par(old_par)                          # restore original settings
dev.off()