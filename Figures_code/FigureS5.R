# -------------------------------------------------------------------------
# Figure S.5
#
# Data available at Applications/airport data
# Manuscript reference : Section S.2.2
# Data application code reference : Applications/DC_temperature.R
# -------------------------------------------------------------------------


## Import main functions
source("src/DC_temperature.R")


## Generate Figure S.5
pdf("FigureS5.pdf", width = 24, height = 8)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings
par(mfrow = c(1, 2))                 # 2 row, 2 columns

col = rainbow_hcl(3)

par(mar=c(5,5,3,2))
plotCI(x= 1:8, y= data.summer$value, uiw = 1.96*data.summer$se, liw = 1.96*data.summer$se,ylab = expression(V[M]),xlab = "Airport",axes=FALSE,scol="black",col=c(rep(col[1],4),rep(col[3],4)),cex=2,pch=19,font.axis=2, cex.axis=2, cex.lab=2, cex.main = 2.5,main="Metric variance of summer temperature distributions")
axis(side=2,cex.axis=2)         ## add default y-axis (ticks+labels)
axis(side=1,at=1:8,  ## add custom x-axis
     label=names.abb.summer,font.axis=1, cex.axis=2, cex.lab=2)
box(bty="l")    

par(mar=c(5,5,3,2))
plotCI(x=1:8, y= data.winter$value, uiw = 1.96*data.winter$se, liw = 1.96*data.winter$se,ylab = expression(V[M]),xlab = "Airport",axes=FALSE,scol="black",col=c(rep(col[1],4),rep(col[3],4)),cex=2,pch=19,font.axis=2, cex.axis=2, cex.lab=2, cex.main = 2.5,main="Metric variance of winter temperature distributions")
axis(side=2,cex.axis=2)         ## add default y-axis (ticks+labels)
axis(side=1,at=1:8,  ## add custom x-axis
     label=names.abb.winter,font.axis=1, cex.axis=2, cex.lab=2)
box(bty="l")

par(old_par)                          # restore original settings
dev.off()
