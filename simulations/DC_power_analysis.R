source("~/src/DC_mainfunctions.R")
source("~/DC_simulations_datagen.R")

#########################
#### Power Analysis  ####
#########################
kappa_pos_param = c(0,seq(0.1, 1, by = 0.1))

k_vec = c(200, 500, 1000)
iter = 500

power_list = list()

for(k in k_vec){
  
  power_vec = c()
  for(kappa in kappa_pos_param){
    pval_vec = c()
    
    for(i in 1:iter){
      dat_gen = dat.gen.power.anal(k, kappa)
      df = dat_gen$df
      L = dat_gen$L
      curv.fit = intrinsic.curv.est(df, L = L)
      pval_vec = c(pval_vec, (curv.fit$pval < 0.05) )
      
    }
    
    power_vec = c(power_vec, mean(pval_vec))
    print(mean(pval_vec))
  }
  
  power_list = c(power_list,power_vec)
}

res_power = unlist(power_list)


t=2
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((t+1))
par(mfrow=c(1,1),mar = c(5,5,3,2))
plot(kappa_pos_param, res_power[1:11], ylim = c(-0.05, 1.05), lty = 1, type = 'l',col=col[1],cex.lab=2,cex.axis=2,cex.main=2,lwd=3,xlab=expression(paste("Alexandrov curvature ", kappa)) , ylab = "Power", main="Intrinsic Curvature Test")
lines(kappa_pos_param, res_power[12:22], col=col[2], lty = 2, lwd = 3)
lines(kappa_pos_param, res_power[23:33], col=col[3], lty = 3, lwd = 3)
legend("right", legend = c("n = 200", "n = 500", "n = 1000"),lty=c(1,2,3),lwd=3, col = col,cex=2)
abline(h = 0.05, lty = 2)
abline(h = 1, lty = 1)
