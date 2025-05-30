# -------------------------------------------------------------------------
# Power and Type-I error analysis
#
# Manuscript reference : Section S.1.3
# Figure reproduced    : Figure S.3
# -------------------------------------------------------------------------

## Import main functions
source("src/DC_mainfunctions.R")
source("simulations/DC_simulations_datagen.R")


## Power Analysis 
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
    
  }
  
  power_list = c(power_list,power_vec)
}

res_power = unlist(power_list)