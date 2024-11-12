source("~/src/DC_mainfunctions.R")

library(dslabs)

mnist <- read_mnist(download = TRUE, destdir = "~/Downloads")

rho_test_image = c()
pval_test_image = c()
labels = 0:9

for (label in labels){
  
  dat_idx = which(mnist$test$labels == label)
  
  df_image = mnist$test$images[dat_idx,]
  
  res_image = intrinsic.curv.est(df_image/255,L=4)
  
  pval_test_image = c(pval_test_image, res_image$pval)
  rho_test_image = c(rho_test_image, res_image$rho)
  
}

table(mnist$train$labels)

rho_train_image = c()
pval_train_image = c()
upper_CI_train_image = c()
lower_CI_train_image = c()

for (label in labels){
  
  dat_idx = which(mnist$train$labels == label)
  
  df_image = mnist$train$images[dat_idx,]
  k=nrow(df_image)
  res_image = intrinsic.curv.est(df_image/255,L=4)
  
  pval_train_image = c(pval_train_image, res_image$pval)
  rho_train_image = c(rho_train_image, res_image$rho)
  
  upper_CI_train_image = c(upper_CI_train_image, res_image$rho + qnorm(0.975)*res_image$sd/sqrt(k))
  lower_CI_train_image = c(lower_CI_train_image, res_image$rho - qnorm(0.975)*res_image$sd/sqrt(k))
  print(label)
  
}

rho_train_image
pval_train_image
upper_CI_train_image
lower_CI_train_image
