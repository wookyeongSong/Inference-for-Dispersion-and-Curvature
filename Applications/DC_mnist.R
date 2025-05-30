# -------------------------------------------------------------------------
# Additional real data analysis: MNIST Dataset
#
# Data available at R package `dslabs`
# Manuscript reference : Section S.2.1
# Table reproduced     : Table S.1
# -------------------------------------------------------------------------

## Import main functions
source("src/DC_mainfunctions.R")


## Import MNIST dataset
mnist <- read_mnist(download = TRUE, destdir = "~/Downloads")


## Image Data Preprocessing
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


## Intrinsic curvature test
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

}


## Table S.1
label_freq <- as.data.frame(table(mnist$train$labels))
colnames(label_freq) <- c("Label", "Count")   # clearer names

tbl_mnist <- data.frame(
  Label      = label_freq$Label,
  Num_sample      = label_freq$Count,
  Lower_CI_95 = lower_CI_train_image,
  Upper_CI_95 = upper_CI_train_image,
  Curvature_Estimate   = rho_train_image,
  stringsAsFactors = FALSE
)

tbl_mnist <- transform(tbl_mnist,
                       Lower_CI_95 = round(Lower_CI_95, 3),
                       Upper_CI_95 = round(Upper_CI_95, 3),
                       Curvature_Estimate = round(Curvature_Estimate,3))

print("MNIST Dataset: Confidence interval for intrinsic curvature")
print(tbl_mnist, row.names = FALSE)

