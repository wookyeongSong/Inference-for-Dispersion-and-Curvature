source("~/Desktop/WK_UCD/Research/varianceproj/src/DC_mainfunctions.R")
source("~/Desktop/WK_UCD/Research/varianceproj/DC_simulations_datagen.R")


####################################################################################
#### High Dimensional Euclidean data with low-dimensional intrinsic dimension  #####
####################################################################################

#########################################################
#### Scenario 1: Low, middle, high noise simulation #####
#########################################################

dim_param    <- c(3, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
noise_param  <- c(0.05, 0.1, 0.2)
k_vec        <- c(200, 500, 1000)
iter         <- 500

# Pre-allocate the list to the correct length
hdim_power_list <- vector("list", length(noise_param))

for (n_idx in seq_along(noise_param)) {
  noise <- noise_param[n_idx]
  
  # We will collect all results for this noise in a numeric vector;
  # length = length(dim_param) * length(k_vec)
  power_hdim <- numeric(length(dim_param) * length(k_vec))
  
  # We will fill 'power_hdim' in chunks corresponding to each k
  start_idx  <- 1
  
  for (k in k_vec) {
    
    # Use sapply to compute power for each dim_p
    power_vec <- sapply(dim_param, function(dim_p) {
      
      # replicate(...) helps avoid manual loops for iter times
      pval_vec <- replicate(iter, {
        df       <- dat.gen.high.dim(k, dim_p, noise)
        curv_fit <- intrinsic.curv.est(df, L = 8)
        as.numeric(curv_fit$pval < 0.05)
      })
      
      mean(pval_vec)  # average of 0/1 to get the power
    })
    
    # Insert the results into the correct location of power_hdim
    end_idx <- start_idx + length(dim_param) - 1
    power_hdim[start_idx:end_idx] <- power_vec
    
    start_idx <- end_idx + 1
  }
  
  hdim_power_list[[n_idx]] <- power_hdim
}

t=2
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((t+1))
par(mfrow=c(1,1),mar = c(5,5,3,2))
plot(dim_param, hdim_power_list[[1]][1:11],xlim = c(4, 100), ylim = c(0, 1.05), lty = 1, type = 'l',col=col[1],cex.lab=2,cex.axis=2,cex.main=2,lwd=3,xlab="Dimension p" , ylab = "Power", main=expression(paste("Fixed Noise ", sigma, " = 0.05 (low)")))
lines(dim_param, hdim_power_list[[1]][12:22], col=col[2], lty = 2, lwd = 3)
lines(dim_param, hdim_power_list[[1]][23:33], col=col[3], lty = 3, lwd = 3)
legend("bottomright", legend = c("n = 200", "n = 500", "n = 1000"),lty=c(1,2,3),lwd=3, col = col,cex=2)
abline(h = 0.05, lty = 2)
abline(h = 1, lty = 1)


plot(dim_param, hdim_power_list[[2]][1:11],xlim = c(4, 100), ylim = c(0, 1.05), lty = 1, type = 'l',col=col[1],cex.lab=2,cex.axis=2,cex.main=2,lwd=3,xlab="Dimension p" , ylab = "Power", main=expression(paste("Fixed Noise ", sigma, " = 0.1 (middle)")))
lines(dim_param, hdim_power_list[[2]][12:22], col=col[2], lty = 2, lwd = 3)
lines(dim_param, hdim_power_list[[2]][23:33], col=col[3], lty = 3, lwd = 3)
legend("bottomright", legend = c("n = 200", "n = 500", "n = 1000"),lty=c(1,2,3),lwd=3, col = col,cex=2)
abline(h = 0.05, lty = 2)
abline(h = 1, lty = 1)


plot(dim_param, hdim_power_list[[3]][1:11],xlim = c(4, 100), ylim = c(0, 1.05), lty = 1, type = 'l',col=col[1],cex.lab=2,cex.axis=2,cex.main=2,lwd=3,xlab="Dimension p" , ylab = "Power", main=expression(paste("Fixed Noise ", sigma, " = 0.2 (high)")))
lines(dim_param, hdim_power_list[[3]][12:22], col=col[2], lty = 2, lwd = 3)
lines(dim_param, hdim_power_list[[3]][23:33], col=col[3], lty = 3, lwd = 3)
legend("bottomright", legend = c("n = 200", "n = 500", "n = 1000"),lty=c(1,2,3),lwd=3, col = col,cex=2)
abline(h = 0.05, lty = 2)
abline(h = 1, lty = 1)

##################################################
#### Scenario 2: Fixed Signal-to-Noise Ratio #####
##################################################

dim_param    <- c(3, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
noise_param  <- 0.3
k_vec        <- c(200, 500, 1000)
iter         <- 100

# Pre-allocate the list to the correct length
hdim_power_list2 <- vector("list", length(noise_param))


noise <- noise_param

# We will collect all results for this noise in a numeric vector;
# length = length(dim_param) * length(k_vec)
power_hdim <- numeric(length(dim_param) * length(k_vec))

# We will fill 'power_hdim' in chunks corresponding to each k
start_idx  <- 1

for (k in k_vec) {
  
  # Use sapply to compute power for each dim_p
  power_vec <- sapply(dim_param, function(dim_p) {
    
    # replicate(...) helps avoid manual loops for iter times
    pval_vec <- replicate(iter, {
      df       <- dat.gen.high.dim.fixed(k, dim_p, noise)
      curv_fit <- intrinsic.curv.est(df, L = 8)
      as.numeric(curv_fit$pval < 0.05)
    })
    
    mean(pval_vec)  # average of 0/1 to get the power
  })
  
  # Insert the results into the correct location of power_hdim
  end_idx <- start_idx + length(dim_param) - 1
  power_hdim[start_idx:end_idx] <- power_vec
  
  start_idx <- end_idx + 1
}

hdim_power_list2[[1]] <- power_hdim

plot(dim_param, hdim_power_list2[[1]][1:11],xlim = c(4, 100), ylim = c(0, 1.05), lty = 1, type = 'l',col=col[1],cex.lab=2,cex.axis=2,cex.main=2,lwd=3,xlab="Dimension p" , ylab = "Power", main=expression(paste("Fixed Signal-to-Noise-Ratio ", sigma, " = 3/(10",sqrt(p),")")))
lines(dim_param, hdim_power_list2[[1]][12:22], col=col[2], lty = 2, lwd = 3)
lines(dim_param, hdim_power_list2[[1]][23:33], col=col[3], lty = 3, lwd = 3)
legend("bottomright", legend = c("n = 200", "n = 500", "n = 1000"),lty=c(1,2,3),lwd=3, col = col,cex=2)
abline(h = 0.05, lty = 2)
abline(h = 1, lty = 1)
