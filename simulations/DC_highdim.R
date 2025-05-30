# -------------------------------------------------------------------------
# Power and Type-I error analysis under High-dimensional ambient Euclidean random objects with low intrinsic dimension
#
# Manuscript reference : Section S.1.4
# Figure reproduced    : Figure S.4
# Simulation code reference : simulations/DC_highdim.R
# -------------------------------------------------------------------------


## Import main functions
source("src/DC_mainfunctions.R")
source("simulations/DC_simulations_datagen.R")


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
    print(power_vec)
    
    start_idx <- end_idx + 1
  }
  
  hdim_power_list[[n_idx]] <- power_hdim
}


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
  print(power_vec)
  
  start_idx <- end_idx + 1
}

hdim_power_list2[[1]] <- power_hdim

