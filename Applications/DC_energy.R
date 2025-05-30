# -------------------------------------------------------------------------
# Data Application: Energy source data analysis
#
# Data available at Applications/energy_data
# Manuscript reference : Section 6.2
# Figure reproduced    : Figure 7
# Data application code reference : Applications/DC_energy.R
# -------------------------------------------------------------------------

## Import main functions
source("src/DC_mainfunctions.R")


############################
## Data Preprocessing 

# Read in the Excel file with multiple sheets
my_data <- lapply(excel_sheets("Applications/energy_data/generation_monthly.xlsx"), read_excel, path = "Applications/energy_data/generation_monthly.xlsx")
# View the data in each sheet
names(my_data) <- excel_sheets("Applications/energy_data/generation_monthly.xlsx")

df_energy = data.frame()
for (k in 6:15){
  
  colnm = my_data[k][[1]][4,]
  df = my_data[k][[1]][-c(1:4),]
  colnames(df) = colnm
  df = df[df[,"STATE"] == "US-Total" & df[,"TYPE OF PRODUCER"] == "Total Electric Power Industry",]
  
  cope_list = list()
  for (month in 1:12){
    
    tmp = df[df["MONTH"]==month ,]
    cope_list[month] = sum(as.numeric(unlist(tmp[tmp["ENERGY SOURCE"]=="Coal" | tmp["ENERGY SOURCE"]=="Petroleum",6])))
    
  }
  
  ng_list = list()
  for (month in 1:12){
    
    tmp = df[df["MONTH"]==month ,]
    ng_list[month] = sum(as.numeric(unlist(tmp[tmp["ENERGY SOURCE"]=="Natural Gas",6])))
    
  }
  
  total_list = list()
  for (month in 1:12){
    
    tmp = df[df["MONTH"]==month ,]
    total_list[month] = sum(as.numeric(unlist(tmp[tmp["ENERGY SOURCE"]=="Total",6])))
    
  }
  
  other_list = list()
  for (month in 1:12){
    
    other_list[month] = total_list[[month]] - cope_list[[month]] - ng_list[[month]]
    
  }
  
  df_tmp = data.frame(
    "CoalPet" = as.numeric(unlist(cope_list)),
    "NaturalGas" = as.numeric(unlist(ng_list)),
    "Other" = as.numeric(unlist(other_list))
  )
  
  df_energy = rbind(df_energy, df_tmp)
  
}


## Compositional Energy Source Data
comp.data = as.matrix(sqrt(df_energy[,1:3]/rowSums(df_energy[,1:3])))


## k: number of observations
k = nrow(comp.data)


## Calculate pairwise spherical distance 
dist.mat.energy = matrix(0,nrow=k,ncol=k)
for(row in 1:k){
  for(col in 1:k){
    
    if (row == col){
      
      dist.mat.energy[row,col] = 0
      
    } else{
      
      dist.mat.energy[row,col] = acos(comp.data[row,] %*% comp.data[col,])
      
    }
    
  }
}


## Intrinsic curvature test results
res_energy = intrinsic.curv.est(distmat = dist.mat.energy, L = 8)
cov.est.energy = res_energy$cov.normalized
Vm.est.energy = res_energy$Vm
Vf.est.energy = res_energy$Vf

## ISOMAP Representation
isomap_result <- isomap(dist.mat.energy, k = 8, ndim = 1)
isomap_rep = as.vector(isomap_result$points)




