source("~/src/DC_mainfunctions.R")


# Read in the Excel file with multiple sheets
my_data <- lapply(excel_sheets("~/energy_data/generation_monthly.xlsx"), read_excel, path = "~/energy_data/generation_monthly.xlsx")
# View the data in each sheet
names(my_data) <- excel_sheets("~/energy_data/generation_monthly.xlsx")

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
    
    other_list[month] = total_list[[month]] - cope_list[[month]] - ng_list[[month]] #- nc_list[[month]]
    
  }
  
  df_tmp = data.frame(
    "CoalPet" = as.numeric(unlist(cope_list)),
    "NaturalGas" = as.numeric(unlist(ng_list)),
    "Other" = as.numeric(unlist(other_list))
  )
  
  df_energy = rbind(df_energy, df_tmp)
  
}

comp.data = as.matrix(sqrt(df_energy[,1:3]/rowSums(df_energy[,1:3])))

# k: number of observations
k = nrow(comp.data)

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

res_energy = intrinsic.curv.est(distmat = dist.mat.energy, L = 8)

isomap_result <- isomap(dist.mat.energy, k = 8, ndim = 1)
isomap_rep = as.vector(isomap_result$points)

time.stamp = seq(0,1,0.1)

geod.sph = data.frame()
for (tt in time.stamp){
  
  temp.geo = find_weight(min(isomap_rep),max(isomap_rep),tt,isomap_rep, h =0.01)
  geod.sph = rbind(geod.sph, find_weighted_barycenter_on_sphere(comp.data, temp.geo))
  
}
colnames(geod.sph) = c("CoalPet","NaturalGas","Other")

colors <- colorRampPalette(c("cadetblue", "coral3"))(nrow(comp.data))
# Create the 3D scatter plot using scatterplot3d and store the returned object
s3d <- scatterplot3d(comp.data[, 1], comp.data[, 2], comp.data[, 3], color = colors, pch = 16,
                     main = "Proportion of energy source",
                     xlab = "sqrt( Coal + Pet )", ylab = "sqrt( Natural Gas )", zlab = "sqrt( Others )",
                     cex.main = 3, cex.axis = 2, cex.lab = 2, box = FALSE)

# Add lines to the scatter plot
# For example, adding lines connecting the points in the order they appear in comp.data
for (i in 1:(nrow(geod.sph) - 1)) {
  s3d$points3d(x = c(geod.sph[i, 1], geod.sph[i + 1, 1]),
               y = c(geod.sph[i, 2], geod.sph[i + 1, 2]),
               z = c(geod.sph[i, 3], geod.sph[i + 1, 3]),
               type = "l", col = "grey30", lty = 1,lwd = 3)  # Line color set to black
}


# Set up layout: main plot on the left (75% width) and color bar on the right (25% width)
# Create the main plot
years <- seq(2012, 2021, length.out = 10)

layout(matrix(1:2, ncol = 2), widths = c(3,1))

# Create the main plot
par(mar = c(5, 5, 4, 1))  # Adjust margins for the main plot
plot(1:120, isomap_rep, col = colors, pch = 16,
     xlab = "Year", ylab = "ISOMAP Representation",
     cex.main = 2, cex.axis = 2, cex.lab = 2, main = "", xaxt = "n")

# Add custom x-axis with year labels
axis(1, at = seq(1, 120, length.out = 4), labels = as.character(years)[c(1,4,7,10)], cex.axis = 2)

# Set up margins for the color bar and plot an empty plot area
par(mar =  c(5, 1, 5, 10))  # Increase the right margin to make space for a wider legend
plot(1:10, (1:10)*10, type="n", bty="n", xaxt = "n", yaxt = "n", ann = FALSE)

# Add a color bar legend using image.plot()
years <- seq(2012, 2021, length.out = 10)
image.plot(legend.only = TRUE, 
           zlim = c(1, 120),  # Scale according to the data size
           col = colorRampPalette(c("cadetblue", "coral3"))(100), 
           axis.args = list(at = seq(1, 120, length.out = 10), labels = as.character(years), cex.axis = 2),
           legend.mar = 6,  # Increased space on the right for the color bar
           legend.width = 8,  # Adjust the width of the legend bar
           legend.args = list(text = "Year", side = 3, line = 1, cex = 2))


## Plot confidence region for the joint distribution of intrinsic metric and Frechet variance
cov.est.energy = res_energy$cov.normalized
Vm.est.energy = res_energy$Vm
Vf.est.energy = res_energy$Vf

t=2
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((t+1))
par(mfrow=c(1,1),mar = c(5,5,3,2))
plot(ellipse(cov.est.energy,center = c(Vm.est.energy,Vf.est.energy),level=0.99),type='l',col=col[1],lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="")
points(Vm.est.energy,Vf.est.energy,pch=16)
lines(ellipse(cov.est.energy,center = c(Vm.est.energy,Vf.est.energy),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.energy,center = c(Vm.est.energy,Vf.est.energy),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

