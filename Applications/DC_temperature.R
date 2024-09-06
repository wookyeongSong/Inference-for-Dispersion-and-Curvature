source("~/src/DC_mainfunctions.R")

##########################################
#### Real Data Analysis: Airport Data ####
##########################################

## Get a list of all CSV files in a directory
file_list <- list.files(path = "~/airport data", pattern = "\\.csv$", full.names = TRUE)

## Initialize an empty list to store data frames
data_list <- list()

## Loop through the list of files and read each CSV file
for (file in file_list) {
  data <- read.csv(file)
  data_list[[file]] <- data
}

n = length(data_list)

## Obtain name of airports
names.airports = c()
for(i in 1:n){
  
  first_file_data <- data_list[[i]]
  names.airports = c(names.airports, (first_file_data[1,2]))
  
}

## Fix NAs
data_list[[9]][122,c(8,9)] = c(62,47)
data_list[[20]][1719,9] = 70

######################################################################
#### Summer period (Jun 21 - Sep 20, 63 years, considering leap days)

## Data Preprocessing
l = 63
summer = c()
for (i in 1:l){
  summer = c(summer, (172+365*(i-1)):(263+365*(i-1)))
}

tmp =60
leapdays = c(tmp)
for(i in 1:15){
  
  tmp = tmp + 365*4 + 1
  leapdays = c(leapdays, tmp)
  
}

df.temp = list()
for(i in 1:n){
  df.temp[[i]] = data_list[[i]][-leapdays,]
  df.temp[[i]] = df.temp[[i]][summer,]
}

## Convert time-series data to exceedance distribution
qSup = seq(0,1,length.out = 201)
Vm.summer = c()
asym.Vm.summer = c()
for(sta in 1:n){

  ## starting from summer, seattle airport
  df.temp.max = df.temp[[sta]][,"TMAX"]
  
  qf.temp.max.station = sapply(1:l, function(i){
    
    yin = df.temp.max[(1+92*(i-1)):(92*i)]
    df = CreateDensity(yin)
    qf = dens2quantile(dens = df$y, dSup = df$x,qSup=qSup)
    qf
    
  })
  
  temperature.curv.fit = functional.curv.est(qSup, qf.temp.max.station)
  Vm.summer = c(Vm.summer,temperature.curv.fit$Vm)
  asym.Vm.summer = c(asym.Vm.summer, sqrt(temperature.curv.fit$cov.normalized[1,1]))
}

low.vari.summer = order(Vm.summer)[4:1]
high.vari.summer = order(Vm.summer)[39:36]

names.summer.whole = names.airports[order(Vm.summer)]
Vm.summer.whole = Vm.summer[order(Vm.summer)]
asym.Vm.summer.whole = asym.Vm.summer[order(Vm.summer)]

names.abb.summer = c("BOI", "RNO", "ATL","SLC", "MSY", "JAX","MIA","PBI")
data.summer = data.frame(value = c(Vm.summer[high.vari.summer],Vm.summer[low.vari.summer]), se = c(asym.Vm.summer[high.vari.summer],asym.Vm.summer[low.vari.summer]))

col = rainbow_hcl(3)

par(mar=c(5,5,3,2))
plotCI(x= 1:8, y= data.summer$value, uiw = 1.96*data.summer$se, liw = 1.96*data.summer$se,ylab = expression(V[M]),xlab = "Airport",axes=FALSE,scol="black",col=c(rep(col[1],4),rep(col[3],4)),cex=2,pch=19,font.axis=2, cex.axis=2, cex.lab=2, cex.main = 2.5,main="Variance of summer temperature")
axis(side=2,cex.axis=2)         ## add default y-axis (ticks+labels)
axis(side=1,at=1:8,  ## add custom x-axis
     label=names.abb.summer,font.axis=1, cex.axis=2, cex.lab=2)
box(bty="l")    

######################################################################
#### Winter period (Dec 22 - Mar 21, 63 years, considering leap days)

## Data Preprocessing
l = 63
winter= c()
for (i in 1:l){
  winter = c(winter, (355+365*(i-1)):(444+365*(i-1)))
}

tmp =60
leapdays = c(tmp)
for(i in 1:15){
  
  tmp = tmp + 365*4 + 1
  leapdays = c(leapdays, tmp)
  
}

df.temp = list()
for(i in 1:n){
  df.temp[[i]] = data_list[[i]][-leapdays,]
  df.temp[[i]] = df.temp[[i]][winter,]
}

qSup = seq(0,1,length.out = 201)
Vm.winter = c()
asym.Vm.winter = c()
for(sta in 1:n){
  
  ## starting from winter, seattle airport
  df.temp.min = df.temp[[sta]][,"TMIN"]
  
  qf.temp.min.station = sapply(1:l, function(i){
    
    yin = df.temp.min[(1+90*(i-1)):(90*i)]
    df = CreateDensity(yin)
    qf = dens2quantile(dens = df$y, dSup = df$x,qSup=qSup)
    qf
    
  })
  
  temperature.curv.fit = functional.curv.est(qSup, qf.temp.min.station)
  Vm.winter = c(Vm.winter,temperature.curv.fit$Vm)
  asym.Vm.winter = c(asym.Vm.winter, sqrt(temperature.curv.fit$cov.normalized[1,1]))
}

low.vari.winter = order(Vm.winter)[4:1]
high.vari.winter = order(Vm.winter)[39:36]

names.abb.winter = c("ANC", "MSP", "ORD","IND", "MIA", "SEA","SAN","LAX")
data.winter = data.frame(value = c(Vm.winter[high.vari.winter],Vm.winter[low.vari.winter]),se = c(asym.Vm.winter[high.vari.winter],asym.Vm.winter[low.vari.winter]))

par(mar=c(5,5,3,2))
plotCI(x=1:8, y= data.winter$value, uiw = 1.96*data.winter$se, liw = 1.96*data.winter$se,ylab = expression(V[M]),xlab = "Airport",axes=FALSE,scol="black",col=c(rep(col[1],4),rep(col[3],4)),cex=2,pch=19,font.axis=2, cex.axis=2, cex.lab=2, cex.main = 2.5,main="Variance of winter temperature")
axis(side=2,cex.axis=2)         ## add default y-axis (ticks+labels)
axis(side=1,at=1:8,  ## add custom x-axis
     label=names.abb.winter,font.axis=1, cex.axis=2, cex.lab=2)
box(bty="l")   
