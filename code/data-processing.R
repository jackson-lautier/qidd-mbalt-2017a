######################################################################################
######################################################################################
######################################################################################
# Data processing scripts to produce the data used in the manuscript:
# 
# "Testing quasi-independence for discrete data subject to left-truncation"
#
# LAUTIER, CHIOU
# 2025
#
# Computer and R version details
# _                           
# platform       x86_64-w64-mingw32               
# arch           x86_64                           
# os             mingw32                          
# crt            ucrt                             
# system         x86_64, mingw32                  
# status                                          
# major          4                                
# minor          5.1                              
# year           2025                             
# month          06                               
# day            13                               
# svn rev        88306                            
# language       R                                
# version.string R version 4.5.1 (2025-06-13 ucrt)
# nickname       Great Square Root         

######################################################################################
######################################################################################
######################################################################################
######################################################################################
# INSTRUCTIONS
#
# supporting files:
# ".\raw-data\MBALT17A_compiled2.csv"
#
# "./code/mbalt-term-time-function.R"
#
#The code must be run sequentially downwards.
#As the new, cleaned files are prepared, they will be saved in a new
#folder 'processed-data' in the wd.
#For data analysis, proceed directly to 'data_analysis.R'.
#
#
######################################################################################
######################################################################################
######################################################################################
######################################################################################
require('lubridate')

source("./code/mbalt-term-time-function.R")

######################################################################################
######################################################################################
######################################################################################
######################################################################################

#where processed data will be stored
dir.create('./processed-data/')

#mbalt 2017 summary statistics
path = "./raw-data/"
mbalt <- read.csv(paste(path,'MBALT17A_compiled2.csv',sep=""))

table(mbalt$originalLeaseTermNumber)

date <- paste(mbalt$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
min(date); max(date)

min(as.numeric(mbalt$lesseeCreditScore), na.rm = TRUE)
mean(as.numeric(mbalt$lesseeCreditScore), na.rm = TRUE)
median(as.numeric(mbalt$lesseeCreditScore), na.rm = TRUE)
max(as.numeric(mbalt$lesseeCreditScore), na.rm = TRUE)

min(mbalt$originalLeaseTermNumber); max(mbalt$originalLeaseTermNumber)

#get table of original loan terms
table(mbalt$originalLeaseTermNumber)

#terms for testing:
lease.terms =
  as.numeric(
    names(
      table(mbalt$originalLeaseTermNumber))[which(table(mbalt$originalLeaseTermNumber) > 400)])

for(l in lease.terms){
  
  lease_term = l
  obs_window = 28
  
  mlease_term <- mbalt[mbalt$originalLeaseTermNumber == lease_term,]
  
  delta = lease_term - max(mlease_term$remainingTermNumber) 
  M = lease_term - min(mlease_term$remainingTermNumber) - delta
  
  e = (M + delta + 1) + obs_window - 1
  
  T_start = M + delta + mlease_term$remainingTermNumber - lease_term 
  Y = M + delta - T_start + 1
  
  mlease_term <- mlease_term[,1:(ncol(mbalt)-8*(28-(e-M-delta)))]
  
  X = vector()
  C = vector()
  for (j in c(1:nrow(mlease_term))) {
    out = term_time(mlease_term[j,])
    X = append(X,out[1])
    C = append(C,out[2])
  }
  
  Xc = X
  
  table(Xc)
  Xc = ifelse(Xc >= lease_term + 1, lease_term + 1, Xc)
  table(Xc)
  
  obs_data <- data.frame("Zi" = Xc,
                         "Yi" = Y,
                         "Di" = C)
  
  #censoring check for omega based on above adjustment
  omega = max(obs_data$Zi)
  
  obs_data[(obs_data$Zi == omega) & (obs_data$Di == 0),]
  obs_data$Di = ifelse((obs_data$Zi == omega) & (obs_data$Di == 0),
                       1,
                       obs_data$Di)
  
  Delta = max(delta, min(obs_data$Yi) - 1)
  M = max(obs_data$Yi) - Delta 
  epsilon = obs_window + (M + Delta)
  tau = epsilon - (M + Delta + 1)
  xi = min(omega, epsilon - 1)
  
  #censoring data check #2
  cens.dat = obs_data[obs_data$Di == 0,]
  cens.dat$check = cens.dat$Zi - cens.dat$Yi
  cens.dat = cens.dat[cens.dat$check == tau,]
  
  no.cens.dat = obs_data[obs_data$Di == 1,]
  obs_data = rbind(no.cens.dat, cens.dat[,c(1:3)])
  
  f.name = paste('./processed-data/mbalt-2017-',
                 l,
                 'mo.csv', sep = "")
  write.csv(obs_data, f.name)
  
  mbalt.2017.l.parameters = data.frame("delta" = Delta,
                                       "m" = M,
                                       "xi" = xi,
                                       "e" = epsilon,
                                       "tau" = tau)
  
  f.name = paste('./processed-data/mbalt-2017-',
                 l,
                 'mo-trapezoid-dim.csv', sep = "")
  write.csv(mbalt.2017.l.parameters,
            f.name)
  
  print(paste(which(l == lease.terms), " of ",
              length(lease.terms), " complete!"))
  
}
