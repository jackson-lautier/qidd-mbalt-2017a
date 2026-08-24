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
# ".\raw-data\2024_FET_Analysis_Public_US_20260112T094254.txt"
#
# "./code/mbalt-term-time-function.R"
# "./code/fet-death-cols.csv
# "./code/2024-date-ref.csv"
#
#The code must be run sequentially downwards.
#As the new, cleaned files are prepared, they will be saved in a new
#folder 'processed-data' in the wd.
#For data analysis, proceed directly to 'data_analysis.R'.
#
#
################################################################################
################################################################################
################################################################################
################################################################################
# FETAL DEATHS
################################################################################
################################################################################
################################################################################
################################################################################

#where processed data will be stored
dir.create('./processed-data/')

require('readr')
require('lubridate')

col.specs = read.csv("./code/fet-death-cols.csv")

col_specs = fwf_positions(
  start = col.specs$start,
  end = col.specs$end,
  col_names = col.specs$col.name
)

data <- read_fwf("./raw-data/2024_FET_Analysis_Public_US_20260112T094254.txt",
                 col_specs)

#write.csv(data, "./data/2024-fet-analysis-public.csv")
#dat <- read.csv("./data/2024-fet-analysis-public.csv")

dat = as.data.frame(data)

#remove unknown gestation ages
dat = dat[ !(dat$COMBGEST == 99), ]

date.ref = read.csv("./code/2024-date-ref.csv")

DOD = c()

for(i in c(1:nrow(dat))){
  
  #step one: pull possible dates
  poss.dates =
    date.ref[ (as.numeric(dat$DOD_MM[i]) == date.ref$dt.mth) &
                (as.numeric(dat$DOD_WK[i]) == date.ref$dt.wk) &
                (as.numeric(dat$DOD_YY[i]) == date.ref$dt.year), ]
  
  #step two: select date at random from possible dates
  set.seed(i)
  rnd.date = sample(poss.dates$full.date, 1)
  
  #step three: append random date
  DOD = append(rnd.date, DOD)
  
}

dat$est.DOD = DOD #estimate date of death

#impute weeks from report period start
report.date = "01-01-2024"
dat$est.wks.rep.date = interval(as.Date(report.date,"%m-%d-%Y"),
                                as.Date(dat$est.DOD, "%m/%d/%Y")) %/% weeks(1)

#get sample of LT lifetimes; if est.wks.rep.date < est.gest.age
dat = dat[ (as.numeric(dat$COMBGEST) > dat$est.wks.rep.date ), ]

#define LT random variable
dat$Yi = as.numeric(dat$COMBGEST) - dat$est.wks.rep.date
dat$Xi = as.numeric(dat$COMBGEST)

write.csv(dat, "./processed-data/2024-fet-analysis-LT-lifetimes.csv")

################################################################################
################################################################################
################################################################################
################################################################################
# CONSUMER AUTO LEASE
################################################################################
################################################################################
################################################################################
################################################################################

source("./code/mbalt-term-time-function.R")

################################################################################
################################################################################
################################################################################
################################################################################

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
