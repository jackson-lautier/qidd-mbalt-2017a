######################################################################################
######################################################################################
######################################################################################
# Data analysis scripts to produce the data used in the manuscript:
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

#dependent files

# supporting files:
# ".\raw-data\MBALT17A_compiled2.csv"
#
# "./code/lrt-formulas-LT-only.R"
# "./code/lrt-formulas-LT-RC.R"

#where any supporting analysis will be stored
dir.create('./results/')



################################################################################
################################################################################
# APPLICATION SECTION
################################################################################
################################################################################
rm(list=ls())

################################################################################
#QIDD tests
################################################################################

#terms for testing:
lease.terms = c(24, 27, 30, 36, 39, 42, 48)

results = matrix(NA, nrow = length(lease.terms), ncol = 12)
colnames(results) = c("lease.term",
                      "n",
                      "Delta",
                      "Xi",
                      "m",
                      "censor.rate",
                      "An.0",
                      "uLRT",
                      "Lambda.n",
                      "deg.free",
                      "crit.value",
                      "result")
rownames(results) = lease.terms

################################################################################
#left-truncation
################################################################################

l = lease.terms[1]

f.name = paste('./data-clean/mbalt-2017-', l, 'mo.csv', sep="")
obs_data = read.csv(f.name)
obs_data = obs_data[,-c(1,4)]
colnames(obs_data) = c("Xi", "Yi")

f.name = paste('./data-clean/mbalt-2017-', l, 'mo-trapezoid-dim.csv', sep="")
trap.param = read.csv(f.name)
trap.param = trap.param[,-1]

Delta = trap.param$delta
omega = trap.param$xi
tau = trap.param$tau
m = min(trap.param$m, max(obs_data$Yi) - Delta)
epsilon = trap.param$e
xi = min(omega, epsilon - 1)

minU = Delta + 1
maxU = min(omega, epsilon - 1)
minV = Delta + 1
maxV = Delta + m

U.support = c(minU:maxU)
V.support = c(minV:maxV)

u = c( (Delta + 1) : omega)
reps = length(c( (Delta + 1) : (Delta + m)))

x_col = c()
y_col = c()
for(u in c( (Delta + 1) : omega)){
  
  for(v in c( (Delta + 1) : (Delta + m))){
    if(v <= u){
      x_col = append(x_col, u)
      y_col = append(y_col, v)
    }
  }
  
}

h_inv = data.frame("X" = x_col,
                   "Y" = y_col)

source('./code/lrt-formulas-LT-only.R')

#null: h_* is true
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
rev_haz_est = sapply( c((Delta + 1):(Delta +m)), bnx)
U_est.n = sapply(c( (Delta + 1) : xi), f_est)
G_est.n = sapply(c((Delta + 1):(Delta + m)), g_est)
l.0.n = log_like_fn_0(U_est.n, G_est.n)

#f.est.n = U_est.n
#g.est.n = G_est.n
alpha.est.n = alpha(U_est.n, G_est.n)

OMEGA_est = c()
for(j in c((Delta + 1):(omega))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    #print(c(j,k))
    OMEGA_est = append(OMEGA_est,
                       sum(( (obs_data$Yi == k ) & (obs_data$Xi == j) )) / nrow(obs_data))
  }
}
l.1.n = log_like_fn(OMEGA_est)

n = nrow(obs_data)
n.0 = sum(OMEGA_est == 0)

#create missing sample
missing.dat = as.data.frame(h_inv[which(OMEGA_est == 0),c(1:2)])
names(missing.dat) = c("Xi", "Yi")
obs_data = rbind(obs_data, missing.dat)

#check l.1
OMEGA_est = c()
for(j in c((Delta + 1):(omega))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    #print(c(j,k))
    OMEGA_est = append(OMEGA_est,
                       sum(( (obs_data$Yi == k ) & (obs_data$Xi == j) )) / nrow(obs_data))
  }
}
l.1.n.star = log_like_fn(OMEGA_est)

l.1.n + n * log( n / (n + n.0) ) + n.0 * log(1 / (n.0 + n))
l.1.n.star

#check l.0
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
rev_haz_est = sapply( c((Delta + 1):(Delta +m)), bnx)
U_est.n.star = sapply(c( (Delta + 1) : xi), f_est)
G_est.n.star = sapply(c((Delta + 1):(Delta + m)), g_est)
l.0.n.star = log_like_fn_0(U_est.n.star, G_est.n.star)

#f.est.n.star = U_est
#g.est.n.star = G_est
alpha.est.n.star = alpha(U_est.n.star, G_est.n.star)

#calculate K2
s1 = -(n + n.0) * log( alpha.est.n.star / alpha.est.n ) -
  (n.0) * log( alpha.est.n )

new.dat = obs_data[c( (n+1) : (n + n.0) ),]
old.dat = obs_data[c(1:n),]

s2 = c()
for(u in c((Delta + 1):(omega))){
  for(v in c((Delta + 1):min(u, Delta + m))){
    
    a1 = sum( (new.dat$Yi == v) & (new.dat$Xi == u) )
    a2 = log( f_X(u, U_est.n.star) )
    a3 = sum( (old.dat$Yi == v) & (old.dat$Xi == u) )
    a4 = log( f_X(u, U_est.n.star) / f_X(u, U_est.n) )
    a5 = ifelse(a3 == 0, 0, a3 * a4)
    s2 = append(s2, a1 * a2 + a5)
  }
}

s3 = c()
for(u in c((Delta + 1):(omega))){
  for(v in c((Delta + 1):min(u, Delta + m))){
    
    a1 = sum( (new.dat$Yi == v) & (new.dat$Xi == u) )
    a2 = log( g_Y(v, G_est.n.star) )
    a3 = sum( (old.dat$Yi == v) & (old.dat$Xi == u) )
    a4 = log( g_Y(v, G_est.n.star) / g_Y(v, G_est.n) )
    a5 = ifelse(a3 == 0, 0, a3 * a4)
    s3 = append(s3, a1 * a2 + a5)
  }
}

#check
l.0.n + s1 + sum(s2) + sum(s3)
l.0.n.star

K1 = n * log( n / (n + n.0) ) + n.0 * log(1 / (n.0 + n))
K2 = s1 + sum(s2) + sum(s3)

deg.free = (omega - Delta - 1/2) * m - (m^2) / 2 - omega + Delta + 1

lambda.n = -2 * (l.0.n - l.1.n) - 2 * (K2 - K1)

#store results
results[as.character(l), "lease.term"] = l
results[as.character(l), "n"] = n
results[as.character(l), "Delta"] = Delta
results[as.character(l), "Xi"] = xi
results[as.character(l), "m"] = m
results[as.character(l), "censor.rate"] = 0
results[as.character(l), "An.0"] = n.0
results[as.character(l), "uLRT"] = -2 * (l.0.n - l.1.n)
results[as.character(l), "Lambda.n"] = lambda.n
results[as.character(l), "deg.free"] = deg.free
results[as.character(l), "crit.value"] = qchisq(0.05, deg.free, lower.tail = FALSE)
results[as.character(l), "result"] = 1 * (lambda.n > qchisq(0.05, deg.free, lower.tail = FALSE))

################################################################################
l = lease.terms[2]

f.name = paste('./data-clean/mbalt-2017-', l, 'mo.csv', sep="")
obs_data = read.csv(f.name)
obs_data = obs_data[,-c(1,4)]
colnames(obs_data) = c("Xi", "Yi")

f.name = paste('./data-clean/mbalt-2017-', l, 'mo-trapezoid-dim.csv', sep="")
trap.param = read.csv(f.name)
trap.param = trap.param[,-1]

Delta = trap.param$delta
omega = trap.param$xi
tau = trap.param$tau
m = min(trap.param$m, max(obs_data$Yi) - Delta)
epsilon = trap.param$e
xi = min(omega, epsilon - 1)

minU = Delta + 1
maxU = min(omega, epsilon - 1)
minV = Delta + 1
maxV = Delta + m

U.support = c(minU:maxU)
V.support = c(minV:maxV)

u = c( (Delta + 1) : omega)
reps = length(c( (Delta + 1) : (Delta + m)))

x_col = c()
y_col = c()
for(u in c( (Delta + 1) : omega)){
  
  for(v in c( (Delta + 1) : (Delta + m))){
    if(v <= u){
      x_col = append(x_col, u)
      y_col = append(y_col, v)
    }
  }
  
}

h_inv = data.frame("X" = x_col,
                   "Y" = y_col)

source('./code/lrt-formulas-LT-only.R')

#null: h_* is true
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
rev_haz_est = sapply( c((Delta + 1):(Delta +m)), bnx)
U_est.n = sapply(c( (Delta + 1) : xi), f_est)
G_est.n = sapply(c((Delta + 1):(Delta + m)), g_est)
l.0.n = log_like_fn_0(U_est.n, G_est.n)

#f.est.n = U_est.n
#g.est.n = G_est.n
alpha.est.n = alpha(U_est.n, G_est.n)

OMEGA_est = c()
for(j in c((Delta + 1):(omega))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    #print(c(j,k))
    OMEGA_est = append(OMEGA_est,
                       sum(( (obs_data$Yi == k ) & (obs_data$Xi == j) )) / nrow(obs_data))
  }
}
l.1.n = log_like_fn(OMEGA_est)

n = nrow(obs_data)
n.0 = sum(OMEGA_est == 0)

#create missing sample
missing.dat = as.data.frame(h_inv[which(OMEGA_est == 0),c(1:2)])
names(missing.dat) = c("Xi", "Yi")
obs_data = rbind(obs_data, missing.dat)

#check l.1
OMEGA_est = c()
for(j in c((Delta + 1):(omega))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    #print(c(j,k))
    OMEGA_est = append(OMEGA_est,
                       sum(( (obs_data$Yi == k ) & (obs_data$Xi == j) )) / nrow(obs_data))
  }
}
l.1.n.star = log_like_fn(OMEGA_est)

l.1.n + n * log( n / (n + n.0) ) + n.0 * log(1 / (n.0 + n))
l.1.n.star

#check l.0
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
rev_haz_est = sapply( c((Delta + 1):(Delta +m)), bnx)
U_est.n.star = sapply(c( (Delta + 1) : xi), f_est)
G_est.n.star = sapply(c((Delta + 1):(Delta + m)), g_est)
l.0.n.star = log_like_fn_0(U_est.n.star, G_est.n.star)

#f.est.n.star = U_est
#g.est.n.star = G_est
alpha.est.n.star = alpha(U_est.n.star, G_est.n.star)

#calculate K2
s1 = -(n + n.0) * log( alpha.est.n.star / alpha.est.n ) -
  (n.0) * log( alpha.est.n )

new.dat = obs_data[c( (n+1) : (n + n.0) ),]
old.dat = obs_data[c(1:n),]

s2 = c()
for(u in c((Delta + 1):(omega))){
  for(v in c((Delta + 1):min(u, Delta + m))){
    
    a1 = sum( (new.dat$Yi == v) & (new.dat$Xi == u) )
    a2 = log( f_X(u, U_est.n.star) )
    a3 = sum( (old.dat$Yi == v) & (old.dat$Xi == u) )
    a4 = log( f_X(u, U_est.n.star) / f_X(u, U_est.n) )
    a5 = ifelse(a3 == 0, 0, a3 * a4)
    s2 = append(s2, a1 * a2 + a5)
  }
}

s3 = c()
for(u in c((Delta + 1):(omega))){
  for(v in c((Delta + 1):min(u, Delta + m))){
    
    a1 = sum( (new.dat$Yi == v) & (new.dat$Xi == u) )
    a2 = log( g_Y(v, G_est.n.star) )
    a3 = sum( (old.dat$Yi == v) & (old.dat$Xi == u) )
    a4 = log( g_Y(v, G_est.n.star) / g_Y(v, G_est.n) )
    a5 = ifelse(a3 == 0, 0, a3 * a4)
    s3 = append(s3, a1 * a2 + a5)
  }
}

#check
l.0.n + s1 + sum(s2) + sum(s3)
l.0.n.star

K1 = n * log( n / (n + n.0) ) + n.0 * log(1 / (n.0 + n))
K2 = s1 + sum(s2) + sum(s3)

deg.free = (omega - Delta - 1/2) * m - (m^2) / 2 - omega + Delta + 1

lambda.n = -2 * (l.0.n - l.1.n) - 2 * (K2 - K1)

#store results
results[as.character(l), "lease.term"] = l
results[as.character(l), "n"] = n
results[as.character(l), "Delta"] = Delta
results[as.character(l), "Xi"] = xi
results[as.character(l), "m"] = m
results[as.character(l), "censor.rate"] = 0
results[as.character(l), "An.0"] = n.0
results[as.character(l), "uLRT"] = -2 * (l.0.n - l.1.n)
results[as.character(l), "Lambda.n"] = lambda.n
results[as.character(l), "deg.free"] = deg.free
results[as.character(l), "crit.value"] = qchisq(0.05, deg.free, lower.tail = FALSE)
results[as.character(l), "result"] = 1 * (lambda.n > qchisq(0.05, deg.free, lower.tail = FALSE))

################################################################################

l = lease.terms[3]

f.name = paste('./data-clean/mbalt-2017-', l, 'mo.csv', sep="")
obs_data = read.csv(f.name)
obs_data = obs_data[,-c(1,4)]
colnames(obs_data) = c("Xi", "Yi")

f.name = paste('./data-clean/mbalt-2017-', l, 'mo-trapezoid-dim.csv', sep="")
trap.param = read.csv(f.name)
trap.param = trap.param[,-1]

Delta = trap.param$delta
omega = trap.param$xi
tau = trap.param$tau
m = min(trap.param$m, max(obs_data$Yi) - Delta)
epsilon = trap.param$e
xi = min(omega, epsilon - 1)

minU = Delta + 1
maxU = min(omega, epsilon - 1)
minV = Delta + 1
maxV = Delta + m

U.support = c(minU:maxU)
V.support = c(minV:maxV)

u = c( (Delta + 1) : omega)
reps = length(c( (Delta + 1) : (Delta + m)))

x_col = c()
y_col = c()
for(u in c( (Delta + 1) : omega)){
  
  for(v in c( (Delta + 1) : (Delta + m))){
    if(v <= u){
      x_col = append(x_col, u)
      y_col = append(y_col, v)
    }
  }
  
}

h_inv = data.frame("X" = x_col,
                   "Y" = y_col)

source('./code/lrt-formulas-LT-only.R')

#null: h_* is true
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
rev_haz_est = sapply( c((Delta + 1):(Delta +m)), bnx)
U_est.n = sapply(c( (Delta + 1) : xi), f_est)
G_est.n = sapply(c((Delta + 1):(Delta + m)), g_est)
l.0.n = log_like_fn_0(U_est.n, G_est.n)

#f.est.n = U_est.n
#g.est.n = G_est.n
alpha.est.n = alpha(U_est.n, G_est.n)

OMEGA_est = c()
for(j in c((Delta + 1):(omega))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    #print(c(j,k))
    OMEGA_est = append(OMEGA_est,
                       sum(( (obs_data$Yi == k ) & (obs_data$Xi == j) )) / nrow(obs_data))
  }
}
l.1.n = log_like_fn(OMEGA_est)

n = nrow(obs_data)
n.0 = sum(OMEGA_est == 0)

#create missing sample
missing.dat = as.data.frame(h_inv[which(OMEGA_est == 0),c(1:2)])
names(missing.dat) = c("Xi", "Yi")
obs_data = rbind(obs_data, missing.dat)

#check l.1
OMEGA_est = c()
for(j in c((Delta + 1):(omega))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    #print(c(j,k))
    OMEGA_est = append(OMEGA_est,
                       sum(( (obs_data$Yi == k ) & (obs_data$Xi == j) )) / nrow(obs_data))
  }
}
l.1.n.star = log_like_fn(OMEGA_est)

l.1.n + n * log( n / (n + n.0) ) + n.0 * log(1 / (n.0 + n))
l.1.n.star

#check l.0
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
rev_haz_est = sapply( c((Delta + 1):(Delta +m)), bnx)
U_est.n.star = sapply(c( (Delta + 1) : xi), f_est)
G_est.n.star = sapply(c((Delta + 1):(Delta + m)), g_est)
l.0.n.star = log_like_fn_0(U_est.n.star, G_est.n.star)

#f.est.n.star = U_est
#g.est.n.star = G_est
alpha.est.n.star = alpha(U_est.n.star, G_est.n.star)

#calculate K2
s1 = -(n + n.0) * log( alpha.est.n.star / alpha.est.n ) -
  (n.0) * log( alpha.est.n )

new.dat = obs_data[c( (n+1) : (n + n.0) ),]
old.dat = obs_data[c(1:n),]

s2 = c()
for(u in c((Delta + 1):(omega))){
  for(v in c((Delta + 1):min(u, Delta + m))){
    
    a1 = sum( (new.dat$Yi == v) & (new.dat$Xi == u) )
    a2 = log( f_X(u, U_est.n.star) )
    a3 = sum( (old.dat$Yi == v) & (old.dat$Xi == u) )
    a4 = log( f_X(u, U_est.n.star) / f_X(u, U_est.n) )
    a5 = ifelse(a3 == 0, 0, a3 * a4)
    s2 = append(s2, a1 * a2 + a5)
  }
}

s3 = c()
for(u in c((Delta + 1):(omega))){
  for(v in c((Delta + 1):min(u, Delta + m))){
    
    a1 = sum( (new.dat$Yi == v) & (new.dat$Xi == u) )
    a2 = log( g_Y(v, G_est.n.star) )
    a3 = sum( (old.dat$Yi == v) & (old.dat$Xi == u) )
    a4 = log( g_Y(v, G_est.n.star) / g_Y(v, G_est.n) )
    a5 = ifelse(a3 == 0, 0, a3 * a4)
    s3 = append(s3, a1 * a2 + a5)
  }
}

#check
l.0.n + s1 + sum(s2) + sum(s3)
l.0.n.star

K1 = n * log( n / (n + n.0) ) + n.0 * log(1 / (n.0 + n))
K2 = s1 + sum(s2) + sum(s3)

deg.free = (omega - Delta - 1/2) * m - (m^2) / 2 - omega + Delta + 1

lambda.n = -2 * (l.0.n - l.1.n) - 2 * (K2 - K1)

#store results
results[as.character(l), "lease.term"] = l
results[as.character(l), "n"] = n
results[as.character(l), "Delta"] = Delta
results[as.character(l), "Xi"] = xi
results[as.character(l), "m"] = m
results[as.character(l), "censor.rate"] = 0
results[as.character(l), "An.0"] = n.0
results[as.character(l), "uLRT"] = -2 * (l.0.n - l.1.n)
results[as.character(l), "Lambda.n"] = lambda.n
results[as.character(l), "deg.free"] = deg.free
results[as.character(l), "crit.value"] = qchisq(0.05, deg.free, lower.tail = FALSE)
results[as.character(l), "result"] = 1 * (lambda.n > qchisq(0.05, deg.free, lower.tail = FALSE))

################################################################################
#right-censoring
################################################################################

l = lease.terms[4]

f.name = paste('./data-clean/mbalt-2017-', l, 'mo.csv', sep="")
obs_data = read.csv(f.name)
obs_data = obs_data[,-c(1)]
colnames(obs_data) = c("Zi", "Yi", "Di")

f.name = paste('./data-clean/mbalt-2017-', l, 'mo-trapezoid-dim.csv', sep="")
trap.param = read.csv(f.name)
trap.param = trap.param[,-1]

Delta = trap.param$delta
omega = trap.param$xi
tau = trap.param$tau
m = min(trap.param$m, max(obs_data$Yi) - Delta)
epsilon = trap.param$e
xi = min(omega, epsilon - 1)

minU = Delta + 1
maxU = min(omega, epsilon - 1)
minV = Delta + 1
maxV = Delta + m

source('./code/lrt-formulas-LT-RC.R')

#create free parameter space reference
R.v = c()
R.u = c()
for(j in c((Delta + 1):(xi))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    if(j <= (k + tau)){
      R.v = append(R.v, k)
      R.u = append(R.u, j)
    }
  }
}


R.p.v = c()
R.p.u = c()
for(j in c((Delta + 1):(xi))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    if(j == (k + tau + 1)){
      R.p.v = append(R.p.v, k)
      R.p.u = append(R.p.u, j)
    }
  }
}

OMEGA.1 = rep(1/(length(R.v) + length(R.p.v)), length(R.v))

df.OMEGA.1 = data.frame("u" = R.u,
                        "v" = R.v,
                        "O1" = OMEGA.1)

OMEGA.2 = rep(1/(length(R.v) + length(R.p.v)), length(R.p.v))

df.OMEGA.2 = data.frame("u" = R.p.u,
                        "v" = R.p.v,
                        "O2" = OMEGA.2)

OMEGA = c(OMEGA.1, OMEGA.2)

U.support = c(minU:maxU)
V.support = c(minV:maxV)

u = c( (Delta + 1) : omega)
reps = length(c( (Delta + 1) : (Delta + m)))

x_col = c()
y_col = c()
for(u in c( (Delta + 1) : omega)){
  
  for(v in c( (Delta + 1) : (Delta + m))){
    if(v <= u){
      x_col = append(x_col, u)
      y_col = append(y_col, v)
    }
  }
  
}

h_inv = data.frame("X" = x_col,
                   "Y" = y_col)

n = nrow(obs_data)

############################################################################
#unadjusted likelihood calculations
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
U_est.0 = sapply(c( (Delta + 1) : xi), f_est)
G_est.0 = sapply(c((Delta + 1):(Delta + m)), g_est)
l_0 = log_like_fn_0(c(U_est.0, G_est.0)) #l0(n)

#alt: h_* not true, free to be all param
OMEGA_est.0 =  c(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v),
                 mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v))

l_A = log_like_fn(OMEGA_est.0) #l1(n)

############################################################################
#create the missing sample, append to obs_data
card.R0 = length(which(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v) == 0))
card.R.prime0 = length(which(mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v) == 0))

if( sum(card.R0, card.R.prime0) > 0){
  
  if(card.R0 > 0){
    R0 = df.OMEGA.1[which(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v) == 0),c("u","v")]
    R0.C = df.OMEGA.1[which(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v) > 0),c("u","v")]
    colnames(R0) = c("Zi", "Yi")
    colnames(R0.C) = c("Zi", "Yi")
    R0$Di = 1
    R0.C$Di = 1
    R0 = R0[,c("Yi", "Zi", "Di")]
    R0.C = R0.C[,c("Yi", "Zi", "Di")]
  }
  
  if(card.R.prime0 > 0){
    R.prime0 = df.OMEGA.2[which(mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v) == 0),c("u","v")]
    R.prime0.C = df.OMEGA.2[which(mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v) > 0),c("u","v")]
    colnames(R.prime0) = c("Zi", "Yi")
    colnames(R.prime0.C) = c("Zi", "Yi")
    R.prime0$Di = 0
    R.prime0.C$Di = 0
    R.prime0 = R.prime0[,c("Yi", "Zi", "Di")]
    R.prime0.C = R.prime0.C[,c("Yi", "Zi", "Di")]
  }
  
  if( (card.R0 > 0) & (card.R.prime0 > 0) ){
    R.prime0.add = R.prime0
    R.prime0.add$Zi = R.prime0.add$Zi - 1
    adj.sample = rbind(obs_data, R0, R.prime0.add)
  }
  
  if( (card.R0 > 0) & (card.R.prime0 == 0) ){
    adj.sample = rbind(obs_data, R0)
  }
  
  if( (card.R0 == 0) & (card.R.prime0 > 0) ){
    R.prime0.add = R.prime0
    R.prime0.add$Zi = R.prime0.add$Zi - 1
    adj.sample = rbind(obs_data, R.prime0.add)
  }
  
  obs_data = adj.sample
  
}

############################################################################
#calculate the adjusted likelihoods plus K1, K2
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
U_est = sapply(c( (Delta + 1) : xi), f_est)
G_est = sapply(c((Delta + 1):(Delta + m)), g_est)
l_0.adj = log_like_fn_0(c(U_est, G_est)) #l0(n.star)

#calculate K1
n0 = card.R0
n0.prime = card.R.prime0

theta.n = c(U_est.0, G_est.0)
theta.n.star = c(U_est, G_est)

s1 = -(n + n0 + n0.prime) * log( alpha(theta.n.star) / alpha(theta.n) ) -
  (n0 + n0.prime) * log( alpha(theta.n))

new.dat = obs_data[c( (n+1) : (n + n0 + n0.prime) ),]
old.dat = obs_data[c(1:n),]

s2 = c()
for(v in c(minV:maxV)){
  
  a1 = sum( new.dat$Yi == v ) * log( g_Y(v, theta.n.star) )
  a2 = sum( old.dat$Yi == v) * log( g_Y(v, theta.n.star) / g_Y(v, theta.n) )
  s2 = append(s2, a1 + a2)
  
}

s3 = c()
for(u in c(minU:maxU)){
  
  a1 = sum( (new.dat$Zi == u) & (new.dat$Di == 1) ) * log( f_X(u, theta.n.star) )
  a2 = sum( (old.dat$Zi == u) & (old.dat$Di == 1) ) * 
    log( f_X(u, theta.n.star) / f_X(u, theta.n) )
  a1 = ifelse(is.na(a1), 0, a1)
  a2 = ifelse(is.na(a2), 0, a2)
  s3 = append(s3, a1 + a2)
  
}

s4 = c()
for(u in c(minU:(maxU-1))){
  
  a1 = sum( (new.dat$Zi == u) & (new.dat$Di == 0) )
  a2 = sum(sapply(c((u + 1):(maxU)), f_X, THETA = theta.n.star ))
  a3 = sum( (old.dat$Zi == u) & (old.dat$Di == 0) )
  a4 = sum(sapply(c((u + 1):(maxU)), f_X, THETA = theta.n ))
  
  a2.s = ifelse(a2 == 0, 0, log(a2))
  a4.s = ifelse( a4 == 0, 0, ifelse( (a2/a4) == 0, 0, log(a2 / a4)))
  
  s4 = append(s4, a1 * a2.s + a3 * a4.s )
  
}

K2 = s1 + sum(s2) + sum(s3) + sum(s4)

#check
l_0 + K2
l_0.adj


#alt: h_* not true, free to be all param
OMEGA_est =  c(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v),
               mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v))

l_A.adj = log_like_fn(OMEGA_est)

K1 = n * log( n / (n + n0 + n0.prime) ) +
  (n0 + n0.prime) * log( 1 / ( n + n0 + n0.prime ))

#check
l_A.adj
l_A + K1

lambda.tau.n = -2 * (l_0 - l_A) - 2 * (K2 - K1)

a = (Delta + m) - (xi - (tau + 1))
deg.free = m * (tau + 2) -
  (a * (a + 1)) / 2 -
  (xi + m - (Delta + 1))

results[as.character(l), "lease.term"] = l
results[as.character(l), "n"] = n
results[as.character(l), "Delta"] = Delta
results[as.character(l), "Xi"] = xi
results[as.character(l), "m"] = m
results[as.character(l), "censor.rate"] = sum(obs_data$Di == 0)/n
results[as.character(l), "An.0"] = n0 + n0.prime
results[as.character(l), "uLRT"] = -2 * (l_0 - l_A)
results[as.character(l), "Lambda.n"] = -2 * (l_0 - l_A) - 2 * (K2 - K1)
results[as.character(l), "deg.free"] = deg.free
results[as.character(l), "crit.value"] = qchisq(0.05, deg.free, lower.tail = FALSE)
results[as.character(l), "result"] = 1 * (lambda.tau.n > qchisq(0.05, deg.free, lower.tail = FALSE))

################################################################################
l = lease.terms[5]

f.name = paste('./data-clean/mbalt-2017-', l, 'mo.csv', sep="")
obs_data = read.csv(f.name)
obs_data = obs_data[,-c(1)]
colnames(obs_data) = c("Zi", "Yi", "Di")

f.name = paste('./data-clean/mbalt-2017-', l, 'mo-trapezoid-dim.csv', sep="")
trap.param = read.csv(f.name)
trap.param = trap.param[,-1]

Delta = trap.param$delta
omega = trap.param$xi
tau = trap.param$tau
m = min(trap.param$m, max(obs_data$Yi) - Delta)
epsilon = trap.param$e
xi = min(omega, epsilon - 1)

minU = Delta + 1
maxU = min(omega, epsilon - 1)
minV = Delta + 1
maxV = Delta + m

#check conditions
obs_data = obs_data[-which( (obs_data$Zi > (obs_data$Yi + tau)) & (obs_data$Di == 1) ),]

source('./code/lrt-formulas-LT-RC.R')

#create free parameter space reference
R.v = c()
R.u = c()
for(j in c((Delta + 1):(xi))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    if(j <= (k + tau)){
      R.v = append(R.v, k)
      R.u = append(R.u, j)
    }
  }
}


R.p.v = c()
R.p.u = c()
for(j in c((Delta + 1):(xi))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    if(j == (k + tau + 1)){
      R.p.v = append(R.p.v, k)
      R.p.u = append(R.p.u, j)
    }
  }
}

OMEGA.1 = rep(1/(length(R.v) + length(R.p.v)), length(R.v))

df.OMEGA.1 = data.frame("u" = R.u,
                        "v" = R.v,
                        "O1" = OMEGA.1)

OMEGA.2 = rep(1/(length(R.v) + length(R.p.v)), length(R.p.v))

df.OMEGA.2 = data.frame("u" = R.p.u,
                        "v" = R.p.v,
                        "O2" = OMEGA.2)

OMEGA = c(OMEGA.1, OMEGA.2)

U.support = c(minU:maxU)
V.support = c(minV:maxV)

u = c( (Delta + 1) : omega)
reps = length(c( (Delta + 1) : (Delta + m)))

x_col = c()
y_col = c()
for(u in c( (Delta + 1) : omega)){
  
  for(v in c( (Delta + 1) : (Delta + m))){
    if(v <= u){
      x_col = append(x_col, u)
      y_col = append(y_col, v)
    }
  }
  
}

h_inv = data.frame("X" = x_col,
                   "Y" = y_col)

n = nrow(obs_data)

############################################################################
#unadjusted likelihood calculations
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
U_est.0 = sapply(c( (Delta + 1) : xi), f_est)
G_est.0 = sapply(c((Delta + 1):(Delta + m)), g_est)
l_0 = log_like_fn_0(c(U_est.0, G_est.0)) #l0(n)

#alt: h_* not true, free to be all param
OMEGA_est.0 =  c(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v),
                 mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v))

l_A = log_like_fn(OMEGA_est.0) #l1(n)

############################################################################
#create the missing sample, append to obs_data
card.R0 = length(which(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v) == 0))
card.R.prime0 = length(which(mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v) == 0))

if( sum(card.R0, card.R.prime0) > 0){
  
  if(card.R0 > 0){
    R0 = df.OMEGA.1[which(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v) == 0),c("u","v")]
    R0.C = df.OMEGA.1[which(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v) > 0),c("u","v")]
    colnames(R0) = c("Zi", "Yi")
    colnames(R0.C) = c("Zi", "Yi")
    R0$Di = 1
    R0.C$Di = 1
    R0 = R0[,c("Yi", "Zi", "Di")]
    R0.C = R0.C[,c("Yi", "Zi", "Di")]
  }
  
  if(card.R.prime0 > 0){
    R.prime0 = df.OMEGA.2[which(mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v) == 0),c("u","v")]
    R.prime0.C = df.OMEGA.2[which(mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v) > 0),c("u","v")]
    colnames(R.prime0) = c("Zi", "Yi")
    colnames(R.prime0.C) = c("Zi", "Yi")
    R.prime0$Di = 0
    R.prime0.C$Di = 0
    R.prime0 = R.prime0[,c("Yi", "Zi", "Di")]
    R.prime0.C = R.prime0.C[,c("Yi", "Zi", "Di")]
  }
  
  if( (card.R0 > 0) & (card.R.prime0 > 0) ){
    R.prime0.add = R.prime0
    R.prime0.add$Zi = R.prime0.add$Zi - 1
    adj.sample = rbind(obs_data, R0, R.prime0.add)
  }
  
  if( (card.R0 > 0) & (card.R.prime0 == 0) ){
    adj.sample = rbind(obs_data, R0)
  }
  
  if( (card.R0 == 0) & (card.R.prime0 > 0) ){
    R.prime0.add = R.prime0
    R.prime0.add$Zi = R.prime0.add$Zi - 1
    adj.sample = rbind(obs_data, R.prime0.add)
  }
  
  obs_data = adj.sample
  
}

############################################################################
#calculate the adjusted likelihoods plus K1, K2
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
U_est = sapply(c( (Delta + 1) : xi), f_est)
G_est = sapply(c((Delta + 1):(Delta + m)), g_est)
l_0.adj = log_like_fn_0(c(U_est, G_est)) #l0(n.star)

#calculate K1
n0 = card.R0
n0.prime = card.R.prime0

theta.n = c(U_est.0, G_est.0)
theta.n.star = c(U_est, G_est)

s1 = -(n + n0 + n0.prime) * log( alpha(theta.n.star) / alpha(theta.n) ) -
  (n0 + n0.prime) * log( alpha(theta.n))

new.dat = obs_data[c( (n+1) : (n + n0 + n0.prime) ),]
old.dat = obs_data[c(1:n),]

s2 = c()
for(v in c(minV:maxV)){
  
  a1 = sum( new.dat$Yi == v ) * log( g_Y(v, theta.n.star) )
  a2 = sum( old.dat$Yi == v) * log( g_Y(v, theta.n.star) / g_Y(v, theta.n) )
  s2 = append(s2, a1 + a2)
  
}

s3 = c()
for(u in c(minU:maxU)){
  
  a1 = sum( (new.dat$Zi == u) & (new.dat$Di == 1) ) * log( f_X(u, theta.n.star) )
  a2 = sum( (old.dat$Zi == u) & (old.dat$Di == 1) ) * 
    log( f_X(u, theta.n.star) / f_X(u, theta.n) )
  a1 = ifelse(is.na(a1), 0, a1)
  a2 = ifelse(is.na(a2), 0, a2)
  s3 = append(s3, a1 + a2)
  
}

s4 = c()
for(u in c(minU:(maxU-1))){
  
  a1 = sum( (new.dat$Zi == u) & (new.dat$Di == 0) )
  a2 = sum(sapply(c((u + 1):(maxU)), f_X, THETA = theta.n.star ))
  a3 = sum( (old.dat$Zi == u) & (old.dat$Di == 0) )
  a4 = sum(sapply(c((u + 1):(maxU)), f_X, THETA = theta.n ))
  
  a2.s = ifelse(a2 == 0, 0, log(a2))
  a4.s = ifelse( a4 == 0, 0, ifelse( (a2/a4) == 0, 0, log(a2 / a4)))
  
  s4 = append(s4, a1 * a2.s + a3 * a4.s )
  
}

K2 = s1 + sum(s2) + sum(s3) + sum(s4)

#check
l_0 + K2
l_0.adj


#alt: h_* not true, free to be all param
OMEGA_est =  c(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v),
               mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v))

l_A.adj = log_like_fn(OMEGA_est)

K1 = n * log( n / (n + n0 + n0.prime) ) +
  (n0 + n0.prime) * log( 1 / ( n + n0 + n0.prime ))

#check
l_A.adj
l_A + K1

lambda.tau.n = -2 * (l_0 - l_A) - 2 * (K2 - K1)

a = (Delta + m) - (xi - (tau + 1))
deg.free = m * (tau + 2) -
  (a * (a + 1)) / 2 -
  (xi + m - (Delta + 1))

results[as.character(l), "lease.term"] = l
results[as.character(l), "n"] = n
results[as.character(l), "Delta"] = Delta
results[as.character(l), "Xi"] = xi
results[as.character(l), "m"] = m
results[as.character(l), "censor.rate"] = sum(obs_data$Di == 0)/n
results[as.character(l), "An.0"] = n0 + n0.prime
results[as.character(l), "uLRT"] = -2 * (l_0 - l_A)
results[as.character(l), "Lambda.n"] = -2 * (l_0 - l_A) - 2 * (K2 - K1)
results[as.character(l), "deg.free"] = deg.free
results[as.character(l), "crit.value"] = qchisq(0.05, deg.free, lower.tail = FALSE)
results[as.character(l), "result"] = 1 * (lambda.tau.n > qchisq(0.05, deg.free, lower.tail = FALSE))

################################################################################
l = lease.terms[6]

f.name = paste('./data-clean/mbalt-2017-', l, 'mo.csv', sep="")
obs_data = read.csv(f.name)
obs_data = obs_data[,-c(1)]
colnames(obs_data) = c("Zi", "Yi", "Di")

f.name = paste('./data-clean/mbalt-2017-', l, 'mo-trapezoid-dim.csv', sep="")
trap.param = read.csv(f.name)
trap.param = trap.param[,-1]

Delta = trap.param$delta
omega = trap.param$xi
tau = trap.param$tau
m = min(trap.param$m, max(obs_data$Yi) - Delta)
epsilon = trap.param$e
xi = min(omega, epsilon - 1)

minU = Delta + 1
maxU = min(omega, epsilon - 1)
minV = Delta + 1
maxV = Delta + m

#check conditions
#obs_data = obs_data[-which((obs_data$Zi != (obs_data$Yi + tau)) & (obs_data$Di == 0) ),]

source('./code/lrt-formulas-LT-RC.R')

#create free parameter space reference
R.v = c()
R.u = c()
for(j in c((Delta + 1):(xi))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    if(j <= (k + tau)){
      R.v = append(R.v, k)
      R.u = append(R.u, j)
    }
  }
}


R.p.v = c()
R.p.u = c()
for(j in c((Delta + 1):(xi))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    if(j == (k + tau + 1)){
      R.p.v = append(R.p.v, k)
      R.p.u = append(R.p.u, j)
    }
  }
}

OMEGA.1 = rep(1/(length(R.v) + length(R.p.v)), length(R.v))

df.OMEGA.1 = data.frame("u" = R.u,
                        "v" = R.v,
                        "O1" = OMEGA.1)

OMEGA.2 = rep(1/(length(R.v) + length(R.p.v)), length(R.p.v))

df.OMEGA.2 = data.frame("u" = R.p.u,
                        "v" = R.p.v,
                        "O2" = OMEGA.2)

OMEGA = c(OMEGA.1, OMEGA.2)

U.support = c(minU:maxU)
V.support = c(minV:maxV)

u = c( (Delta + 1) : omega)
reps = length(c( (Delta + 1) : (Delta + m)))

x_col = c()
y_col = c()
for(u in c( (Delta + 1) : omega)){
  
  for(v in c( (Delta + 1) : (Delta + m))){
    if(v <= u){
      x_col = append(x_col, u)
      y_col = append(y_col, v)
    }
  }
  
}

h_inv = data.frame("X" = x_col,
                   "Y" = y_col)

n = nrow(obs_data)

############################################################################
#unadjusted likelihood calculations
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
U_est.0 = sapply(c( (Delta + 1) : xi), f_est)
G_est.0 = sapply(c((Delta + 1):(Delta + m)), g_est)
l_0 = log_like_fn_0(c(U_est.0, G_est.0)) #l0(n)

#alt: h_* not true, free to be all param
OMEGA_est.0 =  c(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v),
                 mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v))

l_A = log_like_fn(OMEGA_est.0) #l1(n)

############################################################################
#create the missing sample, append to obs_data
card.R0 = length(which(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v) == 0))
card.R.prime0 = length(which(mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v) == 0))

if( sum(card.R0, card.R.prime0) > 0){
  
  if(card.R0 > 0){
    R0 = df.OMEGA.1[which(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v) == 0),c("u","v")]
    R0.C = df.OMEGA.1[which(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v) > 0),c("u","v")]
    colnames(R0) = c("Zi", "Yi")
    colnames(R0.C) = c("Zi", "Yi")
    R0$Di = 1
    R0.C$Di = 1
    R0 = R0[,c("Yi", "Zi", "Di")]
    R0.C = R0.C[,c("Yi", "Zi", "Di")]
  }
  
  if(card.R.prime0 > 0){
    R.prime0 = df.OMEGA.2[which(mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v) == 0),c("u","v")]
    R.prime0.C = df.OMEGA.2[which(mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v) > 0),c("u","v")]
    colnames(R.prime0) = c("Zi", "Yi")
    colnames(R.prime0.C) = c("Zi", "Yi")
    R.prime0$Di = 0
    R.prime0.C$Di = 0
    R.prime0 = R.prime0[,c("Yi", "Zi", "Di")]
    R.prime0.C = R.prime0.C[,c("Yi", "Zi", "Di")]
  }
  
  if( (card.R0 > 0) & (card.R.prime0 > 0) ){
    R.prime0.add = R.prime0
    R.prime0.add$Zi = R.prime0.add$Zi - 1
    adj.sample = rbind(obs_data, R0, R.prime0.add)
  }
  
  if( (card.R0 > 0) & (card.R.prime0 == 0) ){
    adj.sample = rbind(obs_data, R0)
  }
  
  if( (card.R0 == 0) & (card.R.prime0 > 0) ){
    R.prime0.add = R.prime0
    R.prime0.add$Zi = R.prime0.add$Zi - 1
    adj.sample = rbind(obs_data, R.prime0.add)
  }
  
  obs_data = adj.sample
  
}

############################################################################
#calculate the adjusted likelihoods plus K1, K2
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
U_est = sapply(c( (Delta + 1) : xi), f_est)
G_est = sapply(c((Delta + 1):(Delta + m)), g_est)
l_0.adj = log_like_fn_0(c(U_est, G_est)) #l0(n.star)

#calculate K1
n0 = card.R0
n0.prime = card.R.prime0

theta.n = c(U_est.0, G_est.0)
theta.n.star = c(U_est, G_est)

s1 = -(n + n0 + n0.prime) * log( alpha(theta.n.star) / alpha(theta.n) ) -
  (n0 + n0.prime) * log( alpha(theta.n))

new.dat = obs_data[c( (n+1) : (n + n0 + n0.prime) ),]
old.dat = obs_data[c(1:n),]

s2 = c()
for(v in c(minV:maxV)){
  
  a1 = sum( new.dat$Yi == v ) * log( g_Y(v, theta.n.star) )
  a2 = sum( old.dat$Yi == v) * log( g_Y(v, theta.n.star) / g_Y(v, theta.n) )
  s2 = append(s2, a1 + a2)
  
}

s3 = c()
for(u in c(minU:maxU)){
  
  a1 = sum( (new.dat$Zi == u) & (new.dat$Di == 1) ) * log( f_X(u, theta.n.star) )
  a2 = sum( (old.dat$Zi == u) & (old.dat$Di == 1) ) * 
    log( f_X(u, theta.n.star) / f_X(u, theta.n) )
  a1 = ifelse(is.na(a1), 0, a1)
  a2 = ifelse(is.na(a2), 0, a2)
  s3 = append(s3, a1 + a2)
  
}

s4 = c()
for(u in c(minU:(maxU-1))){
  
  a1 = sum( (new.dat$Zi == u) & (new.dat$Di == 0) )
  a2 = sum(sapply(c((u + 1):(maxU)), f_X, THETA = theta.n.star ))
  a3 = sum( (old.dat$Zi == u) & (old.dat$Di == 0) )
  a4 = sum(sapply(c((u + 1):(maxU)), f_X, THETA = theta.n ))
  
  a2.s = ifelse(a2 == 0, 0, log(a2))
  a4.s = ifelse( a4 == 0, 0, ifelse( (a2/a4) == 0, 0, log(a2 / a4)))
  
  s4 = append(s4, a1 * a2.s + a3 * a4.s )
  
}

K2 = s1 + sum(s2) + sum(s3) + sum(s4)

#check
l_0 + K2
l_0.adj


#alt: h_* not true, free to be all param
OMEGA_est =  c(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v),
               mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v))

l_A.adj = log_like_fn(OMEGA_est)

K1 = n * log( n / (n + n0 + n0.prime) ) +
  (n0 + n0.prime) * log( 1 / ( n + n0 + n0.prime ))

#check
l_A.adj
l_A + K1

lambda.tau.n = -2 * (l_0 - l_A) - 2 * (K2 - K1)

a = (Delta + m) - (xi - (tau + 1))
deg.free = m * (tau + 2) -
  (a * (a + 1)) / 2 -
  (xi + m - (Delta + 1))

results[as.character(l), "lease.term"] = l
results[as.character(l), "n"] = n
results[as.character(l), "Delta"] = Delta
results[as.character(l), "Xi"] = xi
results[as.character(l), "m"] = m
results[as.character(l), "censor.rate"] = sum(obs_data$Di == 0)/n
results[as.character(l), "An.0"] = n0 + n0.prime
results[as.character(l), "uLRT"] = -2 * (l_0 - l_A)
results[as.character(l), "Lambda.n"] = -2 * (l_0 - l_A) - 2 * (K2 - K1)
results[as.character(l), "deg.free"] = deg.free
results[as.character(l), "crit.value"] = qchisq(0.05, deg.free, lower.tail = FALSE)
results[as.character(l), "result"] = 1 * (lambda.tau.n > qchisq(0.05, deg.free, lower.tail = FALSE))

################################################################################
l = lease.terms[7]

f.name = paste('./data-clean/mbalt-2017-', l, 'mo.csv', sep="")
obs_data = read.csv(f.name)
obs_data = obs_data[,-c(1)]
colnames(obs_data) = c("Zi", "Yi", "Di")

f.name = paste('./data-clean/mbalt-2017-', l, 'mo-trapezoid-dim.csv', sep="")
trap.param = read.csv(f.name)
trap.param = trap.param[,-1]

Delta = trap.param$delta
omega = trap.param$xi
tau = trap.param$tau
m = min(trap.param$m, max(obs_data$Yi) - Delta)
epsilon = trap.param$e
xi = min(omega, epsilon - 1)

minU = Delta + 1
maxU = min(omega, epsilon - 1)
minV = Delta + 1
maxV = Delta + m

#check conditions
obs_data = obs_data[-which( (obs_data$Zi > (obs_data$Yi + tau)) & (obs_data$Di == 1) ),]

source('./code/lrt-formulas-LT-RC.R')

#create free parameter space reference
R.v = c()
R.u = c()
for(j in c((Delta + 1):(xi))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    if(j <= (k + tau)){
      R.v = append(R.v, k)
      R.u = append(R.u, j)
    }
  }
}


R.p.v = c()
R.p.u = c()
for(j in c((Delta + 1):(xi))){
  for(k in c((Delta + 1):min(j, Delta + m))){
    if(j == (k + tau + 1)){
      R.p.v = append(R.p.v, k)
      R.p.u = append(R.p.u, j)
    }
  }
}

OMEGA.1 = rep(1/(length(R.v) + length(R.p.v)), length(R.v))

df.OMEGA.1 = data.frame("u" = R.u,
                        "v" = R.v,
                        "O1" = OMEGA.1)

OMEGA.2 = rep(1/(length(R.v) + length(R.p.v)), length(R.p.v))

df.OMEGA.2 = data.frame("u" = R.p.u,
                        "v" = R.p.v,
                        "O2" = OMEGA.2)

OMEGA = c(OMEGA.1, OMEGA.2)

U.support = c(minU:maxU)
V.support = c(minV:maxV)

u = c( (Delta + 1) : omega)
reps = length(c( (Delta + 1) : (Delta + m)))

x_col = c()
y_col = c()
for(u in c( (Delta + 1) : omega)){
  
  for(v in c( (Delta + 1) : (Delta + m))){
    if(v <= u){
      x_col = append(x_col, u)
      y_col = append(y_col, v)
    }
  }
  
}

h_inv = data.frame("X" = x_col,
                   "Y" = y_col)

n = nrow(obs_data)

############################################################################
#unadjusted likelihood calculations
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
U_est.0 = sapply(c( (Delta + 1) : xi), f_est)
G_est.0 = sapply(c((Delta + 1):(Delta + m)), g_est)
l_0 = log_like_fn_0(c(U_est.0, G_est.0)) #l0(n)

#alt: h_* not true, free to be all param
OMEGA_est.0 =  c(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v),
                 mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v))

l_A = log_like_fn(OMEGA_est.0) #l1(n)

############################################################################
#create the missing sample, append to obs_data
card.R0 = length(which(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v) == 0))
card.R.prime0 = length(which(mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v) == 0))

if( sum(card.R0, card.R.prime0) > 0){
  
  if(card.R0 > 0){
    R0 = df.OMEGA.1[which(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v) == 0),c("u","v")]
    R0.C = df.OMEGA.1[which(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v) > 0),c("u","v")]
    colnames(R0) = c("Zi", "Yi")
    colnames(R0.C) = c("Zi", "Yi")
    R0$Di = 1
    R0.C$Di = 1
    R0 = R0[,c("Yi", "Zi", "Di")]
    R0.C = R0.C[,c("Yi", "Zi", "Di")]
  }
  
  if(card.R.prime0 > 0){
    R.prime0 = df.OMEGA.2[which(mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v) == 0),c("u","v")]
    R.prime0.C = df.OMEGA.2[which(mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v) > 0),c("u","v")]
    colnames(R.prime0) = c("Zi", "Yi")
    colnames(R.prime0.C) = c("Zi", "Yi")
    R.prime0$Di = 0
    R.prime0.C$Di = 0
    R.prime0 = R.prime0[,c("Yi", "Zi", "Di")]
    R.prime0.C = R.prime0.C[,c("Yi", "Zi", "Di")]
  }
  
  if( (card.R0 > 0) & (card.R.prime0 > 0) ){
    R.prime0.add = R.prime0
    R.prime0.add$Zi = R.prime0.add$Zi - 1
    adj.sample = rbind(obs_data, R0, R.prime0.add)
  }
  
  if( (card.R0 > 0) & (card.R.prime0 == 0) ){
    adj.sample = rbind(obs_data, R0)
  }
  
  if( (card.R0 == 0) & (card.R.prime0 > 0) ){
    R.prime0.add = R.prime0
    R.prime0.add$Zi = R.prime0.add$Zi - 1
    adj.sample = rbind(obs_data, R.prime0.add)
  }
  
  obs_data = adj.sample
  
}

############################################################################
#calculate the adjusted likelihoods plus K1, K2
haz_est = sapply(c( (Delta + 1) : xi), lnx)
if ( is.na(haz_est[length(U.support)]) ){ haz_est[length(U.support)] = 1 }
U_est = sapply(c( (Delta + 1) : xi), f_est)
G_est = sapply(c((Delta + 1):(Delta + m)), g_est)
l_0.adj = log_like_fn_0(c(U_est, G_est)) #l0(n.star)

#calculate K1
n0 = card.R0
n0.prime = card.R.prime0

theta.n = c(U_est.0, G_est.0)
theta.n.star = c(U_est, G_est)

s1 = -(n + n0 + n0.prime) * log( alpha(theta.n.star) / alpha(theta.n) ) -
  (n0 + n0.prime) * log( alpha(theta.n))

new.dat = obs_data[c( (n+1) : (n + n0 + n0.prime) ),]
old.dat = obs_data[c(1:n),]

s2 = c()
for(v in c(minV:maxV)){
  
  a1 = sum( new.dat$Yi == v ) * log( g_Y(v, theta.n.star) )
  a2 = sum( old.dat$Yi == v) * log( g_Y(v, theta.n.star) / g_Y(v, theta.n) )
  s2 = append(s2, a1 + a2)
  
}

s3 = c()
for(u in c(minU:maxU)){
  
  a1 = sum( (new.dat$Zi == u) & (new.dat$Di == 1) ) * log( f_X(u, theta.n.star) )
  a2 = sum( (old.dat$Zi == u) & (old.dat$Di == 1) ) * 
    log( f_X(u, theta.n.star) / f_X(u, theta.n) )
  a1 = ifelse(is.na(a1), 0, a1)
  a2 = ifelse(is.na(a2), 0, a2)
  s3 = append(s3, a1 + a2)
  
}

s4 = c()
for(u in c(minU:(maxU-1))){
  
  a1 = sum( (new.dat$Zi == u) & (new.dat$Di == 0) )
  a2 = sum(sapply(c((u + 1):(maxU)), f_X, THETA = theta.n.star ))
  a3 = sum( (old.dat$Zi == u) & (old.dat$Di == 0) )
  a4 = sum(sapply(c((u + 1):(maxU)), f_X, THETA = theta.n ))
  
  a2.s = ifelse(a2 == 0, 0, log(a2))
  a4.s = ifelse( a4 == 0, 0, ifelse( (a2/a4) == 0, 0, log(a2 / a4)))
  
  s4 = append(s4, a1 * a2.s + a3 * a4.s )
  
}

K2 = s1 + sum(s2) + sum(s3) + sum(s4)

#check
l_0 + K2
l_0.adj


#alt: h_* not true, free to be all param
OMEGA_est =  c(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v),
               mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v))

l_A.adj = log_like_fn(OMEGA_est)

K1 = n * log( n / (n + n0 + n0.prime) ) +
  (n0 + n0.prime) * log( 1 / ( n + n0 + n0.prime ))

#check
l_A.adj
l_A + K1

lambda.tau.n = -2 * (l_0 - l_A) - 2 * (K2 - K1)

a = (Delta + m) - (xi - (tau + 1))
deg.free = m * (tau + 2) -
  (a * (a + 1)) / 2 -
  (xi + m - (Delta + 1))

results[as.character(l), "lease.term"] = l
results[as.character(l), "n"] = n
results[as.character(l), "Delta"] = Delta
results[as.character(l), "Xi"] = xi
results[as.character(l), "m"] = m
results[as.character(l), "censor.rate"] = sum(obs_data$Di == 0)/n
results[as.character(l), "An.0"] = n0 + n0.prime
results[as.character(l), "uLRT"] = -2 * (l_0 - l_A)
results[as.character(l), "Lambda.n"] = -2 * (l_0 - l_A) - 2 * (K2 - K1)
results[as.character(l), "deg.free"] = deg.free
results[as.character(l), "crit.value"] = qchisq(0.05, deg.free, lower.tail = FALSE)
results[as.character(l), "result"] = 1 * (lambda.tau.n > qchisq(0.05, deg.free, lower.tail = FALSE))

write.csv(results, "./results/table-1.csv")
