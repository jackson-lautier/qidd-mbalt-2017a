dn = function(u, v){
  
  n = nrow(obs_data)
  
  dn = sum(( (obs_data$Yi == v) &
               (obs_data$Zi == u) &
               (obs_data$Di == 1))) / n
  
  return(dn)
  
}

cn = function(u, v){
  
  n = nrow(obs_data)
  
  cn = sum(( (obs_data$Yi == v) &
               (obs_data$Zi == (u-1)) &
               (obs_data$Di == 0))) / n
  
  return(cn)
  
}

log_like_fn = function(OMEGA){
  
  n = nrow(obs_data)
  
  cnts =
    c(mapply(dn, u = df.OMEGA.1$u, v = df.OMEGA.1$v),
      mapply(cn, u = df.OMEGA.2$u, v = df.OMEGA.2$v))
  
  l.dat = ifelse(OMEGA == 0, 0, log(OMEGA))
  
  return( sum( n * cnts * l.dat ) )
  
}

#note THETA can be one dimensional 
f_X = function(u, THETA){
  
  U = THETA[(minU - Delta):(maxU - Delta)]
  
  vec_idx = u - Delta
  
  if( ((Delta + 1) <= u) & (u <= (xi)) ){
    return( U[vec_idx] )
  }
  else{
    return(0)
  }
  
}


g_Y = function(v, THETA){
  
  if( ((Delta + 1) <= v) & (v <= (Delta + m)) ){
    vec_idx = (v - Delta) + (maxU - minU) + 1
    return( THETA[vec_idx] )
  }
  
}

alpha = function(THETA){
  
  res = c()
  for(u in c((Delta + 1):(xi))){
    v_end = min(u + 1 - Delta, Delta + m + 1)
    g_sum = sapply(c((Delta + 1):(min(u, Delta + m))), g_Y, THETA = THETA)
    res = append(res,
                 f_X(u, THETA) * sum(g_sum))
    
  }
  
  return(sum(res))
  
}

h_star = function(u , v, THETA){
  
  if( ((Delta + 1) <= u) &
      (u <= xi) &
      ((Delta + 1) <= v) &
      (v <= (Delta + m)) &
      (v <= u)){
    
    return(
      (f_X(u, THETA) * g_Y(v, THETA)) / alpha(THETA)
    )
  }
  else{
    return(0)
  }
  
}

h_bar_star = function(u, v, THETA){
  
  if( ((Delta + 1) <= u) &
      (u <= xi) &
      ((Delta + 1) <= v) &
      (v <= (Delta + m)) &
      (v <= u)){
    
    return(
      (S_X(u + 1, THETA) * g_Y(v, THETA)) / alpha(THETA)
    )
  }
  else{
    return(0)
  }
  
}

S_X = function(u, THETA){
  
  if( ((Delta + 1) <= u) & (u <= (xi)) ){
    return( sum(sapply(c(u:xi), f_X, THETA)) )
  }
  else{
    return(0)
  }
  
}


log_like_fn_0 = function(THETA){
  
  dat.1 = obs_data[obs_data$Di == 1,]
  dat.2 = obs_data[obs_data$Di == 0,]
  
  l1 = c()
  for(i in c(1:nrow(dat.1))){
    l1 = append(l1, h_star(u = dat.1$Zi[i],
                           v = dat.1$Yi[i],
                           THETA = THETA))
  }
  
  if(nrow(dat.2) == 0){
    l2 = 0
  }
  else{
    l2 = c()
    for(i in c(1:nrow(dat.2))){
      l2 = append(l2, h_bar_star(u = dat.2$Zi[i],
                                 v = dat.2$Yi[i],
                                 THETA = THETA))
    }
  }
  #l2 = c()
  #for(i in c(1:nrow(dat.2))){
  #  l2 = append(l2, h_bar_star(u = dat.2$Zi[i],
  #                             v = dat.2$Yi[i],
  #                             THETA = THETA))
  #}
  
  l1 = ifelse(l1 == 0, 0, log(l1))
  l2 = ifelse(l2 == 0, 0, log(l2))
  
  return( sum(l1) + sum(l2) )
  
}

#haz estimator functions from IME paper
Unx <- function(x) {
  ind_x <- ifelse( (obs_data$Yi <= x) & (x <= obs_data$Zi),1,0 )
  return( (1/nrow(obs_data)) * (sum(ind_x)) )
}

fnx <- function(x) {
  ind_x <- ifelse( (obs_data$Zi == x) & (obs_data$Di == 1),1,0 )
  return( (1/nrow(obs_data)) * (sum(ind_x)) )
}

lnx <- function(x) {
  
  return( fnx(x) / Unx(x) )
  
}


gnx <- function(x) {
  return(sum(obs_data$Yi == x) / nrow(obs_data))
}

f_est = function(u){
  
  idx = u - Delta
  
  if( u == (Delta + 1)){
    return( haz_est[idx] )
  }
  if( ((Delta + 1) < u) & (u <= xi ) ){
    return(haz_est[idx] * prod(1 - haz_est[1:(idx - 1)]))
  }
  
}

g_est = function(g){
  
  A = (sum(obs_data$Yi == g)/nrow(obs_data)) / 
    sum(sapply(c(g:xi), f_est))
  
  B = c()
  for(k in c((Delta + 1):(Delta + m))){
    B = append(B,
               (sum(obs_data$Yi == k)/nrow(obs_data)) / 
                 sum(sapply(c(k:xi), f_est)))
  }
  
  return( A / (sum(B)) )
  
}