alpha = function(U, G){
  
  res = c()
  for(u in c((Delta + 1):(omega))){
    g_sum = sapply(c((Delta + 1):(min(u, Delta + m))), g_Y, G = G)
    res = append(res,
                 f_X(u, U) * sum(g_sum))
    
  }
  
  return(sum(res))
  
}

f_X = function(u, U){
  
  vec_idx = u - Delta
  
  if( ((Delta + 1) <= u) & (u <= (omega)) ){
    return( U[vec_idx] )
  }
  else{
    return(0)
  }
  
}

g_Y = function(v, G){
  
  if( ((Delta + 1) <= v) & (v <= (Delta + m)) ){
    vec_idx = v - Delta
    return( G[vec_idx] )
  }
  else{
    return(0)
  }
  
}

#haz estimator functions from E&S paper
Cnx <- function(x) {
  ind_x <- ifelse( (obs_data$Yi <= x) & (x <= obs_data$Xi),1,0 )
  return( (1/nrow(obs_data)) * (sum(ind_x)) )
}

lnx <- function(x) {
  return( (1/nrow(obs_data)) * (1/Cnx(x)) * (sum(obs_data$Xi == x)) )
}

fnx <- function(x) {
  return(sum(obs_data$Xi == x) / nrow(obs_data))
}

bnx <- function(y) {
  return( (1/nrow(obs_data)) * (1/Cnx(y)) * (sum(obs_data$Yi == y)) )
}

gnx <- function(x) {
  return(sum(obs_data$Yi == x) / nrow(obs_data))
}

f_est = function(u){
  
  idx = u - Delta
  
  if( u == (Delta + 1)){
    return( haz_est[idx] )
  }
  if( ((Delta + 1) < u) & (u <= omega ) ){
    return(haz_est[idx] * prod(1 - haz_est[1:(idx - 1)]))
  }
  
}

g_est = function(g){
  
  idx = g - Delta
  
  if( g == (Delta + m)){
    return( rev_haz_est[idx] )
  }
  if( ((Delta + 1) <= g) & (g < Delta + m ) ){
    return(rev_haz_est[idx] * prod(1 - rev_haz_est[(idx + 1):(m)]))
  }
  
}

log_like_fn_0 = function(U, G){
  
  alp = alpha(U, G)
  
  n = nrow(obs_data)
  
  Li = c()
  for(k in c((Delta + 1):(Delta + m))){
    for(j in c(k:(omega))){
      cnt = sum(( (obs_data$Yi == k ) & (obs_data$Xi == j) ))
      val = ( f_X(j, U) * g_Y(k, G) ) / alpha(U, G)
      #if( (cnt > 0) & (val == 0) ){print(cnt); print(val)}
      if( (val == 0) ){
        Li = append(Li, 0 )
      }
      if( (val > 0) )
        Li = append(Li, cnt * log( val ) )
    }
  }
  
  return( sum(Li) )
  
}

# u = c( (Delta + 1) : omega)
# reps = length(c( (Delta + 1) : (Delta + m)))
# 
# x_col = c()
# y_col = c()
# for(u in c( (Delta + 1) : omega)){
#   
#   for(v in c( (Delta + 1) : (Delta + m))){
#     if(v <= u){
#       x_col = append(x_col, u)
#       y_col = append(y_col, v)
#     }
#   }
#   
# }
# 
# h_inv = data.frame("X" = x_col,
#                    "Y" = y_col)

log_like_fn = function(THETA){
  
  n = nrow(obs_data)
  
  Li = c()
  for(j in c((Delta + 1):(omega))){
    for(k in c((Delta + 1):min(j, Delta + m))){
      #print(c(j,k))
      cnt = sum(( (obs_data$Yi == k ) & (obs_data$Xi == j) ))
      val = THETA[which( (h_inv$X == j) & (h_inv$Y == k) )]
      val.adj = ifelse(val == 0, 0, log(val))
      Li = append(Li, cnt * val.adj )
    }
  }
  
  return( sum(Li) )
  
}

#results.M.no.zero = results.M[which(results.M[,"n.0"] > 0),]

#results.M[which( (results.M[,1] > 0) & (results.M[,2] < 0) ),]
