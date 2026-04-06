term_time = function(mbalt_row){
  
  init_age = interval(my(mbalt_row$originationDate), mdy(mbalt_row$reportingPeriodBeginDate)) %/% months(1) + 1
  
  term_vec <- vector()
  for (i in c(1:(e-M-delta))){
    term_vec = append(term_vec, mbalt_row[,paste("TI",i,sep="")])
  }
  t_time = which(!is.na(term_vec))
  if (identical(t_time,integer(0))){
    X = length(term_vec) + init_age - 1
    C = 0
  }
  else {
    X = t_time + init_age - 1
    C = 1
  }
  return(c(X,C))
}