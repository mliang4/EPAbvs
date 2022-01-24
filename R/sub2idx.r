sub2idx = function(subj,ni){
  idx = NULL
  for(i in subj){
    idx = c(idx, ((i-1)*ni+1):(i*ni))
  }
  return(idx)
}

