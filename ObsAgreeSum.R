
### Calculate observed agreement 
ObsAgreeSum <- function(y1, y2){
  
  # 0 = negative agreement
  # 1 = disagreement
  # 2 = positive agreement
  counts = data.frame(table((y1+y2)))
  rownames(counts) = counts[,1]
  
  if(is.na(counts["0",2])){
    neg.agree = 0
  }else{
    neg.agree = counts["0",2] / sum(counts$Freq)
  }
  
  if(is.na(counts["1",2])){
    dis = 0
  }else{
    dis = counts["1",2] /sum(counts$Freq)
  }
  
  if(is.na(counts["2",2])){
    pos.agree = 0
  }else{
    pos.agree = counts["2",2] / sum(counts$Freq)
  }
  agree = pos.agree + neg.agree
  
  obsProbAgree = c(agree,dis,sum(counts$Freq))
  return(obsProbAgree)
}
