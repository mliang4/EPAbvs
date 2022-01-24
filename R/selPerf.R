##################################################
# Evaluate Selection Performance

# input:
#   gammaOut - a vector of posterior inference result (binary) on selection
#   trueGamma - a vector of true inclusion (binary)

##################################################

selPerf = function(estimate, truth){
  Gindex1 = which(truth==1)
  Gindex0 = which(truth==0)
  pos = which(estimate==1)
  neg = which(estimate==0)
  
  FP = sum(pos %in% Gindex0)
  FN = sum(neg %in% Gindex1)
  TP = sum(pos %in% Gindex1)
  TN = sum(neg %in% Gindex0)
  
  # fpr
  fpr = FP/(FP+TN)
  # fnr
  fnr = FN/(FN+TP)
  
  tpr = 1 - fnr
  tnr = 1 - fpr
  
  # Matthews correlation coefficient
  denom = sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN)
  
  mcc = (TP*TN - FP*FN)/denom
  
  return(list(TPR = tpr, TNR = tnr, FPR = fpr, FNR = fnr, MCC = mcc))
}