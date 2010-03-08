################################################
#
#  This program "minmax.ordcorr "checks the first condition of the feasibility
#  of a correlation matrix of ordinal random numbers. A mean vector (as list)
#  needs to be specified. It returns yes/no if also a correlation matrix
#  was given and in either case the Min-Max Correlation Matrix, which has
#  the minimum correlation in the lower triangular matrix and the maximum
#  correlation in the upper triangular matrix.
#
################################################

minmax.ordcorr=function(probs,Cor=0,n=1000000,showX=FALSE)
{
  q=length(probs)
  if (length(Cor)==1) Cor=diag(1,q)
  categ_probs=0
  X=Y=matrix(n,q,data=0)
  order_probs=0
  cumul_probs=list(0)
  rounded=FALSE
  for (i in 1:q)
  {
    categ_probs[i]=length(probs[[i]])
    order_probs=order(probs[[i]]*n-trunc(probs[[i]]*n),decreasing=TRUE)
    probs[[i]]=trunc(probs[[i]]*n,0)
    if ((n-sum(probs[[i]]))>0)
    {
      rounded=TRUE
      for (j in 1:(n-sum(probs[[i]])))
      {
        probs[[i]][order_probs[j]]=probs[[i]][order_probs[j]]+1
      }
    }
    cumul_probs[[i]]=c(0,cumsum(probs[[i]]))
  }
  if (rounded) cat("Warning: Decimals of the means will be rounded after",trunc(log10(n)),"digits! \r\n")
  for (i in 1:q)
  {
    for (j in 1:categ_probs[i])
    {
      if (probs[[i]][j]>0) X[(cumul_probs[[i]][j]+1):cumul_probs[[i]][j+1],i]=j
    }
    Y[,i]=sort(X[,i],decreasing=TRUE)
    probs[[i]]=probs[[i]]/n
  }
  Cor_minmax=diag(1,q)
  yes=TRUE
  for (i in 1:(q-1))
  {
    for (j in (i+1):q)
    {
      if ((0==var(X[,i])) | (0==var(X[,j]))) { Cor_minmax[i,j]=Cor_minmax[j,i]=0 }
      else
      {
        Cor_minmax[i,j]=cor(X[,i],X[,j])
        Cor_minmax[j,i]=cor(X[,i],Y[,j])
      }
      if (length(Cor)>1)
      {
        if ((Cor_minmax[i,j]-Cor[i,j])<0) {yes=FALSE; cat("(",i,",",j,") not feasible (choose lower) \r\n")}
        if ((Cor_minmax[j,i]-Cor[i,j])>0) {yes=FALSE; cat("(",i,",",j,") not feasible (choose higher) \r\n")}
      }
    }
  }
  if (showX)
  {
    returned=list(Cor_minmax,X,yes,probs)
  }
  else
  {
    returned=list(Cor_minmax)
  }
  if (length(Cor)>1)
  {
    if (yes)
    {
      if (showX==FALSE)
      {
        cat("Warning: Only first conditions are checked. ")
        cat("Matrix can still be unfeasible! \r\n");
      }
      cat("First conditions on the specified matrix are fullfilled! \r\n")
    }
    else cat("First conditions on the specified matrix are NOT fullfilled!!! \r\n")
  }
  return(returned)
}



