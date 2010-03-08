################################
#
# This program "check.ordcorr" checks if the stated input for the correlations
# fits with the given marginal probabilities in that way that
# the correlation matrix is feasible. It returns a list with the rounded
# (if neccesary) means, the approximated correlation matrix and the common probabilities
# as output. This output will be a result that is "closest" to the inputs of
# user and is therefore a usefull alternative to his specifications. It also judges
# the entries of the user's correlation matrix to be either
# feasible or not. The function stops if more switches do not
# decrease the distance to the specified correlation matrix or a maximum number
# of switches is reached.
#
################################

check.ordcorr=function(probs,Cor,switchmax=1e+06,n=100,showcommonprob=FALSE,showX=FALSE)
#suppressWarnings(
{
  test=minmax.ordcorr(probs=probs,Cor=Cor,n=n,showX=TRUE)
  X=test[[2]]
  probs_found=test[[4]]                # matters if means were rounded
  q=length(probs)
  delta_sumls=array(data=1/0,dim=c(n,n,q))
  means=vars=0
  for (i in 1:q)
  {
    means[i]=mean(X[,i])
    vars[i]=var(X[,i])*(n-1)/n
  }
  means_prod=means%*%t(means)
  vars_prod=vars%*%t(vars)
  comm_mean_opt=Cor*sqrt(vars_prod)+means_prod
  update_delta_sumls=function(X,delta_sumls,vars_prod,comm_mean_opt,q,n)
  {
    X_diff_help=matrix(data=0,n,n)
    X_diff=array(data=0,dim=c(n,n,q))
    c_=array(data=0,dim=c(n,n,q,q))
    comm_mean=t(X)%*%X/n
    for (g in 1:q)
    {
      X_diff_help[,1:n]=X[,g]
      X_diff[,,g]=X_diff_help-t(X_diff_help)
    }
    for (g in 1:(q-1))
    {
      for (j in (g+1):q)
      {
        c_[,,g,j]=c_[,,j,g]=-1/n*X_diff[,,g]*X_diff[,,j]
      }
    }
    delta_sumls[,,]=0
    for (g in 1:q)
    {
      for (j in 1:q)
      {
        if (j!=g) {delta_sumls[,,g]=delta_sumls[,,g]+2*c_[,,g,j]/vars_prod[g,j]*(comm_mean[g,j]-comm_mean_opt[g,j]+c_[,,g,j]/2) }
      }
    }
    delta_sumls[is.na(delta_sumls)]=0
    2*delta_sumls
  }
  for (count in 1:switchmax)
  {
    delta_sumls=update_delta_sumls(X,delta_sumls,vars_prod,comm_mean_opt,q,n)
    min_pos=which.min(delta_sumls)
    g_s=(min_pos-1) %/% (n^2) + 1
    w_s=(((min_pos-1) %% n^2) %/% n) + 1
    v_s=((min_pos-1) %% n^2) %% n +1
    if (delta_sumls[v_s,w_s,g_s]<0)
    {
      helpv=X[v_s,g_s]             # update X
      X[v_s,g_s]=X[w_s,g_s]
      X[w_s,g_s]=helpv
    } else        # no improvments by switching
    {
      cat("Converged after ",count-1," iterations. \r\n")
      break
    }
  }
  if (count==switchmax) {cat("Procedure interrupted! Maximum of",switchmax,"switches reached! \r\n")}
  Cor_found=cor(X)
  Cor_found[is.na(Cor_found)]=0
  yes=TRUE
  comm_mean=t(X)%*%X/n
  for (i in 1:(q-1))
  {
    for (j in (i+1):q)
    {
      if ((comm_mean[i,j]-comm_mean_opt[i,j])>(q-1)/(2*n)) {yes=FALSE; cat("(",i,",",j,") not feasible (choose higher) \r\n")}
      if ((comm_mean_opt[i,j]-comm_mean[i,j])>(q-1)/(2*n)) {yes=FALSE; cat("(",i,",",j,") not feasible (choose lower) \r\n")}

    }
  }
  if (yes) cat("Stated matrix feasible! \n \n")
  else cat("Stated matrix NOT feasible!!! Try higher n! \n \n")
  if (showcommonprob==TRUE)
  {
    commonprob=array(dim=c(q,q,100,100),data=0)
    for (i in 1:(q-1))
    {
      for (j in (i+1):q)
      {
        for (a in 1:length(probs[[i]]))
        {
          for (b in 1:length(probs[[j]]))
          {
            commonprob[i,j,a,b]=commonprob[j,i,a,b]=sum(X[,i]==a %*% t(X[,j]==b))/n
          }
        }
      }
    }
    if (showX==TRUE) return(list(probs_found,Cor_found,commonprob,X))
    else return(list(probs_found,Cor_found,commonprob))
  } else
  {
    if (showX==TRUE) return(list(probs_found,Cor_found,X))
    else return(list(probs_found,Cor_found))
  }
}#)
