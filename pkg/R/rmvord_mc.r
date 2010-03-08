############################################
#
#  this function produces n multivariate ordinal random numbers
#  via monte-carlo simulation
#
############################################


rmvord_mc=function(probs,Cor,n=1)
{
    q=length(probs)
    categ_probs=0
    cumul_probs=list(0)
    quant_probs=list(0)
    for (i in 1:q)
    {
      categ_probs[i]=length(probs[[i]])
      cumul_probs[[i]]=cumsum(1:categ_probs[i]/10^12+probs[[i]])
      cumul_probs[[i]][categ_probs[i]]=1   # otherwise qnorm computes NA instead of inf
      quant_probs[[i]]=qnorm(p=cumul_probs[[i]],mean=0,sd=1)
    }
    Cor_norm=Cor
    for (i in 1:(q-1))
    {
      for (j in (i+1):q)
      {
        gridd=rep(0,times=21)
        for (steps in 1:21)
        {
          norm=rmvnorm(n=4000,sigma=matrix(2,2,data=c(1,(steps-11)/10,(steps-11)/10,1)))
          norm[,1]=cut(x=norm[,1],breaks=c(-1/0,quant_probs[[i]]),right=FALSE)
          norm[,2]=cut(x=norm[,2],breaks=c(-1/0,quant_probs[[j]]),right=FALSE)
          gridd[steps]=cor(norm[,1],norm[,2])
        }
        f=line(cbind(gridd,(1:21-11)/10))
        cor_opt=coef(f)[2]*Cor[i,j]+coef(f)[1]
        from=cor_opt-0.05
        if (from < -1) from=-1
        if (from + 0.1 > +1) from=0.9
        for (steps in 1:21)
        {
          norm=rmvnorm(n=10000,sigma=matrix(2,2,data=c(1,from+(steps-1)/200,from+(steps-1)/200,1)))
          norm[,1]=cut(x=norm[,1],breaks=c(-1/0,quant_probs[[i]]),right=FALSE)
          norm[,2]=cut(x=norm[,2],breaks=c(-1/0,quant_probs[[j]]),right=FALSE)
          gridd[steps]=cor(norm[,1],norm[,2])
        }
        f=line(cbind(gridd,from+(1:21-1)/200))
        cor_opt=coef(f)[2]*Cor[i,j]+coef(f)[1]
        from=cor_opt-0.005
        if (from < -1) from=-1
        if (from + 0.01 > +1) from=0.99
        for (steps in 1:21)
        {
          norm=rmvnorm(n=60000,sigma=matrix(2,2,data=c(1,from+(steps-1)/2000,from+(steps-1)/2000,1)))
          norm[,1]=cut(x=norm[,1],breaks=c(-1/0,quant_probs[[i]]),right=FALSE)
          norm[,2]=cut(x=norm[,2],breaks=c(-1/0,quant_probs[[j]]),right=FALSE)
          gridd[steps]=cor(norm[,1],norm[,2])
        }
        f=line(cbind(gridd,from+(1:21-1)/2000))
        cor_opt=coef(f)[2]*Cor[i,j]+coef(f)[1]
        Cor_norm[i,j]=Cor_norm[j,i]=cor_opt
      }
    }
    retval=rmvnorm(n=n,sigma=Cor_norm)
    for (i in 1:q)
    {
     retval[,i]=cut(x=retval[,i],breaks=c(-1/0,quant_probs[[i]]),right=FALSE)
    }
    retval
}

