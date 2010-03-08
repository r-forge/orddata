############################################
#
#  this function produces n ordinal random numbers
#  via the analytical method
#
############################################

rmvord=function(probs,Cor,n=1)
{
    q=length(probs)
    categ_probs=0
    cumul_probs=list(0)
    quant_probs=list(0)
    means=0
    vars=0
    var.wt=function(x,w)
    {
      m=weighted.mean(x=x,w=w)
      sum((x[1:length(x)]-m)^2*w[1:length(x)])
    }
    for (i in 1:q)
    {
      categ_probs[i]=length(probs[[i]])
      cumul_probs[[i]]=cumsum(1:categ_probs[i]/10^12+probs[[i]])
      cumul_probs[[i]][categ_probs[i]]=1
      quant_probs[[i]]=qnorm(p=cumul_probs[[i]],mean=0,sd=1)
      means[i]=weighted.mean(x=1:categ_probs[i],w=probs[[i]])
      vars[i]=var.wt(x=1:categ_probs[i],w=probs[[i]])
    }
    Cor_norm=Cor
    for (i in 1:(q-1))
    {
      for (j in (i+1):q)
      {
        gridd=rep(0,times=201)
        for (ci in 1:(categ_probs[i]-1))
        {
          for (cj in 1:(categ_probs[j]-1))
          {
            for (steps in -100:100)
            {
              gridd[101+steps]=gridd[101+steps]+pmvnorm(upper=c(quant_probs[[i]][ci],quant_probs[[j]][cj]),corr=matrix(2,2,data=c(1,steps/100,steps/100,1)))[1]
            }
          }
        }
        f=suppressWarnings(approxfun(y=-100:100/100,x=gridd))
        Cor_norm[i,j]=Cor_norm[j,i]=f(Cor[i,j]*sqrt(vars[i]*vars[j])+means[i]*means[j]-categ_probs[i]*categ_probs[j]+categ_probs[j]*sum(cumul_probs[[i]][1:(categ_probs[i]-1)])+categ_probs[i]*sum(cumul_probs[[j]][1:(categ_probs[j]-1)]))
      }
    }
    retval=rmvnorm(n=n,sigma=Cor_norm)
    for (i in 1:q)
    {
     retval[,i]=cut(x=retval[,i],breaks=c(-1/0,quant_probs[[i]]),right=FALSE)
    }
    retval
}



