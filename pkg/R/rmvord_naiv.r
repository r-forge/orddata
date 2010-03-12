############################################
#
#  this function produces n ordinal random numbers
#  via the native method
#
############################################

rmvord_naiv=function(probs,Cor,n=1,showCor=TRUE)
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
    }
    if(showCor) print(Cor)
    retval=rmvnorm(n=n,sigma=Cor)
    for (i in 1:q)
    {
     retval[,i]=cut(x=retval[,i],breaks=c(-1/0,quant_probs[[i]]),right=FALSE)
    }
    retval
}
