########################################################
#
#  this function generates n random numbers of ordinal distributed
#  variables. It uses the correlation of the ordinal variates (in
#  a matrix) as well as the marginal probs (in a list) as input.
#  The method uses is binarisation.
#
########################################################

rmvord_b=function(n=1,probs,Cor,showCor_b=FALSE)
{
  q=length(probs)
  enum=vector(length=q)
  enum[1:q]=0
  k=0
  margprobs=0
  var_bin=0
  var_ord=0
  m=0
  for (i in 1:q)
  {
    k[i]=length(probs[[i]])
    margprobs[i]=weighted.mean(x=0:(k[i]-1)/(k[i]-1),w=probs[[i]])
    var_bin[i]=margprobs[i]*(1-margprobs[i])
    var_ord[i]=sum(((1:k[i]-1)/(k[i]-1)-margprobs[i])^2*probs[[i]])
    if (k[i]>2)
    {
      m[i]=var_ord[i]/var_bin[i]
    } else m[i]=1
  }
  Cor_b=Cor
  for (i in 1:(q-1))
  {
    for (j in (i+1):q)
    {
      Cor_b[i,j]=Cor_b[j,i]=Cor_b[i,j]/sqrt(m[i]*m[j])
    }
  }
  if (showCor_b) print(Cor_b)
  random_ord=-rmvbin(n=n,margprob=margprobs,bincorr=Cor_b)
  for (i in 1:q)
  {
    if (k[i]>2)
    {
      count=-sum(random_ord[,i])
      random_unif=runif(n=n)
      breakss=vector(length=k[i])
      breakss[1]=0
      breakss[2:k[i]]=probs[[i]][1:(k[i]-1)]*(k[i]-(1:(k[i]-1)))/(k[i]-1)/(1-margprobs[i])
      breakss=1:k[i]/10^12+cumsum(breakss)
      random_ord[random_ord[,i]==0,i]=cut(x=random_unif[1:(n-count)],breaks=breakss,labels=FALSE)
      breakss[2:k[i]]=probs[[i]][2:k[i]]*(1:(k[i]-1)/(k[i]-1))/margprobs[i]
      breakss=1:k[i]/10^12+cumsum(breakss)
      random_ord[random_ord[,i]==-1,i]=cut(x=random_unif[(n-count+1):n],breaks=breakss,labels=FALSE)+1
    } else random_ord[,i]=-random_ord[,i]
  }
  random_ord
}






