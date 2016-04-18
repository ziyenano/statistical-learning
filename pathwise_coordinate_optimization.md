Lasso implemented by Pathwise_Coordinate_Optimization
-----------
```{r}
path_copt<-function(x,y,lambda,tol=1*10^-8)
{
  r=y-mean(y)
  c_x=scale(x,center=T,scale=F)
  ##standardize with mean zero and norm one
  std=apply(c_x,2,function(x) norm(matrix(x),type='f'))
  c_x=c_x%*%diag(1/std)  
  beta_old=matrix(0,ncol(x),1)
  beta=beta_old
  s=1
  eps=0.1
  while(eps>tol)
  { 
    beta_old=beta
    for (i in 1:length(beta))
    {
      ##partial residual
      p_resdiual=r-c_x[,-i]%*%beta[-i]                   
      p_beta=t(c_x[,i])%*%p_resdiual/(c_x[,i]%*%c_x[,i])  
      ##the relation between Lasso estimator and OLSE in independent case
      beta[i]=sign(p_beta)*max((abs(p_beta)-lambda),0)    
    }     
    eps=norm(beta-beta_old)
    s=s+1
  }
  print(s)
  return(beta/std)
}
```
