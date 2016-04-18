Implementation of Least Angle Regression
-------
```{r}
#alpha is called learning rate, and a small value suggested,such as 0.05
lar<-function(x,y,alpha,tol=10^(-8))        
{
#standardize predictors with mean zero and std one or norm one
  r=y-mean(y)   
  r0=r   
  scale_x=scale(x,center=T,scale=T)
#initialize coefficients, correlation, active_set and non-active_set
  beta=matrix(rep(0,ncol(x)))
  rou=abs(cor(r,scale_x))
  active_set=which(rou==max(rou))
  non_active=(1:ncol(x))[-active_set]
  rou=max(rou)
#generate a variable step to record the iteration steps  
  step=0
  while(rou>tol)
  {
#if one variable catches up with the abosulte correaltion with current residual, put it into the active_set
#and remove form non_acitve_set
       if (any(abs(cor(r,scale_x[,non_active]))>=rou))
         {
           index=which(abs(cor(r,scale_x[,non_active]))>=rou)
           active_set=rbind(matrix(active_set),non_active[index])
           non_active=non_active[-index]
          }
#derive the direction of coefficients in active_set, and update resdiual and correaltion
    x_a=as.matrix(scale_x[,active_set])
    beta[active_set]=beta[active_set]+alpha*(solve(t(x_a)%*%x_a)%*%t(x_a)%*%r)
    r=r0-x_a%*%beta[active_set] 
    rou=as.numeric(abs(cor(r,scale_x[,active_set[1]])))
    step=step+1
  }
#return the coefficients with respect to original data(not standardized)
  beta=beta/apply(x,2,sd)
  beta0=as.matrix(mean(y)-apply(x,2,mean)%*%beta)
  #return(rbind(beta0,beta))
  return(list(beta=rbind(beta0,beta),step=step))
}
```
