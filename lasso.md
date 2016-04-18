Lasso implemented by Lars modification
---
```{r}
#this program is very similar to lars, only a snippet added.
#penalized tuning parameters shrink ranging from  0 to 1.
#if  one, it is equivalent to OLSE, and if zero, all parameters are set to zero.
lasso<-function(x,y,alpha,shrink,tol=10^(-8))        
{
  r=y-mean(y)
  r0=r
  scale_x=scale(x,center=T,scale=T)
  beta=matrix(rep(0,ncol(x)))
  rou=abs(cor(r,scale_x))
  active_set=which(rou==max(rou))
  non_active=(1:ncol(x))[-active_set]
  rou=max(rou)
  full_olse=solve(t(scale_x)%*%scale_x)%*%t(scale_x)%*%r0
  step=0
  while(rou>tol)
  { 
    # break the loop, if L1 norm of current coefficients dividing that of OLSE is great than shrink
    if (sum(abs(beta))/sum(abs(full_olse))>=shrink) break
    if (any(abs(cor(r,scale_x[,non_active]))>=rou))
    {
      index=which(abs(cor(r,scale_x[,non_active]))>=rou)
      active_set=rbind(matrix(active_set),non_active[index])
      non_active=non_active[-index]
    }
    x_a=as.matrix(scale_x[,active_set])
```
**lasso modification**
```{r}
    beta_pre=beta[active_set]
    beta[active_set]=beta[active_set]+alpha*(solve(t(x_a)%*%x_a)%*%t(x_a)%*%r)
    #once a nonzero coefficient touches zero, remove it from the active set and recompute the OLSE direction
    if (any((beta[active_set]*beta_pre)<=0 & beta_pre !=0))  
    {
      m_index=which((beta[active_set]*beta_pre)<=0 & beta_pre !=0)
      active_set=active_set[-m_index]
      non_active=rbind(matrix(non_active),matrix(active_set[m_index]))
    }
    x_a=as.matrix(scale_x[,active_set])
    beta[active_set]=beta[active_set]+alpha*(solve(t(x_a)%*%x_a)%*%t(x_a)%*%r)
```
**end modification**
```{r}
    r=r0-x_a%*%beta[active_set] 
    rou=as.numeric(abs(cor(r,scale_x[,active_set[1]])))
    step=step+1
  }
  beta=beta/apply(x,2,sd)
  beta0=as.matrix(mean(y)-apply(x,2,mean)%*%beta)
  return(list(beta=rbind(beta0,beta),step=step))
}
```
