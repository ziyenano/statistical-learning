Lasso implemented by Lars modification
---
```{r}
#hence this program is very similar to lars,only a snippet added
#修改Lars实现lasso
#here variable shrink is penalized tuning parameters,which ranges from  0 to 1.if set 1 the reslut is 
#equivalent to olse,and if set 0,all the slope equals zero.
#参数shrink是一个惩罚调节参数，从0到1，如果取值为1，结果等同于最小二乘，如果取0，所有斜率等于0
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
    #if 1 norm of current coefficients  divide 1 norm of full olse great than shrink,then break the loop
    #which is with respect to lasso
    #如果当前斜率的1范数除以完全最小二乘估计大于shrink，跳出循环，这与lasso相对应
    if (sum(abs(beta))/sum(abs(full_olse))>=shrink) break
    if (any(abs(cor(r,scale_x[,non_active]))>=rou))
    {
      index=which(abs(cor(r,scale_x[,non_active]))>=rou)
      active_set=rbind(matrix(active_set),non_active[index])
      non_active=non_active[-index]
    }
    x_a=as.matrix(scale_x[,active_set])
```
```{r}
   **lasso modification**
    beta_pre=beta[active_set]
    beta[active_set]=beta[active_set]+alpha*(solve(t(x_a)%*%x_a)%*%t(x_a)%*%r)
    #if a nonzero coefficient cross zero,remove from the active set and recopute the olse direction
    #如果一个非零的系数穿过0点，将它从活跃集删除，并重新计算最小二乘的方向
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
