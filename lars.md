Implementation of Least Angle Regression
-------
```{r}
#here alpha is called learning rate,a small value is suggested,such as 0.05
#这里alpha是学习速率，建议给一个比较小的值，比如0.05
lar<-function(x,y,alpha,tol=10^(-8))        
{
#standardize predictors with mean zero and std one,also you can choose norm one
#标准化预测变量，使得均值为0，标准差为1 ，你也可以选择范数为1
  r=y-mean(y)   
  r0=r   
  scale_x=scale(x,center=T,scale=T)
#initialize coefficients,correlation,active_set and non-active_set respectively
#初始化系数,相关系数，活跃集合，与之对应的非活跃集合
  beta=matrix(rep(0,ncol(x)))
  rou=abs(cor(r,scale_x))
  active_set=which(rou==max(rou))
  non_active=(1:ncol(x))[-active_set]
  rou=max(rou)
#generate a variable step to record the iteration steps  
#生成变量step用来记录迭代的步
  step=0
  while(rou>tol)
  {
#if one variable catches up with the abosulte correaltion with current residual,put it into the active_set
#and remove form non_acitve_set
#如果一个变量赶上与当前残差相关系数的绝对值，将它加入活跃集合，并从非活跃集合中删除
       if (any(abs(cor(r,scale_x[,non_active]))>=rou))
         {
           index=which(abs(cor(r,scale_x[,non_active]))>=rou)
           active_set=rbind(matrix(active_set),non_active[index])
           non_active=non_active[-index]
          }
#calculate the direction of coefficients in active_set,the new resdiual and new correaltion
#计算活跃集合中系数前进的方向，并计算新的残差以及新的相关系数
    x_a=as.matrix(scale_x[,active_set])
    beta[active_set]=beta[active_set]+alpha*(solve(t(x_a)%*%x_a)%*%t(x_a)%*%r)
    r=r0-x_a%*%beta[active_set] 
    rou=as.numeric(abs(cor(r,scale_x[,active_set[1]])))
    step=step+1
  }
#get the coefficients with respect to original data(not standardized)
#计算原始数据对应的系数
  beta=beta/apply(x,2,sd)
  beta0=as.matrix(mean(y)-apply(x,2,mean)%*%beta)
  #return(rbind(beta0,beta))
  return(list(beta=rbind(beta0,beta),step=step))
}
```
