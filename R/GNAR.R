#' @title A
#' @name A
#' @description A 
NULL

#' @title Z
#' @name Z
#' @description Z 
NULL

#' @title data
#' @name data
#' @description data 
NULL



#' @title Grouped Network Autoregression using R.
#' @description Grouped Network Autoregression  by two steps estimation using R.
#' @import igraph MASS plyr methods quantmod tseries lmtest vars forecast
#' @param A The adjacency matrix of the data. (numeric)
#' @param Z The  node specific covariates Z should be a vector or matrix.(numeric)
#' @param data The data matrix. (numeric)
#' @param begin The begin day of training set. (numeric)
#' @param end The end day of training set. (numeric)
#' @param K the number of clusters. (numeric)
#' @param t the number of days to forecast. (numeric)
#' @return the forecast and estimation also with the accuracy of the forecast
#' @examples
#' \dontrun{
#' data(A)
#' data(Z)
#' data(data)
#' begin=1; end=503
#' K=2; t=10
#' GNAR(A,Z,data,begin,end,K,t)
#' }
#' @export








GNAR<-function(A,Z,data,begin,end,K,t)
{       
  ### Two step estimation
  Cluster.NAR<-function(Ymat, W,  Z, K, method = "complete",  group = NULL)
  {
    N = nrow(Ymat)
    Time = ncol(Ymat)
    Ymat1 = W%*%Ymat     
    ### 求得回顾数据集
    if (is.null(Z))
      yy_dat = data.frame(id = rep(1:N, each = Time - 1),
                          inter = 1,
                          net = as.vector(t(Ymat1[,-ncol(Ymat)])),
                          lagY = as.vector(t(Ymat[,-ncol(Ymat)])),
                          Y = as.vector(t(Ymat[,-1]))) 
    else
      yy_dat = data.frame(id = rep(1:N, each = Time - 1),
                          inter = 1,
                          net = as.vector(t(Ymat1[,-ncol(Ymat)])),
                          lagY = as.vector(t(Ymat[,-ncol(Ymat)])),
                          Z[rep(1:N, each = Time-1),],
                          Y = as.vector(t(Ymat[,-1])))
    
    ### 对每个节点求b_i估计值
    paras = ddply(yy_dat, .(id), function(x){
      X = as.matrix(x[,2:4])
      invXX = ginv(crossprod(X))                                                                                   ### the response vector
      thetaEst = invXX%*%colSums(X*x$Y)   
      df = data.frame(matrix(c(thetaEst), nrow = 1))
    })
    
    colnames(paras)[-1] = c("intercept", "network", "momentum")
    
    ### 规范化网络和动量参数并计算节点的距离
    para_scale = apply(paras[,3:4], 2, function(x) x/max(abs(x)))#scale(paras[,3:4])
    para_dist = dist(as.matrix(para_scale))
    
    ### 实施聚类算法
    if (method=="kmeans")
    {
      #nar_para = betaOLS(Ymat, W, Z)
      #ini_theta = (sapply(nar_para$theta[2:3], function(x) runif(K, min = x-0.05, max = x+0.05)))
      k_res = kmeans(para_scale, K)
      memb = k_res$cluster
    }
    else
    {
      hcl = hclust(para_dist, method = method)
      memb = cutree(hcl, k = K)
    }
    alpha = table(memb)/length(memb)
    yy_dat$group = rep(memb, each = Time - 1) # 求得不同节点的组
    ### 对每个组重新进行估计
    theta_est = ddply(yy_dat, .(group), function(x){
      X = as.matrix(x[,2:(ncol(x)-2)])
      invXX = ginv(crossprod(X))                                                                                   ### the response vector
      thetaEst = invXX%*%colSums(X*x$Y)   
      df = data.frame(matrix(c(thetaEst), nrow = 1))
    })
    return(list(theta = t(theta_est)[-1,], alpha = alpha, group = memb))
  }
  N=length(data)
  L=end-begin+1
  data1=matrix(0,nrow=N,ncol=L)
  for(i in 1:N) data1[i,]=data[[i]][begin:end]
  for(j in 1:t) data1=cbind(data1,0)
  W=A/rowSums(A)
  W[is.na(W)]=0
  Z=as.matrix(Z)
  cl_res = Cluster.NAR(data1, W, Z, K, method = "complete") # two-step estimation
  theta=cl_res$theta
  group=cl_res$group
  D=NULL
  before=NULL  #最后一个数据的值
  prob=NULL
  forecast=matrix(0,nrow=N,ncol=t)  #预测值
  real=matrix(0,nrow=N,ncol=t)  #真实值
  D=matrix(0,nrow=N,ncol=t)     #是否正确 
  for(i in 1:N)before[i]=exp(data[[i]][end])
  for(j in 1:t){
    Ymat=data1[,1:(L+j-1)]
    WYmat = W%*%Ymat
    N = nrow(Ymat)
    Time = ncol(Ymat)
    if (is.null(Z))
      X = cbind(1, lagWY = as.vector(WYmat), lagY = as.vector(Ymat))
    else
      X = cbind(1, lagWY = as.vector(WYmat), lagY = as.vector(Ymat),
                do.call(rbind, rep(list(Z), Time)))
    Yhat = matrix((X%*%theta)[cbind(1:nrow(X), rep(group, Time))], nrow = N)
    data1[,L+j]=Yhat[,L+j-1]
  }
  for(j in 1:t)  forecast[,j]=data1[,L+j]
  for(j in 1:t){
    for(i in 1:N){
      real[i,j]=exp(data[[i]][end+j])
      if((real[i,j]-before[i])*(forecast[i,j]-before[i])>0)
        D[i,j]=1
      if((real[i,j]-before[i])*(forecast[i,j]-before[i])<=0)
        D[i,j]=0
    }
    prob[j]=sum(D[,j])/length(D[,j])
  }
  M=cbind(before,real,forecast,D)
  colnames(M)=c("before",rep("real",t),rep("forecast",t),rep("D",t))
  rownames(M)=names(A[,1])
  return(list(M,prob))
}
