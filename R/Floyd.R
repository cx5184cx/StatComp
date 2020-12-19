#' @title Using Floyd algorithm to find shortest path.
#' @description Floyd algorithm to find shortest path using R.
#' @importFrom stats cutree dist hclust kmeans
#' @param A The graph with weight.(numeric)
#' @return the shortest cost and the path.
#' @examples
#' \dontrun{
#' a<-matrix(0,7,7)
#' a[1,2]=3;a[1,3]=5; a[2,3]=1;a[2,4]=5;a[2,5]=8
#' a[3,4]=7;a[3,5]=4;a[3,6]=10; a[4,5]=3;a[4,7]=6
#' a[5,6]=1;a[5,7]=2; a[6,7]=2
#' b<-a+t(a)
#' floyd(b)
#' }
#' @export

####定义Floyd函数
floyd<-function(A){
  A[A==0]<-Inf
  n<-nrow(A)
  D<-A
  path<-matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      if(is.finite(D[i,j])==T){path[i,j]=j}
    }
  }
  for(k in 1:n){
    for(i in 1:n){
      for(j in 1:n){
        if(D[i,k]+D[k,j]<D[i,j]){
          D[i,j]=D[i,k]+D[k,j];
          path[i,j]=path[i,k]
        }
      }
    }
  }
  return(list(D=D,path=path))
}
