genW = function(n, type="rook"){
  k = ceiling(sqrt(n))
  W = matrix(nrow=n, ncol=n, rep(0,n*n))
  for (i in 1:n){
    h = (i-1) %/% k + 1
    v = (i-1) %% k + 1
      for (dv in c(-1,0,1)){
        for (dh in c(-1,0,1)){
          diff = abs(dv) + abs(dh)
          v1 = v + dv
          h1 = h + dh
          neib = (diff>0) & (v1 >0) & (h1 >0) & (v1 <=k) & (h1 <=k) & (type=="queen" || diff==1)
          if(neib){
            j  = v1 + (h1-1)*k
            if (j<=n) W[i,j] = 1
          }
        }
      }
  }
  return(W)
}

row.stdrt = function(W){
  for (j in 1:nrow(W)){
    W[j,] = W[j,]/sum(W[j,])
  }
  return(W)
}

W = row.stdrt(genW(100, type="queen"))
log(det(I - rho * W))
sum(log(1-rho*eigen(W,only.values=TRUE)$values))