tst =function(){
n = 10
rho = 0.2
W = genW(n)
SpW = solve(diag(n) - rho*W)
sigma = t(SpW) %*% SpW

n = 2
sigma = matrix(c(1,0.9,0.9,1), n,n)
U = chol(sigma)
L = t(U)
C = solve(L)
round(C %*% sigma %*% t(C),digits=10)
round(L %*% t(L) - sigma,digits=10)

x = rep(0,n)
pmvnorm(lower=rep(-Inf,n),upper=x, sigma=sigma)
z = as.vector(solve(L) %*% x)
pmvnorm(lower=rep(-Inf,n),upper=z)


C = matrix(c(2,0,1,1), 2,2)
x = c(1,2)
pmvnorm(lower=rep(-Inf,2),upper=x, sigma=diag(2))
pmvnorm(lower=rep(-Inf,2),upper=as.vector(C%*%x), sigma=C %*% t(C))

dmvnorm(x=x, sigma=diag(n))
dmvnorm(x=as.vector(C%*%x), sigma=C %*% t(C))*sqrt(det(C %*% t(C)))

mSigma = matrix(c(3,1,1,2), 2,2)
C = chol(mSigma)


pmvnorm(lower=rep(-Inf,n),upper=sqrt(10)*x, sigma=10*diag(n))

pmvnorm(lower=rep(-Inf,n),upper=x, sigma=diag(n))
pmvnorm(lower=rep(-Inf,n),upper=2*x, sigma=4*diag(n))
pnorm(1)
pnorm(2, sd=2)
    
draw2RV = function (func.density, theta, phi){
 x1 <- seq(-5, 5, by = 0.2)
    x2 <- x1
    
    z <- outer(x1, x2, func.density)
    print(sum(z))
    nrz <- nrow(z)
    ncz <- ncol(z)
    nbcol <- 100
    jet.colors <- colorRampPalette( c("blue", "red") )
    color <- jet.colors(nbcol)
    zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    facetcol <- cut(zfacet, nbcol)
    
    persp(x1,x2, z, theta = theta , phi = phi , col=color[facetcol],xlab="X1", ylab="X2", zlab="f(x1,x2)")
}
x1 <- seq(-4, 4, by = 0.1)
x2 <- x1
normal1 <-function(x1,x2){dmvnorm(cbind(x1,x2), rep(0,2), mSigma)}
normal <-function(x1,x2){dmvnorm(t(solve(t(chol(mSigma))) %*% t(cbind(x1,x2))),rep(0,2),diag(2))*det(solve(t(chol(mSigma))))}
outer(x1, x2, normal)/outer(x1, x2, normal1)

Lx = (solve(t(chol(mSigma))) %*% t(cbind(1,2)))
x1 <- seq(-10, 1, by = 0.005)
x2 <- seq(-10, 2, by = 0.005)
sum(outer(x1, x2, normal)*0.000025)
sum(outer(x1, x2, normal1)*0.000025)

pnormal1 <-function(x1,x2){pmvnorm(lower=rep(-Inf,2), upper=c(x1,x2), rep(0,2), sigma=mSigma)}
pnormal1_grad <-function(x1,x2){dmvnorm(x=c(x1,x2), rep(0,2), sigma=mSigma)*chol(mSigma)}
pnormal <-function(x1,x2){pmvnorm(lower=rep(-Inf,2), upper=as.vector(solve(t(chol(mSigma))) %*% c(x1,x2)),rep(0,2), sigma=diag(2))}
x1 <- seq(-5, 1, by = 0.1)
x2 <- seq(-5, 2, by = 0.1)

pnormal(1,2)
pnormal1(1,2)
vnormal = Vectorize(pnormal)
vnormal1 = Vectorize(pnormal1)
vnormal(x1,x2)/vnormal1(x1,x2)

vnormal(1,2)/vnormal1(1,2)

draw2RV(vnormal, 0, 30)
draw2RV(vnormal1, 0, 30)
for(theta in seq(0,180,10)){
    draw2RV(vnormal, theta, 30)
}



}