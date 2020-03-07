x <- seq(0,10,.1)
y <- dgamma(x, 3,2)
plot(x,y)

n <- 200
#Sn <- sum(rgamma(n, 2,2))/n
Sn <- sapply(1:10000, function(x) sum(rgamma(n, 3,2))/n)
mean(Sn)
(3/2)
var(Sn)
3/4/n
hist(Sn)

hist(sqrt(n)*(Sn-(3/2)))
var(sqrt(n)*(Sn-(3/2)))
3/4

hist((Sn-(3/2)))
var((Sn-(3/2)))
3/4/n

#----


p <- seq(0,1,.01)
k <- 200
n <- 400
y <- k*log(p) + (n-k)*log(1-p)
#y <- p^k*(1-p)^(n-k)
Sp <- (k/p) - ((n-k)/(1-p))
Ip <- (n-k)/(1-p)^2+k/p^2
plot(p,y, type='l')
plot(p,Sp)
plot(p,Ip, type='l')

mle <- k/n
Ip.mle <- ((n-k)/(1-mle)^2)+(k/mle^2)
n/((1-mle)*mle)
n^3/(k^2-(k*n))*-1

#-- Exp of Score

#p <- seq(0,1,.01)
k <- seq(0,10,1)
n <- 10
#y <- k*log(p) + (n-k)*log(1-p)
#y <- p^k*(1-p)^(n-k)
fk <- dbinom(k, n, 0.7)
plot(k, fk)
Sp <- (k/0.7) - ((n-k)/(1-0.7))
plot(k, Sp)
sum(Sp*fk)

#-- Var of Score

#p <- seq(0,1,.01)
k <- seq(0,10,1)
n <- 10
#y <- k*log(p) + (n-k)*log(1-p)
#y <- p^k*(1-p)^(n-k)
fk <- dbinom(k, n, 0.7)
plot(k, fk)
Sp <- (k/0.7) - ((n-k)/(1-0.7))
plot(k, Sp)
sum((Sp)^2*fk)
n/(p*(1-0.7))

#-- Exp FI

#p <- seq(0,1,.01)
k <- seq(0,10,1)
n <- 10
#y <- k*log(p) + (n-k)*log(1-p)
#y <- p^k*(1-p)^(n-k)
fk <- dbinom(k, n, 0.7)
plot(k, fk)
Ip <- (k/0.7^2) + ((n-k)/(1-0.7)^2)
plot(k, Ip)
sum(Ip*fk)
n/(0.7*(1-0.7))



#----
opts <- list(
  algorithm = "NLOPT_LN_SBPLX", maxeval = "1000000",
  ftol_rel = .Machine$double.eps^0.5)

n <- 10 # basket size
prob <- 0.7
n.dat <- 100
n.sim <- 100
#dbinom(3, n, prob)
Xi <-sapply(1:n.sim, function(x) rbinom(n.dat, n, prob)) %>% t
#Xi <- rbinom(nsim, n, n.dat)

#dat <- Xi
#p <- 0.5
LikBin <- function(p, dat, basket.size){
  #n <- length(dat)
  out <- stats::dbinom(dat, size=basket.size, prob=p) %>% log %>% sum
  out*-1
}

#p <- seq(0,1, 0.1)
#y <- sapply(1:length(p), function(x)  LikBin(p[x], dat=Xi))
#plot(p, y)

#out <- nloptr(x0=0.5, eval_f = LikBin, lb = 0, ub = 1, opts = opts, dat=Xi)
out <-lapply(1:n.sim, function(x) nloptr(x0=0.5, eval_f = LikBin, lb = 0, ub = 1, opts = opts, dat=Xi[x,], basket.size=n) )

#p.hat <- out$solution
#logLik <- -out$objective
p.hat <- lapply(out, function(x) x$solution) %>% unlist
logLik <-lapply(out, function(x) -x$objective) %>% unlist

# x.space <- seq(0, n, 1)
# y.pred <- stats::dbinom(x.space, size=n, prob=p.hat)
# y.true <- stats::dbinom(x.space, size=n, prob=prob)
# plot(x.space, y.pred)
# points(x.space, y.true, col='red')

x.space <- seq(0, n, 1)
y.true <- stats::dbinom(x.space, size=n, prob=prob)
y.pred <-sapply(1:n.sim, function(x) stats::dbinom(x.space, size=n, prob=p.hat[x])) %>% t

Rn <- apply(y.true*log(y.pred), 1, sum)
Rn0 <- sum(y.true*log(y.true))
Rn0-Rn
KL <- y.true*log(y.true/y.pred)
KL <- apply(KL, 1, sum)

Rn[1]
(y.true*log(y.pred[1,])) %>% sum
sum(y.true)
hist(Rn, breaks=30)
K <- mean(Rn)
#K.hat <- logLik/n.dat
# suppose that 1st sample is the one we use for inference
logLik/n.dat
plot(Rn, logLik/n.dat)
plot(p.hat, logLik/n.dat)
hist(logLik/n.dat)
hist(Rn)
K.hat <- logLik[1]/n.dat

npar <- 1
aic <- (2*npar) - (2*logLik)
-aic/2*n
(K.hat-(npar/n))*100

mean(K.hat-K)
1/n

hist(logLik, breaks=30)


F.inv <- numDeriv::hessian(LikBin, x=prob, dat=Xi[1,], basket.size=n)
#F.inv * solve(F.inv)
solve(F.inv)*
numDeriv::grad(LikBin, x=prob, dat=Xi[1,], basket.size=n)^2


#sapply(1:length(Xi[1,]), function(x) stats::dbinom(Xi[1,x], size=n, prob=prob) %>% log)

#sapply(1:length(Xi[1,]), function(x) numDeriv::grad(dbinom, x=prob, dat=Xi[1,], basket.size=n))

#--- Normal distribution

LikNorm <- function(p, dat){
  out <- stats::dnorm(dat, mean = p[1], sd=p[2]) %>% log %>% sum
  out*-1
}

n.dat <- 10
dat <- rnorm(n.dat, 5, 2)

out <- nloptr(x0=c(1,1), eval_f = LikNorm, lb = c(-10, 0), ub =c(10, 10), opts = opts, dat=dat)

p.hat <- out$solution
logLik <- -out$objective

I.hat <- numDeriv::hessian(LikNorm, x=p.hat, dat=dat)
grad<- numDeriv::grad(LikNorm, x=p.hat, dat=dat)*-1
J.hat.com <- grad %*% t(grad)
d1 <- J.hat.com %*% solve(I.hat) %>% diag %>% sum
round(d1, 5)

#-- Correct
I.hat <- numDeriv::hessian(LikNorm, x=p.hat, dat=dat)
grad.sep <- lapply(1:n.dat, function(x) numDeriv::grad(LikNorm, x=p.hat, dat=dat[x])*-1)
J.hat.sep <- lapply(grad.sep, function(x) x %*% t(x))
J.hat.sep <- Reduce("+", J.hat.sep)
J.hat.sep %*% solve(I.hat) %>% diag %>% sum
