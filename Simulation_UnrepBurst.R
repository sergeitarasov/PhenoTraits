
# simulate data
t <- 1
Qsim <- mk_Qcom(rate1=c(.1,.3), rate2=c(.1,.3))
diag(Qsim) <- -rowSums(Qsim)
Qsim
simL <- simulate_ExpLinLn(Q=Qsim, t=t, n.samples=1000, pi=c(1/4, 1/4, 1/4, 1/4))
#simulate_ExpLinLn(Q=Qsim, t=t, n.samples=2, pi=c(1/4, 1/4, 1/4, 1/4))

# infer
Qcom <- mk_Qcom(rate1=c(1,2), rate2=c(1,2))
p.root=c(1/4, 1/4, 1/4, 1/4)
r1 <- runExpLinLn(Qcom, t=t, p.root=c(1/4, 1/4, 1/4, 1/4), data=simL, start.val=1)


it <- iter_runExpLinLn(t=t, p.root=c(1/4, 1/4, 1/4, 1/4), data=simL, start.val=1)
r1$AIC
lapply(it, function(x) x$AIC)

r1$BIC
bic <- lapply(it, function(x) x$BIC) %>% unlist()
bic.em <- lapply(it, function(x) x$BICemp)  %>% unlist()
det(it[[1]]$Hessian)

plot(bic, bic.em)



#--- All partitions
library("parallel")
library(partitions)
library(dplyr)
Q.list <- mk_QcomPartitions8()

#-------------
# Q.list <- mk_QcomPartitions(8)
#
# xx=setparts(4)
# apply(xx, 2, function(x) length(unique(x)))
#
# npar=lapply(Q.list, function(x) unique(c(x)) %>% length ) %>% unlist
# hist(npar)
# sum(npar<=2)
#_-----------------

# run inference
infr <- mclapply(
  Q.list, function(x) runExpLinLn(x, t=t, p.root=c(1/4, 1/4, 1/4, 1/4), data=simL, start.val=1),
  mc.cores = 10
  )

r1$AIC
lapply(it, function(x) x$AIC)

aic <- lapply(infr, function(x) x$AIC) %>% unlist
min(aic)
which(aic<=r1$AIC)
(which(aic<=r1$AIC) %>% length())/length(aic)
(which(aic<=(r1$AIC-2) ) %>% length())/length(aic)
hist(aic, breaks = 30)
hist(aic-r1$AIC, breaks = 300, xlim =c(-10, 100))

delta <- abs(r1$AIC-min(aic))
(which(aic<=(r1$AIC-delta)) | which(aic>=(r1$AIC+delta) ) )
1-length(which(aic>(r1$AIC+delta) )) /length(aic)

r1$BIC
bic <- lapply(infr, function(x) x$BIC) %>% unlist
min(bic)
which(bic<=r1$BIC)
(which(bic<=r1$BIC) %>% length())/length(bic)
hist(bic-r1$BIC, breaks = 300, xlim =c(-10, 100))

1-length(which(bic>(r1$BIC+10) )) /length(bic)

bic.emp <- lapply(infr, function(x) x$BICemp) %>% unlist
lik <- lapply(infr, function(x) x$loglik) %>% unlist
plot(bic, bic.emp)
plot(bic, lik)
plot(bic, aic)
bic.emp[1:5]
bic[1:5]
exp(-2532.798)

r1$BICemp
max(bic.emp, na.rm = T)

plot(bic, -2*bic.emp)
dd <- (bic)-(-2*bic.emp)
max(dd,na.rm = T)
which(dd==max(dd,na.rm = T))
# [1] 2873
min(dd,na.rm = T)
bic-bic.emp
hist(dd, breaks=30)
dd.m <- which(dd> 2)
rbind((bic), -2*bic.emp)[,2873]

#-- Numerical Integration
library(cubature)        # load the package "cubature"
f <- function(x) { 2/3 * (x[1] + x[2] + x[3]) }  # "x" is vector
adaptIntegrate(f, lowerLimit = c(0, 0, 0), upperLimit = c(0.5, 0.5, 0.5))


# infer
which(dd>10)
rbind((bic), -2*bic.emp)[,133]
Qcom <-Q.list[[133]]
Qcom <- mk_Qcom(rate1=c(1,2), rate2=c(1,2))
p.root=c(1/4, 1/4, 1/4, 1/4)
r1 <- runExpLinLn(Qcom, t=t, p.root=c(1/4, 1/4, 1/4, 1/4), data=simL, start.val=1)
#---
#p <- c(.1)
#--

out1<-data$char1
out2<-data$char2

Qindv <- mk_Qindv(Qcom)
Qrate <- mk_RatePars(Qindv)
Qfin <- Qrate$Qcom

max.log <- r1$loglik

#ExpLinLn_Intg(p=p, out1=out1, out2=out2, Qfin=Qfin, t=t, pi=p.root)
#ExpLinLn(p=p, out1=out1, out2=out2, Qfin=Qfin, t=t, pi=p.root)

int <- adaptIntegrate(ExpLinLn_Intg, lowerLimit = rep(0,3), upperLimit = rep(Inf, 3), tol = 1e-05, fDim = 1, vectorInterface = T,
               out1=out1, out2=out2, Qfin=Qfin, t=t, pi=p.root, max.log=max.log)

int <- cubintegrate(f = ExpLinLn_Intg, lower = rep(0,3), upper = rep(Inf, 3), method = "hcubature", fDim = 1,
                    out1=out1, out2=out2, Qfin=Qfin, t=t, pi=p.root, max.log=max.log)

# int <- pcubature(ExpLinLn_Intg, lowerLimit = rep(0,3), upperLimit = rep(Inf, 3), fDim = 1,
#                        out1=out1, out2=out2, Qfin=Qfin, t=t, pi=p.root, max.log=max.log)

(log(int$integral)+max.log)*-2
r1$BICemp*-2
r1$BIC

#---------------------------
#------- AIC
combn(c(1:4), 2)
expand.grid(c1=c(1:4), c2=c(1:4))
t <- 1
p.root=c(1/4, 1/4, 1/4, 1/4)
#-- true distr
Qsim <- mk_Qcom(rate1=c(.1,.3), rate2=c(.1,.3))
diag(Qsim) <- -rowSums(Qsim)
Qsim
#simL <- simulate_ExpLinLn(Q=Qsim, t=t, n.samples=1000, pi=c(1/4, 1/4, 1/4, 1/4))

getPhyloCTMC_distr <- function(Qsim, t){

  cmbs <- expand.grid(c1=c(1:4), c2=c(1:4)) %>% t

  # construct binary combination ofr each of two states
  v1 <- sapply(cmbs[1,], function(x) {z <- rep(0,4); z[x] <-1; z})
  v2 <- sapply(cmbs[2,], function(x) {z <- rep(0,4); z[x] <-1; z})
  #x <- 1
  distr <- lapply(1:16, function(x) { sum(p.root * ((ExpMat_C(t[1], Qsim) %*% v1[,x]) * (ExpMat_C(t[2], Qsim) %*% v2[,x]))) } ) %>% unlist
  #sum(distr)
  #return(distr)
  return(list(char1=v1, char2=v2, distr=distr))
}

KLD_phylo <- function(D.true, D.pred){
  sum( D.true * log(D.true/D.pred) )
}

KLD_K_phylo <- function(D.true, D.pred){
  #sum( D.true * log(D.true/D.pred) )
  sum(D.true*log(D.pred))
}

# simulation

1. simuate 100
2. run inference for 100
3. get pred distr
4. clac KL smple and get mean


Qcom <- mk_Qcom(rate1=c(1,1), rate2=c(1,1))

Qsim <- mk_Qcom(rate1=c(.1,.1), rate2=c(.1,.1))
#Qsim <- mk_Qcom(rate1=c(.1,.2), rate2=c(.3,.4))
diag(Qsim) <- -rowSums(Qsim)
Qsim

p.root=c(1/4, 1/4, 1/4, 1/4)
n.sim <- 500
t <- 1

# 1 simulate
sim.multi <- lapply(1:n.sim, function(xx) simulate_ExpLinLn(Q=Qsim, t=t, n.samples=100, pi=c(1/4, 1/4, 1/4, 1/4)) )
length(sim.multi)
#expected_KL(sim.multi, Qcom, t=1, pi=p.root)

# many models
Q.list <- mk_QcomPartitions()
lapply(Q.list, function(x) x==Qcom) # 8
Q.list[[8]]
Q.list[[11]]

kl.sim <- mclapply(
  Q.list, function(xx) expected_KL(sim.multi, Qcom=xx, Qsim=Qsim, t=1, pi=p.root),
  mc.cores = 10
)

##-----
#sim.multi[[1]]$char1 %>% dim
dim.dat <- 1000*2
ln<- sapply(kl.sim, function(x) x$loglik) %>% t
K.mean <- sapply(kl.sim, function(x) x$KL_K.mean)
K <- sapply(kl.sim, function(x) x$K) %>% t

best.K <- apply(K, 2, function(x) which(x==max(x))) %>% unlist
delta.K <- apply(K, 2, function(x) x[1]-max(x) )*-1
(delta.K[delta.K>3]) %>% length
max(delta.K)



npar <- sapply(kl.sim, function(x) x$n.par)
dim.dat <- 1
1/dim.dat
K.hat <- ln/dim.dat
mean(K[1,]-K.mean[1])
Exp.K <- sapply(1:nrow(K),  function(x) mean(K[x,]-K.mean[x]) )
plot(Exp.K, (npar/dim.dat) )


aic <- sapply(kl.sim, function(x) x$aic) %>% t
apply(aic, 1, mean)
best.aic <- apply(aic, 2, function(x) which(x==min(x)))
delta.aic <- apply(aic, 2, function(x) x[1]-min(x) )
(delta.aic[delta.aic>2]) %>% length
max(delta.aic)
which(delta.aic==max(delta.aic))
zz <- apply(aic, 2, function(x) which(x[1]-x > 2) %>% length )
zz[zz>0] %>% sum

# tic <- sapply(kl.sim, function(x) x$tic.d) %>% t
# round(tic, 4)
# round(tic, 4)[,1]
tic <- sapply(kl.sim, function(x) x$tic) %>% t
apply(tic, 1, mean)
best.tic <- apply(tic, 2, function(x) which(x==min(x)))
delta.tic <- apply(tic, 2, function(x) x[1]-min(x) )
(delta.tic[delta.tic>3]) %>% length
max(delta.tic)

tic.d <- sapply(kl.sim, function(x) x$tic.d) %>% t
cbind(tic.d, npar)
round(tic.d-npar, 2) %>% abs %>% max

aic.pool <- apply(aic, 2, function(x) x[1]-x ) %>% c()
max(aic.pool)
which(aic.pool> -2) %>% length()
hist(aic.pool, breaks = 40)

kl <- sapply(kl.sim, function(x) x$KL) %>% t
hist(kl[2,])
log(kl)
best.kl <- apply(kl, 2, function(x) which(x==min(x)))
best.kl[best.kl!=1] %>% length

apply(kl, 1, mean)
hist(kl[15,])
plot(aic[,60], kl[,60])
rbind(best.kl, best.aic)

kl.m <- sapply(kl.sim, function(x) x$KL.mean)
which(kl.m ==min(kl.m))
exp(kl.m)
plot(aic[,71], kl.m)
plot(K.mean, kl.m)

bic <- sapply(kl.sim, function(x) x$bic) %>% t
plot(aic[,1], bic[,1])
plot(bic[,8], kl.m)
plot(bic[,2], kl[,2])

best.bic <- apply(bic, 2, function(x) which(x==min(x)))
delta.bic <- apply(bic, 2, function(x) x[1]-min(x) )
which(delta.bic>3) %>% length()
max(delta.bic)

bic.n <- sapply(kl.sim, function(x) x$bic.nsite) %>% t
best.bic.n <- apply(bic.n, 2, function(x) which(x==min(x)))
delta.bic.n <- apply(bic.n, 2, function(x) x[1]-min(x) )
which(delta.bic.n>3) %>% length()
max(delta.bic.n)

edc1 <- sapply(kl.sim, function(x) x$EDC2) %>% t
best.edc1 <- apply(edc1, 2, function(x) which(x==min(x)))
delta.edc1 <- apply(edc1, 2, function(x) x[1]-min(x) )

bic2 <- sapply(kl.sim, function(x) x$BIC2) %>% t
best.bic2 <- apply(bic2, 2, function(x) which(x==min(x)))
delta.bic2 <- apply(bic2, 2, function(x) x[1]-min(x) )
which(delta.bic2>3) %>% length()
max(delta.bic2)

caic <- sapply(kl.sim, function(x) x$CAIC) %>% t
best.caic <- apply(caic, 2, function(x) which(x==min(x)))
delta.caic <- apply(caic, 2, function(x) x[1]-min(x) )
which(delta.caic>3) %>% length()
max(delta.caic)


dim.dat <- 100*2
spply(ln)
-2*ln+(2*np*log(log()))
edc1 <- aic
2*(log(log(1000)))
(log(log(1000)))

bic.pool <- apply(bic, 2, function(x) x[1]-x ) %>% c()
max(bic.pool)
which(bic.pool> -2) %>% length()
which(bic.pool> 2)
hist(bic.pool, breaks = 40)

bic2 <- sapply(kl.sim, function(x) x$bic.nsite) %>% t
best.bic2 <- apply(bic2, 2, function(x) which(x==min(x)))
delta.bic2 <- apply(bic2, 2, function(x) x[1]-min(x) )
max(delta.bic2)

#-- Ln ration
npar
daic <- aic[1,]-aic[15,]
which(daic>2) %>% length(.)/500
hist(daic, breaks=30, freq = F)

apply(aic, 1, function(x) {daic <- aic[1,]-x; which(daic>2) %>% length(.)})
delta.aic

ra <- (ln[1,]-ln[15,])*-2
all(ra>0)
hist(ra, breaks=30, freq = F)

x <- seq(0, 15, 0.1)
y <- dchisq(x, df=3)
lines(x, y)

1-pchisq(8, df=3)

ra0 <- ln[15,]-ln[1,]
which(ra0>4) %>% length(.)/500
which(ra>8) %>% length(.)/500
hist(ra0, breaks = 30)

(-2*ln[15,]+8) - aic[15,]
#---

out1<-sim.multi[[2]]$char1
out2<-sim.multi[[2]]$char2

Q.list[[15]]
Qindv <- mk_Qindv(Q.list[[9]])
Qrate <- mk_RatePars(Qindv)
Qfin <- Qrate$Qcom

max.log <- kl.sim[[9]]$loglik[[2]]
kl.sim[[9]]$bic[[2]]
kl.sim[[15]]$bic[[2]]

int <- adaptIntegrate(ExpLinLn_Intg, lowerLimit = rep(0,3), upperLimit = rep(Inf, 3), tol = 1e-05, fDim = 1, vectorInterface = F,
                      out1=out1, out2=out2, Qfin=Qfin, t=t, pi=p.root, max.log=max.log)


(log(int$integral)+max.log)*-2
kl.sim[[9]]$bic[[2]]
kl.sim[[9]]$bic.nsite[[2]]
#r1$BICemp*-2
#r1$BIC

# infr <- mclapply(
#   Q.list, function(x) runExpLinLn(x, t=t, p.root=c(1/4, 1/4, 1/4, 1/4), data=simL, start.val=1),
#   mc.cores = 10
# )
#expected_KL(sim.multi, Qcom=xx, Qsim=Qsim, t=c(1, 0.5), pi=p.root)
expected_KL <- function(sim.multi, Qsim, Qcom, t, pi){

  n.sim <- length(sim.multi)

  # 2 run inference for 100
  r100 <- lapply(1:n.sim, function(xx)
    runExpLinLn(Qcom, t=t, p.root=pi, data=sim.multi[[xx]], start.val=1) )

  # 3 get pred distr
  D.true <- getPhyloCTMC_distr(Qsim, t=t)
  D.pred <- lapply(r100, function(x) getPhyloCTMC_distr(x$Qfull, t=t) )

  #4. clac KL smple and get mean
  KL <- sapply(1:n.sim, function(x) KLD_phylo(D.true$distr, D.pred[[x]]$distr) )
  KL_K <- sapply(1:n.sim, function(x) KLD_K_phylo(D.true$distr, D.pred[[x]]$distr) )


  list(
  KL=KL,
  KL_part=KL_K,
  #KL.mean=mean(KL),
  #KL_K.mean=mean(KL_K),
  aic = lapply(r100, function(x) x$AIC ) %>% unlist,
  bic = lapply(r100, function(x) x$BIC ) %>% unlist,
  bic.nsite = lapply(r100, function(x) x$BIC.nsites ) %>% unlist,
  loglik = lapply(r100, function(x) x$loglik ) %>% unlist,
  #tic.d=NA, #lapply(r100, function(x) x$TIC.d ) %>% unlist,
  #tic=NA, #lapply(r100, function(x) x$TIC ) %>% unlist,
  EDC1=lapply(r100, function(x) x$EDC1 ) %>% unlist,
  EDC2=lapply(r100, function(x) x$EDC2 ) %>% unlist,
  BIC2=lapply(r100, function(x) x$BIC2 ) %>% unlist,
  CAIC=lapply(r100, function(x) x$CAIC ) %>% unlist,
  Q=lapply(r100, function(x) x$Qfull ),
  n.par=Qcom[Qcom>0] %>% unique() %>% length() %>% rep(., length(KL))
  )

}




hist(ln)
hist(KL)
hist(aic)
plot(aic, KL)
plot(bic, KL)
plot(aic, bic)


#--
# Qcom <- mk_Qcom(rate1=c(1,2), rate2=c(1,2))
# p.root=c(1/4, 1/4, 1/4, 1/4)
# r1 <- runExpLinLn(Qcom, t=t, p.root=c(1/4, 1/4, 1/4, 1/4), data=simL, start.val=1)

# KL using full distr
D.true <- getPhyloCTMC_distr(Qsim)
D.pred <- getPhyloCTMC_distr(r1$Qfull)
KLD_phylo(D.true, D.pred)

# KL using only Qs
D.true1 <- c(ExpMat_C(t=1, Qsim))
D.pred1 <-c(ExpMat_C(t=1, r1$Qfull))
KLD_phylo(D.true1, D.pred1)

#---------------------------------
# Poissom
dt <- rpois(100, 2.5)
hist(dt)

lam=1
PoisLn <- function(lam, dt){
  prod(dpois(dt, lam))
}

PoisLn_log <- function(lam, dt){
  sum(dpois(dt, lam,  log = T))
}

lam <- seq(0, 10, 0.1)
#ln <- PoisLn(lam=1, dt)
ln <-sapply(lam, function(x) PoisLn(x, dt))
plot(lam, ln, type='l')

ln.int <- adaptIntegrate(PoisLn, lowerLimit = c(0), upperLimit = c(Inf), dt=dt)
log(ln.int$integral)

ln.log <-sapply(lam, function(x) PoisLn_log(x, dt))
plot(lam, ln.log, type='l')

max.log <- max(ln.log)+200
exp(max.log)

#lam <- 1
PoisLn_Int <- function(lam, dt, max.log){
 x <- sum(dpois(dt, lam, log = T))-max.log
 exp(x)
}

ln.int <-sapply(lam, function(x) PoisLn_Int(x, dt, max.log))
plot(lam, ln.int, type='l')

intr1 <- adaptIntegrate(PoisLn_Int, lowerLimit = c(0), upperLimit = c(Inf), dt=dt, max.log=max.log)
log(intr1$integral)+max.log




#-- AIC Pois

lam.true <- 2.5
N <- 100
k <- 1
dt <- lapply(1:1, function(x) rpois(N, lam.true))

PoisLn <- function(lam, dt){
  prod(dpois(dt, lam))
}

PoisLn_log <- function(lam, dt){
  sum(dpois(dt, lam,  log = T))
}

lam <- seq(0, 10, 0.1)
ln <-sapply(lam, function(x) PoisLn(x, dt[[1]]))
plot(lam, ln, type='l')

# Max Lik
Mln.par <-lapply(dt, function(x) sum(x)/N) %>% unlist
hist(Mln.par, breaks = 20)

ln.log <-sapply(lam, function(x) PoisLn_log(x, dt[[1]]))
plot(lam, ln.log, type='l')
#Mln.val <- lapply(1:N, function(x)  PoisLn_log(Mln.par[x], dt[[x]]))%>% unlist
#aic <- 2*k-(2*Mln.val)
aic <- 2*k-(2*ln.log)
plot(ln.log, aic)
hist(aic)

# Ped lik
theta <- c(1:100)
pred.distr <-  lapply(1:length(lam), function(x)  dpois(theta, lam[x]) )
plot(1:N, pred.distr[[2]])

true.dist <- dpois(theta, lam.true)
plot(1:N, true.dist)

KLD <- function(true.dist, pred.distr){
  n <- length(pred.distr)
  #x <- 1
  lapply(1:n, function(x) sum( true.dist*log(true.dist/pred.distr[[x]]) )  ) %>% unlist
}

kld <- KLD(true.dist, pred.distr)
plot(lam, kld)
hist(kld)
plot(aic, kld, type='l')
plot(ln.log, kld, type='l')

#   ____________________________________________________________________________
#   Simulate Unreplicated Burst                                             ####


token.maps <- c("00", "01", "10", "11")
names(token.maps) <-  c(1:4)
Q <- mk_Qcom(rate1=c(.1,.1), rate2=c(.1,.1))
diag(Q) <- -rowSums(Q)
Q

sim.UB <- simUnrepBurst(ntips=100, n.trees.sim=1, Upper.otu.bound=0.45, Lower.otu.bound=0.4, Q, p.root=c(1,0,0,0), freq.threshhold=5)

# Plot
i <- 1
tree <- sim.UB[[i]]$tree
chars <- cbind(sim.UB[[i]]$states.binary[[1]], sim.UB[[i]]$states.binary[[2]], sim.UB[[i]]$states.original.recoded)
plot.phylo(tree, edge.width = 3, show.tip.label = T, no.margin = T, label.offset=10, tip.color='white' )
addTipChars(chars,  r=.4,  x.space=3, legend=T)

# Edge Scenarios
tree <- sim.UB[[i]]$tree
chars <- cbind(sim.UB[[i]]$states.binary[[1]], sim.UB[[i]]$states.binary[[2]])
#chars <- cbind(names(sim.UB[[i]]$states.original), sim.UB[[i]]$states.original)
scen <- getEdgeScenarios(chars, focal.states=c(1,1), tree, token.maps)
scen$edge.events

# Plot
EE <- birth_Des(scen)[[1]]
#EE <- scen$edge.events[[1]]
edge.color <-edgePatRecode(EE, from='codes', to='cols')
dat <- scen$data[,2:4]
rownames(dat) <- scen$data[,1]
plot.phylo(scen$tree, edge.color = edge.color, edge.width = 3, show.tip.label = T, no.margin = T, label.offset=10, tip.color='white' )
addTipChars(x=dat, colors, r=.4,  x.space=3, legend=TRUE)
edgelabels(EE, cex = 0.4, frame = "ci", bg = edge.color)


# Make Qs
Qcom <- mk_Qcom(rate1=c(1,1), rate2=c(1,1), names=names(token.maps))

Qindv <- mk_Qindv(Qcom)
Pindv <- mk_Pindv(Qcom, p1=c(0.5,0,0.5,0), p2=c(0.5,.5,0,0), p0=c(1,0,0,0))
#Pindv <- mk_Pindv(Qcom, p1=c(0,0,1,0), p2=c(0,1,0,0), p0=c(1,0,0,0))
root.p <- c(1, 0,0,0)
#root.p <- 'maddfitz'
#root.p <- rep(1/4, 4)

out1 <- runSinba(phy=scen$tree, data=scen$data, edge.pattern=EE,
                root.p = root.p, Qindv=Qindv, Pindv=Pindv,
                p = NULL, ip = NULL, lb = 0, ub = 100,
                cub.method='pcubature', diagn=F, optim.method='subplex', subplex.rounds=1, subplex.n.cores=1L, gensa.its=100)


itre <- iterSearch_Sinba(phy=scen$tree, data=scen$data, edge.pattern=EE, root.p = root.p, Pindv=Pindv,
                         p = NULL, ip = NULL, lb = 0, ub = 100,
                         cub.method='pcubature', diagn=F, optim.method='subplex', subplex.rounds=1, subplex.n.cores=1L, gensa.its=100)

char.or <- cbind(names(sim.UB[[i]]$states.original), sim.UB[[i]]$states.original)
out1a <- runSinba(phy=scen$tree, data=char.or, edge.pattern=NULL,
                 root.p = root.p, Qindv=Qindv, Pindv=Pindv,
                 p = NULL, ip = NULL, lb = 0, ub = 100,
                 cub.method='pcubature', diagn=F, optim.method='subplex', subplex.rounds=1, subplex.n.cores=1L, gensa.its=100)


Qcom <- mk_Qcom(rate1=c(1,2), rate2=c(3,4), names=names(token.maps))
Qindv <- mk_Qindv(Qcom)
Pindv <- mk_Pindv(Qcom, p1=c(0,0,1,0), p2=c(0,1,0,0), p0=c(1,0,0,0))
root.p <- c(1, 0,0,0)

out2 <- runSinba(phy=scen$tree, data=scen$data, edge.pattern=NULL,
                 root.p = root.p, Qindv=Qindv, Pindv=Pindv,
                 p = NULL, ip = NULL, lb = 0, ub = 100,
                 cub.method='pcubature', verbose = TRUE, diagn=F, optim.method='subplex', subplex.rounds=1, subplex.n.cores=1L, gensa.its=100)

#out2$solution.subplex.list
#-----
Qcom.cor <- mk_Qcom(rate1=c(1,1), rate2=c(1,1), names=names(token.maps))
Qcom.cor[Qcom.cor==1] <- c(1:8)

Qindv <- mk_Qindv(Qcom.cor)
Pindv <- mk_Pindv(Qcom.cor, p1=c(0,0,1,0), p2=c(0,1,0,0), p0=c(1,0,0,0))
root.p <- c(1, 0,0,0)

out3 <- runSinba(phy=scen$tree, data=scen$data, edge.pattern=NULL,
                 root.p = root.p, Qindv=Qindv, Pindv=Pindv,
                 p = NULL, ip = NULL, lb = 0, ub = 100,
                 cub.method='pcubature', diagn=F, optim.method='subplex', subplex.rounds=1, subplex.n.cores=1L, gensa.its=100)

out3
Qcom.cor <- iter_new_Qcom(pars=out3$est.pars, Qcom.cor)
#--------------
Qcom <- mk_Qcom(rate1=c(1,1), rate2=c(1,1), names=names(token.maps))

Qindv <- mk_Qindv(Qcom)
Pindv <- mk_Pindv(Qcom, p1=c(0,0,1,0), p2=c(0,1,0,0), p0=c(1,0,0,0))
root.p <- c(1, 0,0,0)
root.p <- 'maddfitz'

#edge.pattern=scen$edge.patterns[[5]]

out <- runSinba(phy=scen$tree, data=chars, edge.pattern=NULL,
                 root.p = root.p, Qindv=Qindv, Pindv=Pindv,
                 p = NULL, ip = NULL, lb = 0, ub = 100,
                 cub.method='pcubature', diagn=F, optim.method='subplex', subplex.rounds=1, subplex.n.cores=1L, gensa.its=100)
#--
Qcom <- mk_Qcom(rate1=c(1,1), rate2=c(1,1), names=names(token.maps))
Qcom[Qcom==1] <- c(1:8)
Qindv <- mk_Qindv(Qcom)
Pindv <- mk_Pindv(Qcom, p1=c(0.5,0,0.5,0), p2=c(0.5,0.5,0,0), p0=c(1,0,0,0))
root.p <- c(1, 0,0,0)
root.p <- rep(1/4,4)

iter <- iterSearch_Sinba(phy=scen$tree, data=scen$data, edge.pattern=NULL,
                 root.p = root.p, Qindv=Qindv, Pindv=Pindv,
                 p = NULL, ip = NULL, lb = 0, ub = 100,
                 cub.method='pcubature', diagn=F, optim.method='subplex', subplex.rounds=1, subplex.n.cores=1L, gensa.its=100)

out <- runSinba(phy=scen$tree, data=scen$data, edge.pattern=scen$edge.patterns[[5]],
                 root.p = root.p, Qindv=Qindv, Pindv=Pindv,
                 p = NULL, ip = NULL, lb = 0, ub = 100,
                 cub.method='pcubature', diagn=F, optim.method='subplex', subplex.rounds=1, subplex.n.cores=1L, gensa.its=100)

#---------------------
extract.clade

Qcom <- mk_Qcom(rate1=c(1,1), rate2=c(1,1), names=names(token.maps))
Qcom[Qcom==1] <- c(1:8)
Qindv <- mk_Qindv(Qcom)

akaike.weights(c(out1a$AIC, out.list[[5]]$AIC))
ll <- lapply(out.list, function(x) x$AIC) %>% unlist
akaike.weights(c(out1a$AIC, ll))



#-----

Qcom.cor <- mk_Qcom(rate1=c(1,1), rate2=c(1,1), names=names(token.maps))
Qcom.cor[1,2] <- 2
Qcom.cor[1,3] <- 2

Qcom.cor[3,1] <- 3
Qcom.cor[2,1] <- 3

Qindv <- mk_Qindv(Qcom.cor)
Pindv <- mk_Pindv(Qcom.cor, p1=c(0,0,1,0), p2=c(0,1,0,0), p0=c(1,0,0,0))
root.p <- c(1, 0,0,0)

out4 <- runSinba(phy=scen$tree, data=scen$data, edge.pattern=NULL,
                 root.p = root.p, Qindv=Qindv, Pindv=Pindv,
                 p = NULL, ip = NULL, lb = 0, ub = 100,
                 cub.method='pcubature', verbose = TRUE, diagn=F, optim.method='subplex', subplex.rounds=1, subplex.n.cores=1L, gensa.its=300)




Qcom.list <- mk_QcomPartitions(num.elenents=4)
scen$edge.events


library(gplots)
y <- matrix(rnorm(500), 100, 5, dimnames=list(paste("g", 1:100, sep=""), paste("t", 1:5, sep="")))
heatmap.2(y) #

hc <- hclust(dist(USArrests), "ave")
plot(hc)
plot(hc, hang = -1)

x <- c(1:5)
hclust(as.dist(x), method = "complete", members = NULL)
## Row- and column-wise clustering
hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete")
## Tree cutting
mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)]
## Plot heatmap
mycol <- colorpanel(40, "darkblue", "yellow", "white") # or try redgreen(75)
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc)
