iterSearch_Sinba <- function(phy, data, edge.pattern, root.p, Pindv, ...){

  Qcom <- mk_Qcom(rate1=c(1,1), rate2=c(1,1))
  Qcom[Qcom==1] <- c(1:8)
  Qindv <- mk_Qindv(Qcom)

  out.list <- list()

  cat('Running iteration: 1')
  out <- runSinba(phy=phy, data=data, edge.pattern=edge.pattern,
                  root.p = root.p, Qindv=Qindv, Pindv=Pindv,
                  ...)

  out.list[[1]] <- out
  print(out)

  i <- 2
  max.par <-  max(Qindv$Qcom)

  while(max.par>1){
    Qcom <- iter_new_Qcom(pars=out$est.pars, Qindv$Qcom)
    Qindv <- mk_Qindv(Qcom)

    cat('Running iteration: ', i)

    out <- runSinba(phy=phy, data=data, edge.pattern=edge.pattern,
                    root.p = root.p, Qindv=Qindv, Pindv=Pindv,
                    ...)


    print(out)
    out.list[[i]] <- out
    i <- i+1
    max.par <-  max(Qindv$Qcom)

  }

  return(out.list)
}



#---


#' Title Find minimal edge pattern node for birth-death
#'
#' @param EdgeScenarios
#'
#' @return
#' @export
#'
#' @examples
birth_Des <-function(EdgeScenarios){
  event <- EdgeScenarios$edge.events
  x <- lapply(event, function(x) any(x=='*1*2' | x=='*2*1') ) %>% unlist

  event.sel <- event[x]
  ed <- lapply(event.sel, function(x) which(x=='*1*2' | x=='*2*1') ) %>% unlist
  nodes <- EdgeScenarios$tree$edge[ed,2]

  Des <- phangorn::Descendants(EdgeScenarios$tree,
                               node=nodes,
                               type="tips")

  lens <- lapply(Des, length) %>% unlist
  min <- which(lens==min(lens))
  event.sel[min]
}

#--


#   ____________________________________________________________________________
#   Linear Exp Simulations                                                  ####


#' Title Calculate analtical distribution of site patterns at two-branch tree
#'
#' @param Qsim rate matrix
#' @param t vector of times for each of 2 branches
#'
#' @return list of site patterns in binary format and probability for each pattern
#' @export
#'
#' @examples
getPhyloCTMC_distr <- function(Qsim, t){

  cmbs <- expand.grid(c1=c(1:4), c2=c(1:4)) %>% t

  # construct binary combination ofr each of two states
  v1 <- sapply(cmbs[1,], function(x) {z <- rep(0,4); z[x] <-1; z})
  v2 <- sapply(cmbs[2,], function(x) {z <- rep(0,4); z[x] <-1; z})
  #x <- 1

  # DECLARE p.root !!!!!!
  p.root <- rep(1/4,4)
  distr <- lapply(1:16, function(x) { sum(p.root * ((ExpMat_C(t[1], Qsim) %*% v1[,x]) * (ExpMat_C(t[2], Qsim) %*% v2[,x]))) } ) %>% unlist
  #sum(distr)
  #return(distr)
  return(list(char1=v1, char2=v2, distr=distr))
}

# kullback leibler divergence (full fomula)
KLD_phylo <- function(D.true, D.pred){
  sum( D.true * log(D.true/D.pred) )
}

# kullback leibler divergence (reduced formula)
KLD_K_phylo <- function(D.true, D.pred){
  #sum( D.true * log(D.true/D.pred) )
  sum(D.true*log(D.pred))
}

#' Title Simulate a two-branch tree with characters
#'
#' @param Q rate matrix
#' @param t vector of times for each of 2 branches
#' @param n.samples N of characters
#' @param pi intial vecotr
#'
#' @return list of two chars in binary mode (y- N states, x- N chars)
#' @export
#'
#' @examples
#t <- 1
simulate_ExpLinLn <- function(Q, t=c(1,1), n.samples=100, pi=c(1/4, 1/4, 1/4, 1/4)){

  # root statets
  ss.pi <- rmultinom(n.samples, size = 1, prob = pi)

  # Probability of seing  states at tips, x-states, y - N chars
  Pr1 <- t(ss.pi) %*% expm::expm((Q*t[1]), method = "Ward77")
  Pr2 <- t(ss.pi) %*% expm::expm((Q*t[2]), method = "Ward77")

  # get charatcers
  out1 <- apply(Pr1, 1, function(x) rmultinom(1, size = 1, prob = x))
  out2 <- apply(Pr2, 1, function(x) rmultinom(1, size = 1, prob = x))

  list(char1=out1, char2=out2)
}

# Likelihood function for two-branches tree
ExpLinLn <- function(p, out1, out2, Qfin, t, pi){

  if (any(is.nan(p)) || any(is.infinite(p))) {
    return(1e+06)
  }

  Q <- matrix(0,4,4)

  Q[] <- c(p, 0)[Qfin]
  diag(Q) <- -rowSums(Q)

  # Ln1 <-expm::expm((Q*t), method = "Ward77") %*% out1
  # Ln2 <-expm::expm((Q*t), method = "Ward77") %*% out2

  #EXP <- ExpMat_C(t, Q)
  EXP1 <- ExpMat_C(t[1], Q)
  EXP2 <- ExpMat_C(t[2], Q)

  Ln1 <-EXP1 %*% out1
  #EXP %*% out1[,2]
  Ln2 <-EXP2 %*% out2

  Ln <- -1*( (pi %*% (Ln1*Ln2)) %>% log() %>% sum )

  if (is.infinite(Ln)) {
    return(1e+06)
  }

  return(Ln)
}

# Likelihood function for two-branches tree to use in numerical integration
# for calculating marginal likelihood
ExpLinLn_Intg <- function(p, out1, out2, Qfin, t, pi, max.log){

  # if (any(is.nan(p)) || any(is.infinite(p))) {
  #   return(-1e+06)
  # }

  Q <- matrix(0,4,4)

  Q[] <- c(p, 0)[Qfin]
  diag(Q) <- -rowSums(Q)

  # Ln1 <-expm::expm((Q*t), method = "Ward77") %*% out1
  # Ln2 <-expm::expm((Q*t), method = "Ward77") %*% out2

  #EXP <- ExpMat_C(t, Q)
  EXP1 <- ExpMat_C(t[1], Q)
  EXP2 <- ExpMat_C(t[2], Q)

  #EXP <-expm::expm((Q*t), method = "Ward77")
  Ln1 <-EXP1 %*% out1
  Ln2 <-EXP2 %*% out2

  Ln <- ( (pi %*% (Ln1*Ln2)) %>% log() %>% sum )

  x <- exp(Ln-max.log)

  # if (is.infinite(Ln)) {
  #   return(-1e+06)
  # }

  return(x)
}


# Likelihood Inference for two-branched tree
# start.val starting value for MLn search
runExpLinLn <- function(Qcom, t, p.root=c(1/4, 1/4, 1/4, 1/4), data, start.val=1, add.TIC=FALSE){

  out1<-data$char1
  out2<-data$char2

  Qindv <- mk_Qindv(Qcom)
  Qrate <- mk_RatePars(Qindv)
  Qfin <- Qrate$Qcom

  opts <- list(
    algorithm = "NLOPT_LN_SBPLX", maxeval = "1000000",
    ftol_rel = .Machine$double.eps^0.5
  )


  np <- max(Qfin) -1
  ip <- rep(start.val, np)
  lower <- rep(0, np)
  upper <- rep(100, np)

  out <- nloptr(
    x0 = ip,
    eval_f = ExpLinLn, lb = lower, ub = upper,
    opts = opts,
    out1=out1, out2=out2, Qfin=Qfin, t=t, pi=p.root
  )


  loglik <- -out$objective
  est.pars <- out$solution
  solution <- matrix(est.pars[Qfin], dim(Qfin))

  #----- Empirical BIC
  #empirical I.hat
  h <- numDeriv::hessian(func = ExpLinLn, x = est.pars,
                         out1=out1, out2=out2, Qfin=Qfin, t=t, pi=p.root
  )

  BICemp <- loglik - (log(det(h))/2)


  # #- claculate Tic
  TIC.d = NA
  TIC= NA
  if (add.TIC){
    grad.sep <- lapply(1:ncol(out1), function(x)  numDeriv::grad(func = ExpLinLn, x = est.pars, out1=out1[,x], out2=out2[,x], Qfin=Qfin, t=t, pi=p.root) * -1 )
    J.hat.sep <- lapply(grad.sep, function(x) x %*% t(x))
    J.hat.sep <- Reduce("+", J.hat.sep)
    TIC.d <- J.hat.sep %*% solve(h) %>% diag %>% sum
    TIC <- -2 * loglik+(2*TIC.d)
  }

  # #-- Tic Correct for Normal
  # #I.hat <- numDeriv::hessian(LikNorm, x=p.hat, dat=dat)
  # # grad.sep <- lapply(1:n.dat, function(x) numDeriv::grad(LikNorm, x=p.hat, dat=dat[x])*-1)
  # # J.hat.sep <- lapply(grad.sep, function(x) x %*% t(x))
  # # J.hat.sep <- Reduce("+", J.hat.sep)
  # # J.hat.sep %*% solve(I.hat) %>% diag %>% sum
  # #----

  #--- claculate EDC
  # EDC1 <- -2*loglik + (np*log(log(2*ncol(out1))))
  # EDC2 <- -2*loglik + (2*np*log(log(2*ncol(out1))))


  #---------
  Qfull <- solution
  Qfull[is.na(Qfull)] <- 0
  diag(Qfull)<- -rowSums(Qfull)

  X <- list(
            loglik=loglik,
            n.par=np,
            AIC = -2 * loglik + 2 * np,
            #AICc = -2 * loglik + (2 * np * (ncol(out1) / (ncol(out1) - np - 1))),

            BIC= (np*log(2*ncol(out1)))-2*loglik,
            BIC.nsites= (np*log(ncol(out1)))-2*loglik,
            BIC2= (2*np*log(2*ncol(out1)))-2*loglik,
            #BIC.emp=BICemp,

            #CAIC= (np*log((2*ncol(out1)) + 1))-2*loglik,
            #Hessian=h,
            TIC.d = TIC.d,
            TIC= TIC,
            #EDC1=EDC1,
            #EDC2=EDC2,
            Q=solution,
            Qcom=Qcom,
            Qfull= Qfull,
            est.pars=est.pars)


  return(X)

}




# Iterative Model search for two-branched tree
iter_runExpLinLn <- function(t, p.root=c(1/4, 1/4, 1/4, 1/4), data, start.val=1){

  Qcom <- mk_Qcom(rate1=c(1,1), rate2=c(1,1))
  Qcom[Qcom==1] <- c(1:8)
  #Qindv <- mk_Qindv(Qcom)

  out.list <- list()

  cat('Running iteration: 1\n')
  out <-runExpLinLn(Qcom, t=t, p.root=p.root, data=data, start.val=start.val)
  # out <- runSinba(phy=phy, data=data, edge.pattern=edge.pattern,
  #                 root.p = root.p, Qindv=Qindv, Pindv=Pindv,
  #                 ...)

  out.list[[1]] <- out
  print(out)

  #i <- 2
  max.par <-  max(Qcom)

  while(max.par>1){
    Qcom <- iter_new_Qcom(pars=out$est.pars, Qcom)
    Qindv <- mk_Qindv(Qcom)

    cat('Running iteration: ', i, '\n')

    out <-runExpLinLn(Qcom, t=t, p.root=p.root, data=data, start.val=start.val)

    print(out)
    out.list[[i]] <- out
    i <- i+1
    max.par <-  max(Qindv$Qcom)

  }

  return(out.list)
}

# kl.sim <- mclapply(
#   Q.list, function(xx) infer2branchTrees(sim.multi, Qcom=xx, Qsim=Qsim, t=t, pi=p.root),
#   mc.cores = 4
# )
#
# sim.multi.a <- sim.multi
#sim.multi

#Qcom <- Q.list[[2]]

#' Title Run batch inference on two-branch tree
#'
#' @param sim.multi
#' @param Qsim
#' @param Qcom
#' @param t
#' @param pi
#'
#' @return
#' @export
#'
#' @examples
infer2branchTrees <- function(sim.multi, Qsim, Qcom, t, pi, add.TIC=FALSE){ # former names expected_KL

  n.sim <- length(sim.multi)

  # 2 run inference for 100
  r100 <- lapply(1:n.sim, function(xx)
    runExpLinLn(Qcom, t=t, p.root=pi, data=sim.multi[[xx]], start.val=1, add.TIC=add.TIC) )

  # 3 get pred distr
  #r100[[4]]
  D.true <- getPhyloCTMC_distr(Qsim, t=t)
  D.pred <- lapply(r100, function(x) getPhyloCTMC_distr(x$Qfull, t=t) )

  #4. clac KL smple and get mean
  KL <- sapply(1:n.sim, function(x) KLD_phylo(D.true$distr, D.pred[[x]]$distr) )
  KL_K <- sapply(1:n.sim, function(x) KLD_K_phylo(D.true$distr, D.pred[[x]]$distr) )

  # flatten parameters in Q to pass to output
  Q.flat <- lapply(r100, function(x) {zz <- c(x$Q); zz[!is.na(zz)] } )


  list(
    KL=KL,
    KL_part=KL_K,
    #KL.mean=mean(KL),
    #KL_K.mean=mean(KL_K),
    aic = lapply(r100, function(x) x$AIC ) %>% unlist,
    bic = lapply(r100, function(x) x$BIC ) %>% unlist,
    bic.nsite = lapply(r100, function(x) x$BIC.nsites ) %>% unlist,
    BIC2=lapply(r100, function(x) x$BIC2 ) %>% unlist,
    loglik = lapply(r100, function(x) x$loglik ) %>% unlist,
    tic.d= lapply(r100, function(x) x$TIC.d ) %>% unlist,
    tic= lapply(r100, function(x) x$TIC ) %>% unlist,
    #EDC1=lapply(r100, function(x) x$EDC1 ) %>% unlist,
    #EDC2=lapply(r100, function(x) x$EDC2 ) %>% unlist,
    #CAIC=lapply(r100, function(x) x$CAIC ) %>% unlist,
    #Q=lapply(r100, function(x) x$Qfull ),
    Q=Q.flat,
    n.par=Qcom[Qcom>0] %>% unique() %>% length() %>% rep(., length(KL))
  )

}

