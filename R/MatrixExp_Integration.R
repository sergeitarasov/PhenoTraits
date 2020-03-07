#   ____________________________________________________________________________
#   Matrix Exponentials and Integration                                     ####




# ExpMat_C(double t, arma::mat A)
exp_Q <- function(t, Q) {
  out <- expm::expm(Q * t, method = "Ward77")
  return(c(out))
}


# Exp_2Q_C(double t1, double t2, double Tau, arma::mat Q1, arma::mat Q2)
# Exp(Q1*t2-t1)*Exp(Q2*Tau-t2)
exp_2Q <- function(t1, t2, Tau, Q1, Q2) {
  out <- expm::expm(Q1 * (t2 - t1), method = "Ward77") %*% expm::expm(Q2 * (Tau - t2), method = "Ward77")
  return(c(out))
}


# Exp_2Q_scen3_C(double t1, double Tau, arma::mat Q1, arma::mat Q2)
# Exp(Q1*t1)*Exp(Q2*Tau-t1)
exp_2Q_scen3 <- function(t1, Tau, Q1, Q2) {
  out1 <- expm::expm(Q1 * t1, method = "Ward77")
  out2 <- expm::expm(Q2 * (Tau - t1), method = "Ward77")
  out <- out1 %*% out2
  return(c(out))
}

#---- Vectorization

exp_Q_vec <- Vectorize(exp_Q, vectorize.args = c("t"), SIMPLIFY = T)
exp_2Q_vec <- Vectorize(exp_2Q, vectorize.args = c("t1", "t2"), SIMPLIFY = T)
exp_2Q_scen3_vec <- Vectorize(exp_2Q_scen3, vectorize.args = c("t1"), SIMPLIFY = T)



#--------------------------- Fisrt Integral

#--------------
#--- R code
#--------------

# Int[0, t2]dt1: Exp(Q1*t2-t1)*Exp(Q2*Tau-t2)
integrate_Exp <- function(t2, Q1, Q2, Tau, f.Dim = 4, cub.method = "pcubature") {
  int <- cubintegrate(
    f = exp_2Q_vec, lower = 0, upper = t2, method = cub.method, fDim = f.Dim,
    t2 = t2, Tau = Tau, Q1 = Q1, Q2 = Q2
  )

  return(int$integral)
}


# Integral for *1 or *2; Scen #1 and #2
# 1/2*integrateOneBirth(Tau=2, Q$Q1)
integrateOneBirth <- function(Tau, Q1, cub.method = "pcubature") {
  f.Dim <- ncol(Q1)^2
  out <- cubintegrate(f = exp_Q_vec, lower = 0, upper = Tau, method = cub.method, fDim = f.Dim, Q = Q1)

  out1 <- matrix(out$integral, ncol = ncol(Q1))
  colnames(out1) <- rownames(out1) <- rownames(Q1)
  return(out1)
}


# Integral for 1*2; Scen#3
# 1/2*integrateScen3(Tau=2, Q1, Q1)
integrateScen3 <- function(Tau, Q1, Q2, cub.method = "pcubature") {
  f.Dim <- ncol(Q1)^2

  out <- cubintegrate(f = exp_2Q_scen3_vec, lower = 0, upper = Tau, method = cub.method, fDim = f.Dim, Tau = Tau, Q = Q1, Q2 = Q2)

  out1 <- matrix(out$integral, ncol = ncol(Q1))
  colnames(out1) <- rownames(out1) <- rownames(Q1)
  return(out1)
}

#--- R code Vecorized for the Second Integral

# first Integral Vec
integrate_Exp_vec <- Vectorize(integrate_Exp, vectorize.args = "t2", SIMPLIFY = T)

#--------------
#--- C++ code
#--------------

integrate_Exp_C <- function(t2, Q1, Q2, Tau, f.Dim = 4) {
  int <- pcubature(f = vExp_2Q_t1_C, lowerLimit = 0, upperLimit = t2, fDim = f.Dim, t2 = t2, Tau = Tau, Q1 = Q1, Q2 = Q2, tol = 1e-05, vectorInterface = FALSE)
  # int <-cubintegrate(f = vExp_2Q_t1_C, lower = 0, upper = t2, method = cub.method, fDim=f.Dim,
  #                    t2=t2, Tau=Tau, Q1=Q1, Q2=Q2)
  return(int$integral)
}

# Integral for *1 or *2; Scen #1 and #2
integrateOneBirth_C <- function(Tau, Q1, cub.method = "pcubature") {
  f.Dim <- ncol(Q1)^2

  out <- cubintegrate(f = vExpMat_C, lower = 0, upper = Tau, method = cub.method, fDim = f.Dim, Q = Q1)

  out1 <- matrix(out$integral, ncol = ncol(Q1))
  colnames(out1) <- rownames(out1) <- rownames(Q1)
  return(out1)
}

# Integral for 1*2; Scen#3
integrateScen3_C <- function(Tau, Q1, Q2, cub.method = "pcubature") {
  f.Dim <- ncol(Q1)^2

  out <- cubintegrate(f = vExp_2Q_scen3_C, lower = 0, upper = Tau, method = cub.method, fDim = f.Dim, Tau = Tau, Q = Q1, Q2 = Q2)

  out1 <- matrix(out$integral, ncol = ncol(Q1))
  colnames(out1) <- rownames(out1) <- rownames(Q1)
  return(out1)
}



#--- C++ code Vecorized for the Second Integral

integrate_Exp_C_vec <- function(t2, Q1, Q2, Tau, f.Dim = 4) {
  int <- sapply(t2, function(x) integrate_Exp_C(x, Q1, Q2, Tau, f.Dim = f.Dim))
  return(int)
}


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### MicroBenchmark                                                          ####

# microbenchmark(
#   exp_Q_vec(c(1,2), Q$Qcom),
#   vExpMat_C(c(1,2), Q$Qcom)
# )
#
# microbenchmark(
# integrate_Exp_vec(t2=2, Q1=Q$Qcom, Q2=Q$Q1, Tau=5, f.Dim=4, cub.method='pcubature'),
# integrate_Exp_C_vec(t2=2, Q1=Q$Qcom, Q2=Q$Q1, Tau=5, f.Dim=4)
# )
#
# microbenchmark(
# integrateOneBirth(Tau=1, Q1=Q$Qcom, cub.method='pcubature'),
# integrateOneBirth_C(Tau=1, Q1=Q$Qcom, cub.method='pcubature')
# )
#
# microbenchmark(
# integrateScen3(Tau=1, Q1=Q$Qcom, Q2=Q$Q1, cub.method='pcubature'),
# integrateScen3_C(Tau=1, Q1=Q$Qcom, Q2=Q$Q1, cub.method='pcubature')
# )


#--------------------------- Second Integral

integrate2_Exp <- function(Q1, Q2, Tau, cub.method = "pcubature") {
  f.Dim <- ncol(Q1)^2

  out <- cubintegrate(
    f = integrate_Exp_vec, lower = 0, upper = Tau,
    method = cub.method, fDim = f.Dim,
    Tau = Tau, Q1 = Q1, Q2 = Q2, f.Dim = f.Dim
  )

  out1 <- matrix(out$integral, ncol = ncol(Q1))
  colnames(out1) <- rownames(out1) <- rownames(Q1)
  return(out1)
}

integrate2_Exp_C <- function(Q1, Q2, Tau, cub.method = "pcubature") {
  f.Dim <- ncol(Q1)^2

  out <- cubintegrate(
    f = integrate_Exp_C_vec, lower = 0, upper = Tau,
    method = cub.method, fDim = f.Dim,
    Tau = Tau, Q1 = Q1, Q2 = Q2, f.Dim = f.Dim
  )

  out1 <- matrix(out$integral, ncol = ncol(Q1))
  colnames(out1) <- rownames(out1) <- rownames(Q1)
  return(out1)
}


# microbenchmark(
#   integrate2_Exp(Q1=Q$Qcom, Q2=Q$Q1, Tau=1, cub.method='pcubature'),
#   integrate2_Exp_C(Q1=Q$Qcom, Q2=Q$Q1, Tau=1, cub.method='pcubature')
# )


# # this integral sums to 1
# microbenchmark(
# 1/5*1/5*2*integrate2_Exp(Q1=Q1, Q2=Q2, Tau=5),
# 1/5*1/5*2*integrate2_Exp_C2(Q1=Q1, Q2=Q2, Tau=5)
# )
#
# Tau <- 1
# microbenchmark(
#   (Pindv$P1 %*% (1/Tau^2*integrate2_Exp(Q1=Q$Q1, Q2=Q$Qcom, Tau=Tau))) + (Pindv$P2 %*% (1/Tau^2*integrate2_Exp(Q1=Q$Q2, Q2=Q$Qcom, Tau=Tau))),
#   (Pindv$P1 %*% (1/Tau^2*integrate2_Exp_C(Q1=Q$Q1, Q2=Q$Qcom, Tau=Tau) ) ) + (Pindv$P2 %*% (1/Tau^2*integrate2_Exp_C(Q1=Q$Q2, Q2=Q$Qcom, Tau=Tau))),
#
#   times = 10L
# )
