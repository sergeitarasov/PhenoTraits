#   ____________________________________________________________________________
#   Create Rate Matrices                                                    ####



#' Title
#'
#' @param rate1
#' @param rate2
#' @param names
#'
#' @return
#' @export
#'
#' @examples
mk_Qcom <- function(rate1 = c(1, 2), rate2 = c(3, 4), names = NULL) {
  char.state <- c(0, 1)
  Q1s <- init_char_matrix(char.state, rate1, diag.as = 0)
  Q2s <- init_char_matrix(char.state, rate2, diag.as = 0)

  out <- comb2matrices(Q1s, Q2s, controlling.state = NULL, name.sep = "", diag.as = "", non.rate.as = NULL)
  if (!is.null(names))
    colnames(out) <- rownames(out) <- names
  else colnames(out) <- rownames(out) <-c(1:4)

  return(out)
}

#mk_QcomPartitions(num.elenents=4)
#' Title
#'
#' @param num.elenents
#'
#' @return
#' @export
#'
#' @examples
mk_QcomPartitions <- function(num.elenents=4, names=NULL){
  Part <- setparts(num.elenents)
  lapply(1:ncol(Part), function(x) mk_Qcom(rate1=Part[c(1,2), x], rate2=Part[c(3,4), x], names=names))
}

#' Title Get all parametrizations of Q with 8 elements
#'
#' @return
#' @export
#'
#' @examples
mk_QcomPartitions8 <- function(){

  Part <- setparts(8)
  Qcom <- mk_Qcom(rate1=c(1,1), rate2=c(1,1))
  Qcom[Qcom==1] <- c(1:8)
  Qcom[Qcom==0] <- max(Qcom)+1
  #dim(Part)

  Q.z <- matrix(0,4,4)
  #Q.z[] <- c(Part[,2],0)[Qcom]
  Part.list <- lapply(seq_len(ncol(Part)), function(i) Part[,i])
  Q.list <- lapply(Part.list, function(x) {Q.z <- matrix(0,4,4); Q.z[] <- c(x,0)[Qcom]; colnames(Q.z) <- rownames(Q.z) <- c(1:4); Q.z} )
  return(Q.list)

}


#' Title
#'
#' @param rate1
#' @param rate2
#' @param names
#'
#' @return
#' @export
#'
#' @examples
mk_Qindv <- function(Qcom) {
  Q1 <- Q2 <- Qcom
  Q1[] <- Q2[] <- 0

  Q1[1, 3] <- Qcom[1, 3]
  Q1[3, 1] <- Qcom[3, 1]

  Q2[1, 2] <- Qcom[1, 2]
  Q2[2, 1] <- Qcom[2, 1]

  return(list(Q1 = Q1, Q2 = Q2, Qcom = Qcom))
}


#' Title
#'
#' @param rate1
#' @param rate2
#' @param names
#'
#' @return
#' @export
#'
#' @examples
mk_QindvZero <- function(Qindv) {
  lapply(Qindv, function(x) {
    x[] <- 0
    x
  })
}


#' Title
#'
#' @param rate1
#' @param rate2
#' @param names
#'
#' @return
#' @export
#'
#' @examples
mk_Pindv <- function(Qcom, p1 = c(0.5, 0, 0.5, 0), p2 = c(0.5, 0.5, 0, 0), p0 = c(1, 0, 0, 0)) {
  P1 <- P2 <- P0 <- Qcom
  P1[] <- P2[] <- P0[] <- 0

  P1[1, ] <- p1
  P2[1, ] <- p2
  P0[1, ] <- p0

  return(list(P1 = P1, P2 = P2, P0 = P0))
}

is.sequential <- function(x, eps = 1e-8) {
  if (length(x) && isTRUE(abs(x[1] - floor(x[1])) < eps)) {
    all(abs(diff(x) - 1) < eps)
  } else {
    FALSE
  }
}


#' Title
#'
#' @param rate1
#' @param rate2
#' @param names
#'
#' @return
#' @export
#'
#' @examples
mk_RatePars <- function(Qindv) {
  uniq.rt <- lapply(Qindv, function(x) c(x)) %>%
    unlist() %>%
    unique()
  uniq.rt <- uniq.rt[uniq.rt > 0]

  zero.rate <- max(uniq.rt + 1)

  if (!is.sequential(sort(uniq.rt))) stop("Rate paremeters are not sequential")

  out <- lapply(Qindv, function(x) {
    x[x == 0] <- zero.rate
    x
  })

  return(out)
}


#' Title
#'
#' @param rate1
#' @param rate2
#' @param names
#'
#' @return
#' @export
#'
#' @examples
mk_OneRate <- function(Qindv) {
  out <- lapply(Qindv, function(x) {
    x[x > 0] <- 1
    x[x == 0] <- 2
    x
  })

  return(out)
}


#' Title Create new Qcom by merging paramter estimates that are similar, and setting 0 estmates to 0
#'
#' @param pars
#' @param Qcom.cor Qcom
#'
#' @return
#' @export
#'
#' @examples
#iter_new_Qcom(pars=out$est.pars, Qcom.cor=Qindv$Qcom)
iter_new_Qcom <- function(pars, Qcom.cor){

  # remove zero pars
  pars0 <- which(pars==0)

  if (length(pars0>0)){

    for (i in 1:length(pars0)){
      Qcom.cor[Qcom.cor==pars0[i]] <- 0
    }
  } else {
    #-- distance
    dd <- dist(pars)
    ij <- finv(which(dd==min(dd)), dd)
    ij <- ij[1,]

    max.state <- max(Qcom.cor)
    Qcom.cor[Qcom.cor==ij[1]] <- max.state+1
    Qcom.cor[Qcom.cor==ij[2]] <- max.state+1
  }

  # make sequential
  all <- Qcom.cor[Qcom.cor>0]
  unq <- Qcom.cor[Qcom.cor>0] %>% unique()
  unq.new <- c(1:length(unq))
  new <- unq.new [match(all, unq)]
  Qcom.cor[Qcom.cor>0] <- new

  return(Qcom.cor)
}


## 1D index to 2D index
finv <- function (k, dist_obj) {
  if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object")
  n <- attr(dist_obj, "Size")
  valid <- (k >= 1) & (k <= n * (n - 1) / 2)
  k_valid <- k[valid]
  j <- rep.int(NA_real_, length(k))
  j[valid] <- floor(((2 * n + 1) - sqrt((2 * n - 1) ^ 2 - 8 * (k_valid - 1))) / 2)
  i <- j + k - (2 * n - j) * (j - 1) / 2
  cbind(i, j)
}

