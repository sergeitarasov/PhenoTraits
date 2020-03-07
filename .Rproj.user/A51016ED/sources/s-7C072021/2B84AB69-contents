
#' @useDynLib PhenoTraits
#' @importFrom Rcpp sourceCpp
#' @importFrom nloptr nloptr
#' @importFrom GenSA GenSA
#' @importFrom parallel mclapply
#' @importFrom numDeriv hessian
#' @importFrom corpcor pseudoinverse
#' @importFrom plotrix draw.circle
#' @importFrom partitions setparts

usethis::use_package("ape", type = "Depends")
usethis::use_package("RcppArmadillo", "LinkingTo")
usethis::use_package("phangorn", "Import")
#usethis::use_package("nloptr", "Import")
#usethis::use_package("GenSA", "Import")

usethis::use_package("cubature", type = "Depends")
usethis::use_package("dplyr", type = "Depends")
usethis::use_package("stringr", type = "Depends")
usethis::use_package("phytools", type ="Depends")


#corpcor
#parallel
#numDeriv
# library("phytools")
# library("phangorn")
#
# library("dplyr")
# library("ape")
# setwd("~/Documents/Char_coding/Maddison_problem")
# source("SMM_functions.R")
# source("make_edge_scenarios_Functions.R")
#
# setwd("~/Documents/Macbbok_VT/RADT/Package/PhyHiFi/R")
# source("mod_ray_disc.R")
# library("corHMM")
# library("microbenchmark")
# library("RcppArmadillo")
# library("cubature")
# library("stringr")


#--- Devs
# devtools::document()

# usethis::use_testthat()
# devtools::test()

#data.sim<-readRDS(file='Chars_tree_example.rds')
#save(data.sim, file="data/data_sim.rda")
#----




#' Tree and two traits for running Sinba
#'
#' A list with two characters and a tree
#'
#' @docType data
#'
#' @usage
#' data(data.sim)
"data.sim"



# x <- scen
#' Title Print edge scenarios object
#' @export
#'
print.edge_scenarios <- function(x, ...) {
  cat("\n/--------------------------------\n")
  ntips <- Ntip(x$tree)
  nscen <- length(x$edge.patterns)

  cat("Edge Scenario object:

")
  paste0("1 Phylogenetic tree with ", ntips, " tips.
") %>% cat()
  paste0(nscen, " Edge Scenarios.
") %>% cat()
  cat("
")

  out.tok <- rbind(x$token.maps, names(x$token.maps))
  out.tok <- data.frame(out.tok)
  colnames(out.tok) <- NULL
  rownames(out.tok) <- c("original", "present")
  cat("Coding:")
  print(out.tok)

  counts <- table(x$data[, 2])
  levels <- names(counts)

  cat("\n")
  cat("State distribution:\n")
  cat("States:", levels, "\n", sep = "\t")
  cat("Counts:", counts, "\n", sep = "\t")
  cat("\\--------------------------------\n")
}


#print.sinba(out)
#x <- outABS2
#class(out)
#x <- out

#' Title Print sinba object
#'
#' @export
#'
#' @examples
print.sinba <- function(x,...){

  cat("\n/---------------------------------------------------\n")

  #token.maps, if all chars present
  #edge scenarios
  # root
  taxa <- x$data
  counts <- table(taxa[, 2])
  levels <- names(counts)

  ntips=Ntip(x$phy)
  output<-data.frame(x$loglik, x$AIC, x$AICc, ntips, row.names="")
  names(output)<-c("-lnL","AIC","AICc","ntax")
  cat('Inference with runSinba:\n')
  cat("\n")
  print(output)
  cat("\n")

  mm <- match(sort(unique(x$edge.pattern)), as.numeric(patternList()[3,]) )
  cat(paste0('Uniques edge codes: ', paste0(patternList()[2,mm ], collapse = ', '), '.\n' ))

  #---
  diagn <-diagn1 <- NULL
  if(any(x$eigval<0)){
    #If any eigenvalue is less than 0 then the solution is not the maximum likelihood solution
    diagn <- c("The objective function may be at a saddle point.")
  } else if (!all(x$eigval<0)){
    diagn <- c("Arrived at a reliable solution.")
  }

  if(x$bound.hit){
    diagn1 <-c("At least one rate parameter equals the boundary value (lb or ub).  This may be a non-optimal solution. Consider changing boundary values.\n")
  }
  #---
  if(length(paste(diagn, diagn1))>0){
    cat('Performance (', x$optim.method, '): ', paste(diagn, diagn1), '\n', sep = '')
  } else
    cat('Performance (', x$optim.method, '): ', 'diagnostics was off.\n', sep = '')


  cat("\n")

  param.est<- x$solution
  cat("Rate matrix (root.p = ", paste0(x$root.p,  collapse = ' '), '):\n', sep='')

  print(param.est)
  cat("\n")

  if(!is.null(x$token.maps)){
    out.tok <- rbind(x$token.maps, names(x$token.maps))
    out.tok <- data.frame(out.tok)
    colnames(out.tok) <- NULL
    rownames(out.tok) <- c("original", "present")
    cat("Coding:")
    print(out.tok)
  }



  counts <- table(x$data[, 2])
  levels <- names(counts)

  cat("\n")
  cat("State distribution:\n")
  cat("States:", levels, "\n", sep = "\t")
  cat("Counts:", counts, "\n", sep = "\t")

  tk <-x$data[,2]
  tk <- str_split(tk, '&' ) %>% unlist %>%  unique()

  if (!all(colnames(Qindv$Qcom) %in%  tk))
    cat(("NOT all states from Q are present in the character.\n"))

  if (all(colnames(Qindv$Qcom) %in%  tk))
    cat(("ALL states from Q are present in the character.\n"))

  cat("\\---------------------------------------------------\n")

}




