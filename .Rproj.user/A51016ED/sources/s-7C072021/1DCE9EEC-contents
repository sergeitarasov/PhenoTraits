#   ____________________________________________________________________________
#   Simulate Tree and Traits                                               ####


# tr.Darw <- simDarwin(ntips=100, n.trees.sim=10, Upper.otu.bound=0.6, Lower.otu.bound=0.4)
#
# tree <- tr.Darw[[1]]
# dat <- getEdgeScenarios(cbind(tree$states,tree$states), focal.states=c(1,1), tree, token.maps)
#
# plot.phylo(dat$tree, edge.width = 3,
#            show.tip.label = T, no.margin = F, label.offset=10, tip.color='white' )
# addTipChars(dat$data[,3:4],  r=.3,  x.space=3, legend=T)


#' Title Simulate Darwin Scenario
#'
#' @param ntips
#' @param n.trees.sim
#' @param Upper.otu.bound
#' @param Lower.otu.bound
#'
#' @return
#' @export
#'
#' @examples
simDarwin <- function(ntips=100, n.trees.sim=10, Upper.otu.bound=0.6, Lower.otu.bound=0.4){

  Upper.otu.bound<-ntips*Upper.otu.bound
  Lower.otu.bound<-ntips*Lower.otu.bound

  tree.out <- vector('list', n.trees.sim)
  tree.count <- 0

  while(tree.count<n.trees.sim){
    # simulate tree using pure birth process
    tree<-pbtree(n=ntips, scale=100, b=1, d=0)

    tips.per.node <- phangorn::Descendants(tree,
                                           node=c((ntips+1):(ntips+Nnode(tree))),
                                           type="tips")

    tips.per.node.count <- lapply(tips.per.node, length) %>% unlist

    count.boolean <- tips.per.node.count < Upper.otu.bound & tips.per.node.count > Lower.otu.bound

    if(any(count.boolean)){

      # randomly choose a clade
      rnd <- NA
      rnd <- runif(1, 1, length(tips.per.node[count.boolean])) %>% round(.,0)
      tips <- tips.per.node[count.boolean][[rnd]]

      states <- rep(0, length(tree$tip.label))
      states[tips] <- 1
      names(states) <- tree$tip.label
      tree$states <- states

      # recode statets
      XX <- cbind(paste0(tree$states, tree$states))
      XX <-recode(XX, "00"='1', "11"='4')
      names(XX) <- tree$tip.label
      tree$states.recode <- XX
      #---

      tree.count <- tree.count+1
      print(paste0('Getting tree ', tree.count))
      tree.out[[tree.count]] <- tree

    }
  }

  return(tree.out)
}



# Q <- matrix(c(-0.1, 0.1,
#     0.1, -0.1), 2,2, byrow = T)
# colnames(Q) <- rownames(Q) <- c(0,1)
#
# sim.UB <- simUnrepBurst(ntips=100, n.trees.sim=10, Upper.otu.bound=0.6, Lower.otu.bound=0.4, Q, p.root=c(1,0), freq.threshhold=3)
#
# i <- 6
# tree <- sim.UB[[i]]$sim.trinary[[1]]
# chars <- cbind(sim.UB[[i]]$sim.trinary[[1]]$states, sim.UB[[i]]$sim.trinary[[2]]$states, sim.UB[[i]]$sim.original[[1]]$states, sim.UB[[i]]$sim.original[[2]]$states)
# plot.phylo(tree, edge.width = 3, show.tip.label = T, no.margin = F, label.offset=10, tip.color='white' )
#addTipChars(chars,  r=.4,  x.space=3, legend=T)

#' Title Simulate Unreolicated Burst
#'
#' @param ntips
#' @param n.trees.sim
#' @param Upper.otu.bound
#' @param Lower.otu.bound
#' @param Q
#' @param p.root
#' @param freq.threshhold frequency of cocuring of each state inside the clade where chars are born.This trees that have this value < freq.threshhold are filtered out.
#'
#' @return
#' @export
#'
#' @examples
# simUnrepBurst <- function(ntips=100, n.trees.sim=10, Upper.otu.bound=0.6, Lower.otu.bound=0.4, Q, p.root=c(1,0), freq.threshhold=3){
#
#   names(p.root) <- rownames(Q)
#   Upper.otu.bound<-ntips*Upper.otu.bound
#   Lower.otu.bound<-ntips*Lower.otu.bound
#
#   tree.out <- vector('list', n.trees.sim)
#   tree.count <- 0
#
#   while(tree.count<n.trees.sim){
#
#     tree<-pbtree(n=ntips, scale=100, b=1, d=0)
#     hist <- hist2 <- hist3 <- sim.history(tree, Q, anc=p.root, nsim=2)
#
#     tips.per.node <- phangorn::Descendants(tree,
#                                            node=c((ntips+1):(ntips+Nnode(tree))),
#                                            type="tips")
#
#     tips.per.node.count <- lapply(tips.per.node, length) %>% unlist
#     count.boolean <- tips.per.node.count < Upper.otu.bound & tips.per.node.count > Lower.otu.bound
#
#     if(any(count.boolean)){
#
#       # randomly choose a clade
#       rnd <- NULL
#       rnd <- runif(1, 1, length(tips.per.node[count.boolean])) %>% round(.,0)
#
#       tips <- tips.per.node[count.boolean][[rnd]]
#
#       hist2[[1]]$states[-tips] <- 0
#       hist3[[1]]$states[-tips] <- 2
#
#       hist2[[2]]$states[-tips] <- 0
#       hist3[[2]]$states[-tips] <- 2
#
#       freq_1 <- table(hist2[[1]]$states[tips])
#       freq_2 <- table(hist2[[2]]$states[tips])
#       freq <- rbind(freq_1, freq_2)
#
#       if (all(freq_1>freq.threshhold) & all(freq_1>freq.threshhold)){
#         tree.count <- tree.count+1
#         print(paste0('Getting tree ', tree.count))
#         tree.out[[tree.count]] <- list(sim.original=hist, sim.binary=hist2, sim.trinary=hist3, internal.state.freq = freq)
#       }
#     }
#   }
#
#   return(tree.out)
# }

#_------------------------------


# token.maps <- c("00", "01", "10", "11")
# names(token.maps) <-  c(1:4)
# Q <- mk_Qcom(rate1=c(.1,.1), rate2=c(.1,.1), names=c(1:4))
# diag(Q) <- -rowSums(Q)
# Q
#
# sim.UB <- simUnrepBurst(ntips=100, n.trees.sim=10, Upper.otu.bound=0.6, Lower.otu.bound=0.4, Q, p.root=c(1,0,0,0), freq.threshhold=3)
#
# i <- 2
# tree <- sim.UB[[i]]$tree
# chars <- cbind(sim.UB[[i]]$states.binary[[1]], sim.UB[[i]]$states.binary[[2]], sim.UB[[i]]$states.original, sim.UB[[i]]$states.original.recoded)
# plot.phylo(tree, edge.width = 3, show.tip.label = T, no.margin = F, label.offset=10, tip.color='white' )
# addTipChars(chars,  r=.4,  x.space=3, legend=T)


#' Title Simulate Unreplicated Burst
#'
#' @param ntips
#' @param n.trees.sim
#' @param Upper.otu.bound
#' @param Lower.otu.bound
#' @param Q
#' @param p.root
#' @param freq.threshhold
#'
#' @return
#' @export
#'
#' @examples
simUnrepBurst <- function(ntips=100, n.trees.sim=10, Upper.otu.bound=0.6, Lower.otu.bound=0.4, Q, p.root=c(1,0,0,0), freq.threshhold=5){

  names(p.root) <- rownames(Q)
  Upper.otu.bound<-ntips*Upper.otu.bound
  Lower.otu.bound<-ntips*Lower.otu.bound

  tree.out <- vector('list', n.trees.sim)
  tree.count <- 0

  while(tree.count<n.trees.sim){

    tree<-pbtree(n=ntips, scale=100, b=1, d=0)
    #hist <- hist2 <- hist4 <- sim.history(tree, Q, anc=p.root, nsim=2)
    hist  <- sim.history(tree, Q, anc=p.root, nsim=1)

    tips.per.node <- phangorn::Descendants(tree,
                                           node=c((ntips+1):(ntips+Nnode(tree))),
                                           type="tips")

    tips.per.node.count <- lapply(tips.per.node, length) %>% unlist
    count.boolean <- tips.per.node.count < Upper.otu.bound & tips.per.node.count > Lower.otu.bound

    if(any(count.boolean)){

      # randomly choose a clade
      rnd <- NULL
      rnd <- runif(1, 1, length(tips.per.node[count.boolean])) %>% round(.,0)

      tips <- tips.per.node[count.boolean][[rnd]]

      states.original.recoded <- states.original <- hist$states
      #states.original.recoded <- states.original <- list(hist[[1]]$states, hist[[2]]$states)

      # make Unrep Burst
      states.original.recoded[-tips] <- 1
      #states.original.recoded[[2]][-tips] <- 1

      states.binary <-list()
      # Recode into binary
      # char 1
      states.binary[[1]] <- dplyr::recode(states.original.recoded, '1'='0', '2'='0', '3'='1', '4'='1')
      names(states.binary[[1]]) <- names(states.original.recoded)
      # char 2
      states.binary[[2]] <-dplyr::recode(states.original.recoded, '1'='0', '2'='1', '3'='0', '4'='1')
      names(states.binary[[2]]) <- names(states.original.recoded)

      # Claculate frequencies
      #freq_1 <- table(states.original[[1]][tips])
      #freq_2 <- table(states.original[[2]][tips])
      freq <- table(states.original[tips])

        if (all(freq > freq.threshhold) & length(freq)==4){
          tree.count <- tree.count+1
          print(paste0('Getting tree ', tree.count))
          tree.out[[tree.count]] <- list(sim.original=hist, tree=tree, states.original=states.original,
                                         states.original.recoded = states.original.recoded, states.binary=states.binary, internal.state.freq = freq)
        }

    }
  }

  return(tree.out)
}


