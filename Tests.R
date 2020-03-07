library("corHMM")
library("phytools")

setwd("~/Documents/Courses/Bergen_2019/PhyloCourse_Bergen_2019/PhyloCourse_Bergen-2019/R")
source('write_nexus_RB.R')

#   ____________________________________________________________________________
#   Simulate Tree and Traits                                               ####

# simulate tree using pure birth process
tree<-pbtree(n=100, scale=100, b=1, d=0)
plot(tree)

# simulate character evolution

# make rate matrix Q
Q <- mk_Qcom(rate1=c(.1,.2), rate2=c(.1,.2), names=NULL)
diag(Q) <- -rowSums(Q)
#colnames(Q) <- rownames(Q) <- c(1:4)
Q
# simulate character evolution on tree using Q
hist <- sim.history(tree, Q, nsim=1)
plot(hist)


#   ____________________________________________________________________________
#   Comparing cor and mine                                         ####

# Corr

Q.Asym <-Q*10
Q.Asym[Q.Asym<=0] <- NA

# Inference
taxa <- cbind(hist$tip.label, hist$states)
Recon_Q.Asim <- rayDISC(hist, taxa, rate.mat=Q.Asym, node.states="marginal",
                        model="ARD", root.p=c(1,0,0,0))

# infered rate matrix
Recon_Q.Asim
class(Recon_Q.Asim)
print.corhmm()

#------ Mine

Qcom.int <- mk_Qcom(rate1=c(1,2), rate2=c(1,2), names=NULL)
token.maps <- colnames(Qcom.int)
names(token.maps) <-  c(1:nrow(Qcom.int))
token.maps

char1 <-substr(taxa[,2], 1, 1)
char2 <-substr(taxa[,2], 2, 2)
chars <- cbind(char1, char2)

scen <- getEdgeScenarios(chars, focal.states=c(1,1), tree, token.maps)

#---
Qcom <- mk_Qcom(rate1=c(1,2), rate2=c(1,2), names=names(token.maps))
Qindv <- mk_Qindv(Qcom)
Pindv <- mk_Pindv(Qcom, p1=c(0.5,0,0.5,0), p2=c(0.5,0.5,0,0), p0=c(1,0,0,0))
root.p <- rep(1/4,4)


outM <- runSinba(phy=scen$tree, data=scen$data,  root.p = rep(1/4,4), Qindv=Qindv, Pindv=Pindv, edge.pattern=scen$edge.patterns[[1]],
                  p = NULL, ip = NULL, lb = 0, ub = 100,
                  cub.method='pcubature', verbose = TRUE, do.thorough=T)


#   ____________________________________________________________________________
#   Darwin Scen                                         ####

#------- 4 state corr Hidden-Extinct state models

M1 <- matrix(0, nrow=4, ncol=4)
rownames(M1) <- colnames(M1) <- c(1:4)
M1[M1==0] <- NA


Q4E.corr <- M1
Q4E.corr [1,2] <- 1
Q4E.corr [1,3] <-1
Q4E.corr [2,4] <-2
Q4E.corr [3,4] <-2
Q4E.corr

Q4E.n.corr <- M1
Q4E.n.corr [1,2] <- 1
Q4E.n.corr [1,3] <-1
Q4E.n.corr [2,4] <-1
Q4E.n.corr [3,4] <-1
Q4E.n.corr

# store simulations
Q4E.results <- matrix(0, ncol = 2, nrow = length(tree.out))
colnames(Q4E.results) <- c('Q4E.corr', 'Q4E.n.corr')

i <- 1
for (i in 1:length(tree.out)){
  tree <-tree.out[[i]]
  taxa <- cbind(names(tree$states.recode), tree$states.recode)
  taxa[,2] <- dplyr::recode(taxa[,2], '1'='1', '4'='4' )

  Q4E.results[i,1] <- greedy_rayDISC(tree, taxa, rate.mat=Q4E.corr , node.states="marginal", model="ARD", root.p=c(1,0,0,0), do.thorough=TRUE, ub=10000)$AIC
  Q4E.results[i,2] <- greedy_rayDISC(tree, taxa, rate.mat=Q4E.n.corr, node.states="marginal", model="ARD", root.p=c(1,0,0,0), do.thorough=TRUE, ub=10000)$AIC
}

Q4E.results
apply(Q4E.results, 1, function(x) which(x==min(x)))

Q4E.results[,1]-Q4E.results[,2]
#Q4E.results[,3]-Q4E.results[,2]

#------ Mine

#-- Case when Ln=1
Qcom.int <- mk_Qcom(rate1=c(1,2), rate2=c(1,2), names=NULL)
token.maps <- colnames(Qcom.int)
names(token.maps) <-  c(1:nrow(Qcom.int))
token.maps

tree <-tree.out[[2]]
char1 <-tree$states
char2 <-tree$states
chars <- cbind(char1, char2)

scen <- getEdgeScenarios(chars, focal.states=c(1,1), tree, token.maps)

dotTree(tree, x=tree$states, ftype="i",fsize=0.7)
plot(scen$tree)
#nodelabels()
edgelabels(scen$edge.events[[1]], c(1:nrow(tree$edge)), cex = 0.7)

#---
Qcom <- mk_Qcom(rate1=c(1,0), rate2=c(1,0), names=names(token.maps))
Qindv <- mk_Qindv(Qcom)
Pindv <- mk_Pindv(Qcom, p1=c(0,0,1,0), p2=c(0,1,0,0), p0=c(1,0,0,0))
root.p <- c(1, 0,0,0)


outABS <- infer_TBM(phy=scen$tree, data=scen$data,  root.p = root.p, Qindv=Qindv, Pindv=Pindv, edge.pattern=scen$edge.patterns[[1]],
                    p = NULL, ip = NULL, lb = 0, ub = 100,
                    cub.method='pcubature', verbose = TRUE, do.thorough=T)

exp(outABS$loglik)


#-- Case with absorbing model

Qcom.int <- mk_Qcom(rate1=c(1,2), rate2=c(1,2), names=NULL)
token.maps <- colnames(Qcom.int)
names(token.maps) <-  c(1:nrow(Qcom.int))
token.maps

tree <-tree.out[[2]]
char1 <-tree$states
char2 <-tree$states
chars <- cbind(char1, char2)

scen <- getEdgeScenarios(chars, focal.states=c(1,1), tree, token.maps)
class(scen)

dotTree(tree, x=tree$states, ftype="i",fsize=0.7)
plot(scen$tree)
#nodelabels()
edgelabels(scen$edge.events[[1]], c(1:nrow(tree$edge)), cex = 0.7)

#---
Qcom <- mk_Qcom(rate1=c(1,0), rate2=c(1,0), names=names(token.maps))
Qcom[2,1] <- 1
Qcom[3,1] <- 1

Qindv <- mk_Qindv(Qcom)
Pindv <- mk_Pindv(Qcom, p1=c(0,0,1,0), p2=c(0,1,0,0), p0=c(1,0,0,0))
root.p <- c(1, 0,0,0)


outABS2 <- infer_TBM(phy=scen$tree, data=scen$data,  root.p = root.p, Qindv=Qindv, Pindv=Pindv, edge.pattern=scen$edge.patterns[[1]],
                     p = NULL, ip = NULL, lb = 0, ub = 100,
                     cub.method='pcubature', verbose = TRUE, do.thorough=T)

exp(outABS2$loglik)


#--- Case3
Qcom <- mk_Qcom(rate1=c(1,1), rate2=c(1,1), names=names(token.maps))


Qindv <- mk_Qindv(Qcom)
Pindv <- mk_Pindv(Qcom, p1=c(0,0,1,0), p2=c(0,1,0,0), p0=c(1,0,0,0))
root.p <- c(1, 0,0,0)


outABS3 <- infer_TBM(phy=scen$tree, data=scen$data,  root.p = root.p, Qindv=Qindv, Pindv=Pindv, edge.pattern=scen$edge.patterns[[1]],
                     p = NULL, ip = NULL, lb = 0, ub = 100,
                     cub.method='pcubature', verbose = TRUE, do.thorough=T)

exp(outABS3$loglik)


#-----

outABS2 <- runSinba(phy=scen$tree, data=scen$data,
                    root.p = root.p, Qindv=Qindv, Pindv=Pindv, edge.pattern=scen$edge.patterns[[1]],
                    p = NULL, ip = NULL, lb = 0, ub = 100,
                    cub.method='pcubature', verbose = TRUE, do.thorough=T)

outABS2 <-runSinba(phy=scen$tree, data=scen$data,
                    root.p = root.p, Qindv=Qindv, Pindv=Pindv,
                    p = NULL, ip = NULL, lb = 0, ub = 100,
                    cub.method='pcubature', verbose = TRUE, do.thorough=T)

outABS2 <-runSinba(phy=scen$tree, data=scen$data,  edge.pattern=scen$edge.patterns[[1]],
         root.p = root.p, Qindv=Qindv, Pindv=Pindv,
         p = NULL, ip = NULL, lb = 0, ub = 100,
         cub.method='pcubature', verbose = TRUE, do.thorough=T)

out <-runSinba(Edge.Scenarios=scen,
                   root.p = NULL, Qindv=Qindv, Pindv=Pindv,
                   p = NULL, ip = NULL, lb = 0, ub = 100,
                   cub.method='pcubature', verbose = TRUE, do.thorough=T)


out
e1 <- edgePatRecode(scen$edge.events[[1]], from='codes', to='ids')

e2 <- edgePatRecode(e1, from='ids', to='codes')
cat(paste0('Uniques edge scenarios (codes): ', pasteunique(e2), '.\n' ))

out <-runSinba(phy=scen$tree, data=scen$data,  edge.pattern=scen$edge.patterns[[1]],
                   root.p = root.p, Qindv=Qindv, Pindv=Pindv,
                   p = NULL, ip = NULL, lb = 0, ub = 100,
                   cub.method='pcubature', verbose = TRUE, diagn=T, optim.method='subplex', subplex.rounds=1, subplex.n.cores=10L, gensa.its=50)


#----------------
tree<-pbtree(n=10, scale=100, b=1, d=0)
plot(tree)

# simulate character evolution

# make rate matrix Q
Q <- mk_Qcom(rate1=c(.1,.2), rate2=c(.1,.2), names=NULL)
diag(Q) <- -rowSums(Q)
#colnames(Q) <- rownames(Q) <- c(1:4)
Q
# simulate character evolution on tree using Q
hist <- sim.history(tree, Q, nsim=1)
plot(hist)
tree <- hist

#char1 <-substr(tree$states, 1, 1)
char1 <-rep(0, 10)
char1[1:2] <- 1
names(char1) <- names(tree$states)
char2 <-char1
chars <- cbind(char1, char2)

dotTree(tree, char1)

scen <- getEdgeScenarios(chars, focal.states=c(1,1), tree, token.maps)

plot.phylo(scen$tree)
edgelabels(scen$edge.events[[1]], c(1:nrow(tree$edge)), cex = 0.7)


out <-runSinba(phy=tree, data=scen$data,  edge.pattern=scen$edge.patterns[[1]],
                   root.p = c(1,0,0,0), Qindv=Qindv, Pindv=Pindv,
                   p = NULL, ip = NULL, lb = 0, ub = 100,
                   cub.method='pcubature', verbose = TRUE, do.thorough=T)


out

tr <- tree
neworder <- reorder(tr, order = "pr", index.only = TRUE)
tr$edge <- tr$edge[neworder, ]
tr$edge.length <- tr$edge.length[neworder]
attr(tr, "order") <- "pruningwise"
ee <- edge.pattern[neworder]

plot.phylo(tree)
plot.phylo(tr)
tree$edge
tr$edge

plot.phylo(out$phy)
edgelabels(scen$edge.events[[1]], c(1:nrow(tree$edge)), cex = 0.7)

dotTree(out$phy, out$data[,2])
