#   ____________________________________________________________________________
#   Simulate Unreplicated Burst                                             ####


token.maps <- c("00", "01", "10", "11")
names(token.maps) <-  c(1:4)
Q <- mk_Qcom(rate1=c(.1,.1), rate2=c(.1,.1))
diag(Q) <- -rowSums(Q)
Q

# Simulate Darwin scenario (unreplicated burst)
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

# Vary i to observe various edge patterns
i=2
EE <- scen$edge.events[[i]]
edge.color <-edgePatRecode(EE, from='codes', to='cols')
dat <- scen$data[,2:4]
rownames(dat) <- scen$data[,1]
plot.phylo(scen$tree, edge.color = edge.color, edge.width = 3, show.tip.label = T, no.margin = T, label.offset=10, tip.color='white' )
addTipChars(x=dat, colors, r=.4,  x.space=3, legend=TRUE)
edgelabels(EE, cex = 0.4, frame = "ci", bg = edge.color)




#---------- Run inference
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

