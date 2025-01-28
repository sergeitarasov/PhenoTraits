# #--------------


# pat <- patternList()
# edge.color <- pat[4,][match(scen$edge.events[[1]], pat[2,])]

patternList()
EE <- scen$edge.events[[3]]
edge.color <-edgePatRecode(EE, from='codes', to='cols')
plot.phylo(scen$tree, edge.color = edge.color, edge.width = 3)
edgelabels(EE, cex = 0.7, frame = "ci", bg = edge.color)

patternList()
EE <- scen$edge.events[[3]]
edge.color <-edgePatRecode(EE, from='codes', to='cols')

dat <- scen$data[,3:4]
rownames(dat) <- scen$data[,1]

plot.phylo(scen$tree, edge.color = edge.color, edge.width = 3, show.tip.label = T, no.margin = F, label.offset=10 )
addTipChars(x=dat, colors, r=1,  x.space=3, legend=TRUE)
edgelabels(EE, cex = 0.7, frame = "ci", bg = edge.color)

dat<- scen$data[,2:3]
names(dat) <- scen$data[,1]
plot.phylo(scen$tree, edge.color = edge.color, edge.width = 3, show.tip.label = T, no.margin = F, label.offset=10 )
addTipChars(x=dat, colors, r=1,  x.space=3, legend=TRUE)
edgelabels(EE, cex = 0.7, frame = "ci", bg = edge.color)

plot.phylo(scen$tree, edge.color = edge.color, edge.width = 3, show.tip.label = T, no.margin = F, label.offset=10 )
addTipChars(x=scen$data[,2:4], colors, r=1,  x.space=3, legend=TRUE)
edgelabels(EE, cex = 0.7, frame = "ci", bg = edge.color)
#----

token.maps <- c("11", "12", "21", "22")
names(token.maps) <-  c(1:4)
token.maps

tree <-data.sim$tree
char1 <- c(rep(1,5), rep(2,5))
names(char1) <-names(data.sim$char1)
char2 <- char1
chars <- cbind(char1, char2)

scen1 <- getEdgeScenarios(chars, focal.states=c(2,2), tree, token.maps)
scen1$data

token.maps <- c("00", "01", "10", "11")
names(token.maps) <-  c(1:4)
token.maps

tree <-data.sim$tree
char1 <-data.sim$char1
char2 <-data.sim$char2
chars <- cbind(char1, char2)

scen <- getEdgeScenarios(chars, focal.states=c(1,1), tree, token.maps)

dotTree(tree, x=chars, ftype="i",fsize=0.7, method = 'plot.phylo')
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



out <- runSinba(phy=scen$tree, data=scen$data, edge.pattern=scen$edge.patterns[[1]],
                root.p = NULL, Qindv=Qindv, Pindv=Pindv,
                p = NULL, ip = NULL, lb = 0, ub = 100,
                cub.method='pcubature', verbose = TRUE, diagn=FALSE, optim.method='subplex', subplex.rounds=1, subplex.n.cores=1L, gensa.its=100)

out <- runSinba(phy=scen$tree, data=scen$data, edge.pattern=NULL,
                root.p = NULL, Qindv=Qindv, Pindv=Pindv,
                p = NULL, ip = NULL, lb = 0, ub = 100,
                cub.method='pcubature', verbose = TRUE, diagn=T, optim.method='subplex_annl', subplex.rounds=1, subplex.n.cores=1L, gensa.its=100)


out <- runSinba(phy=scen$tree, data=scen$data, edge.pattern=NULL,
                root.p = NULL, Qindv=Qindv, Pindv=Pindv,
                p = NULL, ip = NULL, lb = 0, ub = 100,
                cub.method='pcubature', verbose = TRUE, diagn=T, optim.method='subplex_annl', subplex.rounds=100, subplex.n.cores=10L, gensa.its=100)


out$solution.subplex.list
out$solution.se
out$solution.low95
out$solution.up95

exp(out$loglik)
str(out)

out2 <- runSinba(phy=scen$tree, data=scen$data, edge.pattern=scen$edge.patterns[[1]],
                 root.p = NULL, Qindv=Qindv, Pindv=Pindv,
                 p = NULL, ip = NULL, lb = 0, ub = 100,
                 cub.method='pcubature', verbose = TRUE, diagn=FALSE, optim.method='subplex', subplex.rounds=5, subplex.n.cores=5L, gensa.its=100)

out2$solution.subplex.list[1,] %>% round(., 3)

out3 <- runSinba(phy=scen$tree, data=scen$data, edge.pattern=scen$edge.patterns[[1]],
                 root.p = NULL, Qindv=Qindv, Pindv=Pindv,
                 p = NULL, ip = NULL, lb = 0, ub = 100,
                 cub.method='pcubature', verbose = TRUE, diagn=FALSE, optim.method='subplex_annl', subplex.rounds=1, subplex.n.cores=1L, gensa.its=10)


