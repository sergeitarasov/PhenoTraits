
#   ____________________________________________________________________________
#   Make Edge Scenarios                                                     ####



mrca_edges <- function(char1, tree, focal.state = 1) {
  char1.1 <- which(char1 == focal.state)
  mrca.char1 <-
    phangorn::mrca.phylo(tree, node = char1.1, full = FALSE)
  mrca.char1 <-
    c(
      mrca.char1,
      phangorn::Ancestors(tree, node = mrca.char1, type = c("all"))
    )
  mrca.edges <- match(mrca.char1, tree$edge[, 2])
  mrca.edges <- mrca.edges[!is.na(mrca.edges)]

  # root is -1
  if (length(mrca.edges) == 0) {
    return(-1)
  }

  # add root
  mrca.edges <- c(-1, mrca.edges)
  return(mrca.edges)
}

# getEdgeProps_one(focal.comb, tree)
getEdgeProps_one <- function(focal.comb, tree) {
  edge.props <- vector(mode = "list", length = nrow(tree$edge))

  if (!any(focal.comb == -1)) {
    ANC <-
      phangorn::Ancestors(tree, c(tree$edge[focal.comb[1, 1], 2], tree$edge[focal.comb[1, 2], 2]), type = c("all"))
    ANC <- lapply(ANC, length) %>% unlist()
    char.order <- order(ANC)

    # i <- 1
    for (i in char.order) {
      # edge.props[[focal.comb[1,i]]] <- c(edge.props[[focal.comb[1,i]]], -1*i)
      edge.props[[focal.comb[1, i]]] <-
        c(edge.props[[focal.comb[1, i]]], paste0("*", i))
      des <-
        phangorn::Descendants(tree, c(tree$edge[focal.comb[1, i], 2]), type = c("all"))

      des <- match(des, tree$edge[, 2])
      # edge.props[[des]] <- paste0(edge.props[[des]], ' c_', i)
      edge.props[des] <- lapply(edge.props[des], function(x)
        c(x, i))
    }

    # root is -1
  } else if (length(which(focal.comb == -1)) == 1) {
    not.i <- which(focal.comb == -1)
    edge.props <- lapply(edge.props, function(x)
      x <- not.i)

    i <- which(focal.comb != -1)
    # edge.props[[focal.comb[1,i]]] <- c(edge.props[[focal.comb[1,i]]], -1*i)
    edge.props[[focal.comb[1, i]]] <-c(edge.props[[focal.comb[1, i]]], paste0("*", i))

    des <-phangorn::Descendants(tree, c(tree$edge[focal.comb[1, i], 2]), type = c("all"))

    des <- match(des, tree$edge[, 2])
    # edge.props[[des]] <- paste0(edge.props[[des]], ' c_', i)
    edge.props[des] <- lapply(edge.props[des], function(x)
      c(x, i))
  } else {
    edge.props <- lapply(edge.props, function(x)
      x <- c(1, 2))
  }


 edge.props <-
    lapply(edge.props, function(x)
      if (is.null(x)) {
        x <- 0
      } else {
        x
      })
  edge.props.conc <-
    lapply(edge.props, function(x)
      paste0(x, collapse = "")) %>% unlist()

  return(edge.props.conc)
}

# chars=cbind(char1.masked, char2.masked)
getEdgeProps <- function(chars, focal.states = c(1, 1), tree) {
  mrca.edges1 <- mrca_edges(chars[, 1], tree, focal.state = focal.states[1])
  mrca.edges2 <- mrca_edges(chars[, 2], tree, focal.state = focal.states[2])

  combs <- expand.grid(mrca.edges1, mrca.edges2, stringsAsFactors = F)

  out <- vector(length = nrow(combs), mode = "list")
  # i <- 2
  for (i in 1:nrow(combs)) {
    focal.comb <- combs[i, ]
    out[[i]] <- getEdgeProps_one(focal.comb, tree)
  }

  return(out)
}


# chars <- cbind(char1, char2)
# chars[1,1] <- '0&1'
# chars[3,1] <- '0&1'
# getEdgeScenarios(chars, focal.states=c(1,1), tree, token.maps)

#' Title Create edge scenarios
#'
#' @param chars
#' @param focal.states
#' @param tree
#' @param token.maps
#'
#' @return
#' @export
#'
#' @examples
getEdgeScenarios <- function(chars, focal.states = c(1, 1), tree, token.maps) {

  # Preperation
  tree <- reorder(tree, "pruningwise")

  taxa <- rownames(chars)
  chars <- chars[match(taxa, tree$tip.label), ]

  # mask ambiguos chars as focal.states
  char1.masked <- chars[, 1]
  char1.masked[str_which(char1.masked, "&")] <- as.character(focal.states[1])

  char2.masked <- chars[, 2]
  char2.masked[str_which(char2.masked, "&")] <- as.character(focal.states[2])

  #---
  Births <- getEdgeProps(cbind(char1.masked, char2.masked), focal.states, tree)
  edge.patterns <- get_edgePattern(Births)

  #---
  chars <- apply(chars, 1, function(x) recode2chars_row(x, token.maps)) %>% cbind(taxa = rownames(chars), c12 = ., chars)
  rownames(chars) <- NULL
  #---
  out <- list(edge.events = Births, edge.patterns = edge.patterns, token.maps = token.maps, tree = tree, data = chars)

  class(out) <-c(class(out), 'edge_scenarios')

  return(out)
}


#' Title Get a table of edge scenarios and their ids
#'
#' @return
#' @export
#'
#' @examples
patternList <- function() {
  codes <- c("*1", "*2", "1*2", "2*1", "*1*2", "*2*1", "21", "12", "1", "2", "0")
  names <- c("P1Q1", "P2Q2", "Q1Qc", "Q2Qc", "P1Q1Qc", "P2Q2Qc", "Qc", "Qc", "Q1", "Q2", "P00")
  ids <- c(1:4, 5, 5, 6, 6, 7:9)
  cols <- c("#E41A1C", "#E41A1C", "#E41A1C", "#E41A1C", "#E41A1C", "#E41A1C", "#4DAF4A", "#4DAF4A", "#377EB8", "#377EB8", "#999999")

  rbind(names, codes, ids, cols)

  # library("RColorBrewer")
  # mypalette<-brewer.pal(9,"Set1")
  # cols4 <- c("#999999", "#4DAF4A",  "#E41A1C", "#377EB8")
  # col4.ids <- c(rep(3,6), 2,2, 4,4, 1)
  # pat <- patternList()
  # pat <-rbind(pat, cols4[col4.ids])
  # row.names(pat)[4] <- 'col'
}


#' Title Recode edge props given patternList()
#'
#' @param vec
#' @param from
#' @param to
#'
#' @return
#' @export
#'
#' @examples
#' # edgePatRecode(vec, from='codes', to='ids')
edgePatRecode <- function(vec, from='codes', to='ids'){
  pat <- patternList()
  from <- which(from==rownames(pat))
  to <- which(to==rownames(pat))

  out <- pat[to,][match(vec, pat[from,])]
  if (to==3) out <- as.numeric(out)
  return(out)
}



get_edgePattern.character <- function(edge.pattern) {
  pat.ids <- patternList()[3, ]
  pat.vals <- patternList()[2, ]

  as.numeric(pat.ids[match(edge.pattern, pat.vals)])
}
# get_edgePattern_one(Births[[3]])

get_edgePattern.list <- function(edge.pattern) {
  lapply(edge.pattern, function(x) get_edgePattern.character(x))
}

get_edgePattern <- function(edge.pattern, ...) {
  UseMethod("get_edgePattern", edge.pattern)
}

#--- Recode 2 chars

# Qcom <- mk_Qcom(rate1=c(1,2), rate2=c(3,4), names=NULL)
# token.maps <- colnames(Qcom)
# names(token.maps) <-  c(1:nrow(Qcom))
#
# dt <- getEdgeScenarios(chars, focal.states=c(1,1), tree)
# chars <- dt$data

# row.focal <- chars[1,]
# recode2chars_row(row.focal, token.maps)
# apply(chars, 1, function(x) recode2chars_row(x, token.maps))

recode2chars_row <- function(row.focal, token.maps) {
  X <- str_split(row.focal, "&")
  X <- expand.grid(X, stringsAsFactors = FALSE)
  X <- apply(X, 1, function(x) paste(x, collapse = ""))
  X <- token.maps[match(X, token.maps)] %>% names()
  paste(X, collapse = "&")
}

#' @title Initialize binary matrices
#' @param char.state names for character states
#' @param rate.param names for the rate parameters
#' @param diag.as values to pas to the main diagonal elements
#' @return matrix
#' @export
init_char_matrix <- function(char.state, rate.param, diag.as = NA) {
  n.state <- length(char.state)
  Q <- matrix(ncol = n.state, nrow = n.state, byrow = TRUE)
  Q[xor(lower.tri(Q, diag = FALSE), upper.tri(Q, diag = FALSE))] <- rate.param
  Q <- t(Q)
  diag(Q) <- diag.as
  rownames(Q) <- colnames(Q) <- as.character(char.state)
  return(Q)
}

#' @title Combining two matrices
#' @description Combining two matrices. The parametric schem of matrice is defined by nattural
#' numbers; different numbers = different rate parameters
#' @param M1 matrix; if dependency true thenM1 controls M2
#' @param M2 matrix; if dependency true then: M2 depends on those states of M1 specified in controlling.state
#' @param controlling.state state(s) of M1 that switches on/off matrix M2
#' @param name.sep separator for state names
#' @param diag.as hpopulate main diagonal with
#' @return Matrix
#' @export
#' @examples
#' M1 <- matrix(c(-1, 1, 2, -2), 2, 2, byrow = TRUE)
#' rownames(M1) <- colnames(M1) <- c("0", "1")
#' M2 <- matrix(c(-3, 3, 4, -4), 2, 2, byrow = TRUE)
#' rownames(M2) <- colnames(M2) <- c("0", "1")
#' comb2matrices(M1, M2, controlling.state = NULL)
#' comb2matrices(M1, M2, controlling.state = 2)
comb2matrices <- function(M1, M2, controlling.state = NULL, name.sep = "", diag.as = "", non.rate.as = NULL) {
  if (!is.null(controlling.state)) {
    matrix.diag <- rep(0, ncol(M1))
    matrix.diag[controlling.state] <- 1
    matrix.diag <- diag(matrix.diag)
  }

  if (is.null(controlling.state)) {
    matrix.diag <- diag(nrow(M1))
  }

  M_kr <- (M1 %x% diag(nrow(M2))) + (matrix.diag %x% M2)


  # getting colnames

  col <- paste(colnames(kronecker(M1, diag(nrow(M2)), make.dimnames = T)),
    colnames(kronecker(diag(nrow(M1)), M2, make.dimnames = T)),
    sep = ""
  )
  col <- gsub("::", name.sep, col, fixed = T)

  # merging two names
  rownames(M_kr) <- colnames(M_kr) <- col
  if (diag.as != "") diag(M_kr) <- diag.as
  if (!is.null(non.rate.as)) M_kr[which(M_kr <= 0)] <- non.rate.as

  return(M_kr)
}



##  ............................................................................
##  Add character states to tree plot                                       ####


#' Title Add chars to plot
#'
#' @param x
#' @param colors
#' @param r
#' @param x.space
#' @param legend
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' # patternList()
#' # EE <- scen$edge.events[[3]]
#' # edge.color <-edgePatRecode(EE, from='codes', to='cols')
#' #
#' # dat <- scen$data[,3:4]
#' # rownames(dat) <- scen$data[,1]
#' #
#' # plot.phylo(scen$tree, edge.color = edge.color, edge.width = 3, show.tip.label = T, no.margin = F, label.offset=10 )
#' # addTipChars(x=dat, colors, r=1,  x.space=3, legend=TRUE)
#' # edgelabels(EE, cex = 0.7, frame = "ci", bg = edge.color)
addTipChars <- function(x, r=1,  x.space=5, legend=TRUE, ...){

  if (is.data.frame(x) | is.character(x))
    x <- as.matrix(x)
  if (hasArg(colors)) {
    colors <- list(...)$colors
  } else {
    ss <- unique(as.vector(x))
    colors <- setNames(palette()[1:length(ss)], ss)
  }

  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  x.tip <- obj$xx[1:obj$Ntip]
  y.tip <- obj$yy[1:obj$Ntip]

  for (i in 1:ncol(x)) {
    nulo <- mapply(draw.circle, x = x.tip + 1.2 * strwidth("W") +
                     x.space * (i - 1), y = y.tip, col = colors[as.character(x[,i])], MoreArgs = list(nv = 20, radius = r))
  }

  if (legend) {
    phytools::add.simmap.legend(colors = colors, prompt = FALSE,
                                vertical = FALSE, shape = "circle", x = -0.45,
                                y = -0.06)
  }

}


