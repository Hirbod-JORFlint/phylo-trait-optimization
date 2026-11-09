library(ape)



set.seed(42)

nTips <- 500

my_tree <- rtree(
  nTips,
  tip.label=paste0("Tip",1:nTips)
)


traits <- matrix(
  runif(nTips*10,0,5),
  nrow=nTips,
  dimnames=list(
    my_tree$tip.label,
    paste0("Trait",1:10)
  )
)

nTips <- 10

my_tree <- rtree(
  nTips,
  tip.label = paste0("Tip", 1:nTips)
)

traits <- matrix(
  runif(nTips * 10, 0, 5),
  nrow = nTips,
  dimnames = list(
    my_tree$tip.label,
    paste0("Trait", 1:10)
  )
)

get_layout_cost <- function(tree, traits) {
  n <- length(tree$tip.label)
  
  # 1. Vectorized Traversal (replaces recursive DFS)
  # Extracts tips in the exact order they appear in the edge matrix.
  # *CRITICAL*: Tree must be in "cladewise" order (see note below).
  is_tip <- tree$edge[, 2] <= n
  ord <- tree$tip.label[tree$edge[is_tip, 2]]
  
  # 2. Optimized Distance Calculation
  X <- traits[ord, , drop = FALSE]
  
  # diff() is a native C-level function that does exactly what X[-1,] - X[-nrow,] does, but faster.
  sum(sqrt(rowSums(diff(X)^2)))
}

get_layout_cost2 <- function(tree, traits, scale_traits = TRUE) {
  n <- length(tree$tip.label)
  
  # Make sure edge order corresponds to the displayed/rotated cladewise order
  tree <- reorder(tree, order = "cladewise")
  
  is_tip <- tree$edge[, 2] <= n
  ord_nodes <- tree$edge[is_tip, 2]
  ord_labels <- tree$tip.label[ord_nodes]
  
  # Important: scale exactly the same way as make_trait_dist()
  X_all <- traits[tree$tip.label, , drop = FALSE]
  
  if (scale_traits) {
    X_all <- scale(X_all)
  }
  
  X <- X_all[ord_labels, , drop = FALSE]
  
  sum(sqrt(rowSums(diff(X)^2)))
}

# ---------------------------------------------------------------
# children[[node]] = integer vector of child node ids, for every
# node (tip or internal) in an ape "phylo" object.
# ---------------------------------------------------------------
get_children_list <- function(phy) {
  n_tip  <- length(phy$tip.label)
  n_node <- phy$Nnode
  children <- vector("list", n_tip + n_node)
  for (i in seq_len(nrow(phy$edge))) {
    p  <- phy$edge[i, 1]
    ch <- phy$edge[i, 2]
    children[[p]] <- c(children[[p]], ch)
  }
  children
}

# ---------------------------------------------------------------
# olo_phylo(phy, d)
#
#   phy : a fully bifurcating ape "phylo" object (n tips). Branch
#         lengths are irrelevant -- only topology matters.
#   d   : n x n distance matrix / dist object over the TRAITS,
#         with row/col i corresponding to phy$tip.label[i].
#
# Returns list(order, cost):
#   order = integer vector (permutation of 1:n) giving the leaf
#           order that minimizes summed adjacent trait-distance,
#           subject to the tree's topology.
#   cost  = the corresponding minimal total adjacent dissimilarity.
#
# Algorithm: Bar-Joseph dynamic programming. For every internal
# node v with children L, R, and for every pair of leaves (u,w)
# with u under L and w under R, define
#     M(u,w) = min_{b in leaves(L), c in leaves(R)}
#                 [ M_L(u,b) + M_R(c,w) + d(b,c) ]
# i.e. the cheapest way to lay out v's leaves starting at u and
# ending at w, where b/c are whichever leaves end up adjacent
# across the L/R seam. Because L's and R's leaf blocks can each be
# internally flipped independently, this exactly enumerates every
# topology-respecting leaf arrangement.
# ---------------------------------------------------------------
olo_phylo <- function(phy, d) {
  
  n_tip <- length(phy$tip.label)
  if (!is.matrix(d)) d <- as.matrix(d)
  stopifnot(nrow(d) == n_tip, ncol(d) == n_tip)
  
  children <- get_children_list(phy)
  memo <- vector("list", n_tip + phy$Nnode)
  
  dp <- function(node) {
    if (!is.null(memo[[node]])) return(memo[[node]])
    kids <- children[[node]]
    
    if (is.null(kids)) {                       # leaf / base case
      res <- list(leaves = node,
                  table = matrix(0, 1, 1, dimnames = list(node, node)))
      memo[[node]] <<- res
      return(res)
    }
    if (length(kids) != 2) {
      stop("Node ", node, " has ", length(kids), " children; the tree ",
           "must be fully bifurcating. Resolve polytomies first with ",
           "ape::multi2di().")
    }
    
    left  <- dp(kids[1])
    right <- dp(kids[2])
    Lleaves <- left$leaves;  Rleaves <- right$leaves
    LT <- left$table;        RT <- right$table
    
    # Re-embed each child's table as a full symmetric matrix over its
    # own leaves, so we can look up M(u,b) no matter which of the
    # child's own two sub-children u and b originally belonged to.
    L_full <- matrix(Inf, length(Lleaves), length(Lleaves),
                     dimnames = list(Lleaves, Lleaves))
    if (length(Lleaves) == 1) L_full[1,1] <- 0 else {
      L_full[rownames(LT), colnames(LT)] <- LT
      L_full[colnames(LT), rownames(LT)] <- t(LT)
    }
    R_full <- matrix(Inf, length(Rleaves), length(Rleaves),
                     dimnames = list(Rleaves, Rleaves))
    if (length(Rleaves) == 1) R_full[1,1] <- 0 else {
      R_full[rownames(RT), colnames(RT)] <- RT
      R_full[colnames(RT), rownames(RT)] <- t(RT)
    }
    
    tab  <- matrix(Inf, length(Lleaves), length(Rleaves),
                   dimnames = list(Lleaves, Rleaves))
    ptrB <- matrix(NA_integer_, length(Lleaves), length(Rleaves))
    ptrC <- matrix(NA_integer_, length(Lleaves), length(Rleaves))
    
    for (ui in seq_along(Lleaves)) {
      u <- Lleaves[ui]
      Lrow <- L_full[as.character(u), ]
      for (wi in seq_along(Rleaves)) {
        w <- Rleaves[wi]
        Rcol <- R_full[, as.character(w)]
        best <- Inf; bestB <- NA_integer_; bestC <- NA_integer_
        for (bi in seq_along(Lleaves)) {
          Lcost <- Lrow[bi]
          if (!is.finite(Lcost)) next
          b <- Lleaves[bi]
          total_vec <- Lcost + Rcol + d[as.character(b), as.character(Rleaves)]
          ci <- which.min(total_vec)
          if (total_vec[ci] < best) { best <- total_vec[ci]; bestB <- b; bestC <- Rleaves[ci] }
        }
        tab[ui, wi] <- best; ptrB[ui, wi] <- bestB; ptrC[ui, wi] <- bestC
      }
    }
    
    res <- list(leaves = c(Lleaves, Rleaves), table = tab,
                Lleaves = Lleaves, Rleaves = Rleaves,
                ptrB = ptrB, ptrC = ptrC, left = left, right = right)
    memo[[node]] <<- res
    res
  }
  
  root <- unique(setdiff(phy$edge[,1], phy$edge[,2]))
  stopifnot(length(root) == 1)
  top <- dp(root)
  
  # backtrack from the globally best (u,w) boundary pair at the root
  get_order <- function(node_res, start, end) {
    if (length(node_res$leaves) == 1) return(node_res$leaves)
    Lleaves <- node_res$Lleaves; Rleaves <- node_res$Rleaves
    if (start %in% Lleaves && end %in% Rleaves) {
      ui <- match(start, Lleaves); wi <- match(end, Rleaves)
      b <- node_res$ptrB[ui, wi]; c_ <- node_res$ptrC[ui, wi]
      c(get_order(node_res$left, start, b), get_order(node_res$right, c_, end))
    } else if (start %in% Rleaves && end %in% Lleaves) {
      rev(get_order(node_res, end, start))
    } else stop("invalid boundary pair")
  }
  
  min_cost <- min(top$table)
  idx <- which(top$table == min_cost, arr.ind = TRUE)[1, ]
  best_order <- get_order(top, top$Lleaves[idx[1]], top$Rleaves[idx[2]])
  
  list(order = as.integer(best_order), cost = min_cost)
}

# ---------------------------------------------------------------
# make_trait_dist(traits, phy): build a correctly ordered/labelled
# distance matrix from a trait table whose rownames are taxon names
# (in ANY row order) matching phy$tip.label.
# ---------------------------------------------------------------
make_trait_dist <- function(traits, phy, method = "euclidean", scale_traits = TRUE) {
  stopifnot(all(phy$tip.label %in% rownames(traits)))
  traits <- traits[phy$tip.label, , drop = FALSE]
  if (scale_traits) traits <- scale(traits)
  d <- as.matrix(dist(traits, method = method))
  rownames(d) <- colnames(d) <- as.character(seq_len(nrow(d)))
  d
}

# ---------------------------------------------------------------
# reorder_tree_olo(phy, d): run the DP and return a topology-
# identical phylo object whose tips are rotated into the OLO order
# (uses ape::rotateConstr, so no plotting package is required).
# ---------------------------------------------------------------
reorder_tree_olo <- function(phy, d) {
  res <- olo_phylo(phy, d)
  ord_labels <- phy$tip.label[res$order]
  phy_ord <- rotateConstr(phy, constraint = ord_labels)
  list(tree = phy_ord, order = ord_labels, cost = res$cost)
}

d <- make_trait_dist(traits, my_tree, method = "euclidean", scale_traits = TRUE)
out <- reorder_tree_olo(my_tree, d)

cat("Minimal total adjacent trait-dissimilarity:", out$cost, "\n")
cat("Optimal tip order:\n"); print(out$order)

# out$tree is a normal phylo object -- plot it as usual:
plot(out$tree, main = "Tree with OLO-optimized tip order")

cat("OLO tree cost",get_layout_cost2(out$tree,traits))

# Experimental reorder algorithm



# =====================================================================
# Utility: children list
# =====================================================================

get_children_list <- function(phy) {
  n_tip  <- length(phy$tip.label)
  n_node <- phy$Nnode
  children <- vector("list", n_tip + n_node)
  
  for (i in seq_len(nrow(phy$edge))) {
    p  <- phy$edge[i, 1]
    ch <- phy$edge[i, 2]
    children[[p]] <- c(children[[p]], ch)
  }
  
  children
}

# =====================================================================
# Distance matrix from traits
# For an actual subquadratic implementation this should be replaced by
# an on-demand geometric distance oracle. For testing O(N^2) ideas,
# explicit construction is acceptable because it is already O(N^2).
# =====================================================================

make_trait_dist <- function(traits, phy, method = "euclidean", scale_traits = TRUE) {
  stopifnot(all(phy$tip.label %in% rownames(traits)))
  
  traits <- traits[phy$tip.label, , drop = FALSE]
  
  if (scale_traits) {
    traits <- scale(traits)
  }
  
  d <- as.matrix(dist(traits, method = method))
  rownames(d) <- colnames(d) <- as.character(seq_len(nrow(d)))
  
  d
}

# =====================================================================
# Small brute-force checker for N <= 10 or 12
# =====================================================================

brute_olo <- function(phy, d) {
  if (!is.matrix(d)) d <- as.matrix(d)
  
  children <- get_children_list(phy)
  n_tip <- length(phy$tip.label)
  
  root <- unique(setdiff(phy$edge[, 1], phy$edge[, 2]))
  stopifnot(length(root) == 1)
  
  gen_orders <- function(node) {
    if (node <= n_tip) {
      return(list(node))
    }
    
    kids <- children[[node]]
    stopifnot(length(kids) == 2)
    
    L <- gen_orders(kids[1])
    R <- gen_orders(kids[2])
    
    out <- list()
    
    for (l in L) {
      for (r in R) {
        out[[length(out) + 1]] <- c(l, r)
        out[[length(out) + 1]] <- c(r, l)
      }
    }
    
    out
  }
  
  orders <- gen_orders(root)
  
  costs <- vapply(
    orders,
    function(o) {
      if (length(o) <= 1) return(0)
      sum(d[cbind(o[-length(o)], o[-1])])
    },
    numeric(1)
  )
  
  ib <- which.min(costs)
  
  list(
    order = orders[[ib]],
    cost  = costs[ib]
  )
}

# =====================================================================
# Detect rooted caterpillar tree
# =====================================================================

is_caterpillar <- function(phy) {
  children <- get_children_list(phy)
  n_tip <- length(phy$tip.label)
  
  root <- unique(setdiff(phy$edge[, 1], phy$edge[, 2]))
  stopifnot(length(root) == 1)
  
  rec <- function(node) {
    if (node <= n_tip) return(TRUE)
    
    kids <- children[[node]]
    if (length(kids) != 2) return(FALSE)
    
    l_leaf <- kids[1] <= n_tip
    r_leaf <- kids[2] <= n_tip
    
    if (l_leaf && r_leaf) return(TRUE)
    if (l_leaf && !r_leaf) return(rec(kids[2]))
    if (!l_leaf && r_leaf) return(rec(kids[1]))
    
    FALSE
  }
  
  rec(root)
}

# =====================================================================
# Exact O(N^2) algorithm for caterpillar trees
# This is a true quadratic fallback for an important tree family.
# =====================================================================

olo_caterpillar_exact <- function(phy, d, verbose = TRUE) {
  if (!is.matrix(d)) d <- as.matrix(d)
  
  children <- get_children_list(phy)
  n_tip <- length(phy$tip.label)
  
  root <- unique(setdiff(phy$edge[, 1], phy$edge[, 2]))
  stopifnot(length(root) == 1)
  
  if (!is_caterpillar(phy)) {
    stop("olo_caterpillar_exact() requires a rooted caterpillar tree.")
  }
  
  # -----------------------------------------------------------------
  # Recover attachment order from base cherry to root.
  # Example: (((a,b),c),d) -> c(a,b,d) = a,b,c,d
  # -----------------------------------------------------------------
  
  attach_order <- function(node) {
    if (node <= n_tip) return(node)
    
    kids <- children[[node]]
    stopifnot(length(kids) == 2)
    
    l_leaf <- kids[1] <= n_tip
    r_leaf <- kids[2] <= n_tip
    
    if (l_leaf && r_leaf) {
      return(kids)
    }
    
    if (l_leaf && !r_leaf) {
      return(c(attach_order(kids[2]), kids[1]))
    }
    
    if (!l_leaf && r_leaf) {
      return(c(attach_order(kids[1]), kids[2]))
    }
    
    stop("Not a caterpillar.")
  }
  
  ord <- attach_order(root)
  N <- length(ord)
  
  if (N == 1) {
    return(list(order = ord, cost = 0))
  }
  
  if (N == 2) {
    return(list(order = ord, cost = d[ord[1], ord[2]]))
  }
  
  # -----------------------------------------------------------------
  # DP:
  # After processing ord[1..k], one endpoint is ord[k],
  # the other endpoint is ord[j], j < k.
  #
  # F[k, j] = minimum cost for that endpoint pair.
  #
  # Transition to q = ord[k+1]:
  #   carry: attach q to current endpoint ord[k]
  #   switch: attach q to other endpoint ord[j]
  # -----------------------------------------------------------------
  
  F <- matrix(Inf, nrow = N, ncol = N)
  prev_switch <- matrix(NA_integer_, nrow = N, ncol = N)
  
  F[2, 1] <- d[ord[2], ord[1]]
  
  for (k in 2:(N - 1)) {
    q <- ord[k + 1]
    
    # carry states: other endpoint remains j
    if (k >= 3) {
      F[k + 1, 1:(k - 1)] <- F[k, 1:(k - 1)] + d[ord[k], q]
    }
    
    # switch state: new other endpoint is ord[k]
    vals <- F[k, 1:(k - 1)] + d[ord[1:(k - 1)], q]
    best_j <- which.min(vals)
    
    F[k + 1, k] <- vals[best_j]
    prev_switch[k + 1, k] <- best_j
  }
  
  j_end <- which.min(F[N, 1:(N - 1)])
  cost <- F[N, j_end]
  
  # -----------------------------------------------------------------
  # Reconstruct one optimal deque order.
  # We backtrack to discover, for every inserted leaf ord[k],
  # which old endpoint it was attached to.
  # -----------------------------------------------------------------
  
  target_endpoint <- numeric(N)
  
  k <- N
  j <- j_end
  
  while (k > 2) {
    if (j == k - 1) {
      # switch transition
      pj <- prev_switch[k, j]
      target_endpoint[k] <- ord[pj]
      j <- pj
    } else {
      # carry transition
      target_endpoint[k] <- ord[k - 1]
      # j unchanged
    }
    
    k <- k - 1
  }
  
  # Replay insertions forward.
  # Initial deque orientation can be chosen arbitrarily; reversal gives same cost.
  deque <- ord[1:2]
  
  for (k in 3:N) {
    tgt <- target_endpoint[k]
    
    if (deque[1] == tgt) {
      deque <- c(ord[k], deque)
    } else if (deque[length(deque)] == tgt) {
      deque <- c(deque, ord[k])
    } else {
      stop("Reconstruction failure in caterpillar DP.")
    }
  }
  
  if (verbose) {
    message("Caterpillar exact O(N^2) DP used.")
    message("N = ", N, ", cost = ", cost)
  }
  
  list(order = as.integer(deque), cost = cost)
}

# =====================================================================
# Term-Envelope profile algorithm for general binary trees
#
# This implements the exact algebraic identity:
#
#   (alpha; p, q) o (beta; r, s)
#     = (alpha + beta + D(q, r); p, s)
#
# and the reversed orientation.
#
# The only pruning used is:
#   1. exact endpoint-pair minimization;
#   2. safe triangle dominance:
#        if alpha_i + D(p_i,p_j) + D(q_i,q_j) <= alpha_j,
#        term i dominates term j for all future metric contexts.
#
# No unsafe heuristic pruning is used.
#
# Complexity is proportional to the total number of surviving terms.
# The hidden unproved theorem is that this total is O(N^2) or better.
# =====================================================================

olo_term_envelope <- function(
    phy,
    d,
    verbose = TRUE,
    prune_interval = 20000,
    max_dominance_terms = 1500,
    hard_term_limit = 1e6,
    tol = 1e-12
) {
  if (!is.matrix(d)) d <- as.matrix(d)
  
  # IMPORTANT: avoid using D, because D is a base R function.
  distm <- d
  
  n_tip <- length(phy$tip.label)
  N <- n_tip
  
  stopifnot(nrow(distm) == N, ncol(distm) == N)
  
  children <- get_children_list(phy)
  memo <- vector("list", n_tip + phy$Nnode)
  
  root <- unique(setdiff(phy$edge[, 1], phy$edge[, 2]))
  stopifnot(length(root) == 1)
  
  empty_terms <- function() {
    data.frame(
      alpha  = numeric(0),
      p      = integer(0),
      q      = integer(0),
      left   = integer(0),
      right  = integer(0),
      orient = integer(0),
      stringsAsFactors = FALSE
    )
  }
  
  # -----------------------------------------------------------------
  # Collapse by endpoint pair: keep best alpha for each (p,q).
  # -----------------------------------------------------------------
  
  collapse_terms <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(empty_terms())
    
    key <- df$p * (N + 1L) + df$q
    
    ord <- order(key, df$alpha)
    df <- df[ord, , drop = FALSE]
    key <- key[ord]
    
    keep <- !duplicated(key)
    df <- df[keep, , drop = FALSE]
    
    rownames(df) <- NULL
    df
  }
  
  # -----------------------------------------------------------------
  # Safe dominance pruning using triangle inequality.
  #
  # Term i dominates term j if:
  #
  #   alpha_i + distm(p_i,p_j) + distm(q_i,q_j) <= alpha_j
  #
  # Then for every external u,v:
  #
  #   alpha_i + distm(u,p_i) + distm(v,q_i)
  #     <= alpha_j + distm(u,p_j) + distm(v,q_j)
  #
  # This is safe for any metric.
  # -----------------------------------------------------------------
  
  dominance_prune <- function(df) {
    m <- nrow(df)
    if (m <= 1) return(df)
    
    # Skip expensive quadratic dominance checks when the term set is large.
    # Skipping safe pruning does not affect exactness; it only affects compression.
    if (m > max_dominance_terms) {
      return(df)
    }
    
    keep <- rep(TRUE, m)
    
    for (i in seq_len(m)) {
      if (!keep[i]) next
      
      vals <- df$alpha +
        distm[df$p, df$p[i]] +
        distm[df$q, df$q[i]]
      
      vals[i] <- Inf
      
      if (any(vals <= df$alpha[i] - tol)) {
        keep[i] <- FALSE
      }
    }
    
    df <- df[keep, , drop = FALSE]
    rownames(df) <- NULL
    df
  }
  
  # -----------------------------------------------------------------
  # Add candidate terms and periodically prune.
  # -----------------------------------------------------------------
  
  add_candidates <- function(acc, cand) {
    if (nrow(cand) == 0) return(acc)
    
    acc <- rbind(acc, cand)
    
    if (nrow(acc) > prune_interval) {
      acc <- collapse_terms(acc)
      acc <- dominance_prune(acc)
    }
    
    if (nrow(acc) > hard_term_limit) {
      stop(
        "Term envelope exceeded hard_term_limit = ", hard_term_limit,
        ". Profile compression failed or hard_term_limit is too low."
      )
    }
    
    acc
  }
  
  # -----------------------------------------------------------------
  # Compose two child profiles.
  # -----------------------------------------------------------------
  
  compose_terms <- function(L, R) {
    acc <- empty_terms()
    
    nL <- nrow(L)
    nR <- nrow(R)
    
    if (nL == 0 || nR == 0) return(empty_terms())
    
    # Orientation 1: L then R
    #
    # L term: (alpha_L; p_L, q_L)
    # R term: (alpha_R; p_R, q_R)
    # seam cost distm(q_L, p_R)
    # new term: (alpha_L + alpha_R + distm(q_L,p_R); p_L, q_R)
    
    for (i in seq_len(nL)) {
      lp <- L$p[i]
      lq <- L$q[i]
      la <- L$alpha[i]
      
      cand <- data.frame(
        alpha  = la + R$alpha + distm[lq, R$p],
        p      = rep(as.integer(lp), nR),
        q      = R$q,
        left   = rep(as.integer(i), nR),
        right  = seq_len(nR),
        orient = 1L,
        stringsAsFactors = FALSE
      )
      
      acc <- add_candidates(acc, cand)
    }
    
    # Orientation 2: R then L
    #
    # R term: (alpha_R; p_R, q_R)
    # L term: (alpha_L; p_L, q_L)
    # seam cost distm(q_R, p_L)
    # new term: (alpha_R + alpha_L + distm(q_R,p_L); p_R, q_L)
    
    for (j in seq_len(nR)) {
      rp <- R$p[j]
      rq <- R$q[j]
      ra <- R$alpha[j]
      
      cand <- data.frame(
        alpha  = ra + L$alpha + distm[rq, L$p],
        p      = rep(as.integer(rp), nL),
        q      = L$q,
        left   = seq_len(nL),
        right  = rep(as.integer(j), nL),
        orient = 2L,
        stringsAsFactors = FALSE
      )
      
      acc <- add_candidates(acc, cand)
    }
    
    acc <- collapse_terms(acc)
    acc <- dominance_prune(acc)
    
    acc
  }
  
  # -----------------------------------------------------------------
  # Bottom-up DP.
  # -----------------------------------------------------------------
  
  dp <- function(node) {
    if (!is.null(memo[[node]])) {
      return(memo[[node]])
    }
    
    if (node <= n_tip) {
      terms <- data.frame(
        alpha  = 0,
        p      = as.integer(node),
        q      = as.integer(node),
        left   = NA_integer_,
        right  = NA_integer_,
        orient = 0L,
        stringsAsFactors = FALSE
      )
      
      res <- list(
        is_leaf = TRUE,
        leaf    = node,
        terms   = terms
      )
      
      memo[[node]] <<- res
      return(res)
    }
    
    kids <- children[[node]]
    if (length(kids) != 2) {
      stop("Tree must be fully bifurcating.")
    }
    
    left_res  <- dp(kids[1])
    right_res <- dp(kids[2])
    
    terms <- compose_terms(left_res$terms, right_res$terms)
    
    if (verbose) {
      message(
        "Node ", node,
        ": left terms = ", nrow(left_res$terms),
        ", right terms = ", nrow(right_res$terms),
        ", surviving terms = ", nrow(terms)
      )
    }
    
    res <- list(
      is_leaf = FALSE,
      left    = kids[1],
      right   = kids[2],
      terms   = terms
    )
    
    memo[[node]] <<- res
    res
  }
  
  top <- dp(root)
  
  if (nrow(top$terms) == 0) {
    stop("Term envelope became empty unexpectedly.")
  }
  
  best_idx <- which.min(top$terms$alpha)
  best_cost <- top$terms$alpha[best_idx]
  
  # -----------------------------------------------------------------
  # Reconstruct order from selected root term.
  # -----------------------------------------------------------------
  
  get_order <- function(node, term_idx) {
    res <- memo[[node]]
    
    if (res$is_leaf) {
      return(res$leaf)
    }
    
    term <- res$terms[term_idx, ]
    
    if (term$orient == 1L) {
      left_order  <- get_order(res$left,  term$left)
      right_order <- get_order(res$right, term$right)
      return(c(left_order, right_order))
    } else if (term$orient == 2L) {
      right_order <- get_order(res$right, term$right)
      left_order  <- get_order(res$left,  term$left)
      return(c(right_order, left_order))
    } else {
      stop("Invalid orientation in term reconstruction.")
    }
  }
  
  best_order <- get_order(root, best_idx)
  
  term_counts <- vapply(
    memo,
    function(z) {
      if (is.null(z) || is.null(z$terms)) {
        0L
      } else {
        nrow(z$terms)
      }
    },
    integer(1)
  )
  
  list(
    order       = as.integer(best_order),
    cost        = best_cost,
    term_counts = term_counts
  )
}

# =====================================================================
# Wrapper: use exact O(N^2) caterpillar algorithm when possible,
# otherwise use the Term-Envelope candidate.
# =====================================================================

reorder_tree_olo_fast <- function(
    phy,
    d,
    verbose = TRUE,
    prune_interval = 20000,
    max_dominance_terms = 1500,
    hard_term_limit = 1e6
) {
  if (is_caterpillar(phy)) {
    res <- olo_caterpillar_exact(phy, d, verbose = verbose)
  } else {
    res <- olo_term_envelope(
      phy,
      d,
      verbose = verbose,
      prune_interval = prune_interval,
      max_dominance_terms = max_dominance_terms,
      hard_term_limit = hard_term_limit
    )
  }
  
  ord_labels <- phy$tip.label[res$order]
  
  phy_ord <- tryCatch(
    rotateConstr(phy, constraint = ord_labels),
    error = function(e) {
      warning("rotateConstr failed; returning original tree.")
      phy
    }
  )
  
  list(
    tree  = phy_ord,
    order = ord_labels,
    cost  = res$cost
  )
}

# =====================================================================
# Run
# =====================================================================
test_sizes <- c(50, 100, 150, 200, 250, 300)

for (n in test_sizes) {
  cat("\n\n==================== N =", n, "====================\n")
  
  set.seed(42)
  
  my_tree <- rtree(
    n,
    tip.label = paste0("Tip", 1:n)
  )
  
  traits <- matrix(
    runif(n * 10, 0, 5),
    nrow = n,
    dimnames = list(
      my_tree$tip.label,
      paste0("Trait", 1:10)
    )
  )
  
  d <- make_trait_dist(
    traits,
    my_tree,
    method = "euclidean",
    scale_traits = TRUE
  )
  
  t0 <- Sys.time()
  
  out <- reorder_tree_olo_fast(
    my_tree,
    d,
    verbose = FALSE,
    prune_interval = 50000,
    max_dominance_terms = 1500,
    hard_term_limit = 5e6
  )
  
  t1 <- Sys.time()
  
  layout_cost <- get_layout_cost2(out$tree, traits, scale_traits = TRUE)
  
  cat("N:", n, "\n")
  cat("Term-envelope cost:", out$cost, "\n")
  cat("Layout cost:       ", layout_cost, "\n")
  cat("Time:              ", format(t1 - t0), "\n")
}

d <- make_trait_dist(
  traits,
  my_tree,
  method = "euclidean",
  scale_traits = TRUE
)

out <- reorder_tree_olo_fast(
  my_tree,
  d,
  verbose = TRUE,
  prune_interval = 20000,
  max_dominance_terms = 1500,
  hard_term_limit = 1e6
)

cat("Minimal total adjacent trait-dissimilarity:", out$cost, "\n")
cat("Optimal tip order:\n")
print(out$order)
cat("experiment tree cost",get_layout_cost2(out$tree,traits))
# =====================================================================
# Optional brute-force verification for small trees
# =====================================================================

if (nTips <= 10) {
  bf <- brute_olo(my_tree, d)
  cat("Brute-force cost:", bf$cost, "\n")
  cat("Fast cost:      ", out$cost, "\n")
  stopifnot(abs(bf$cost - out$cost) < 1e-10)
}