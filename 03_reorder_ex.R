# =====================================================================
# Full diagnostic script:
# Original OLO vs experimental Term-Envelope OLO
# =====================================================================

library(ape)

# =====================================================================
# Basic utilities
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

get_layout_cost2 <- function(tree, traits, scale_traits = TRUE) {
  n <- length(tree$tip.label)
  
  tree <- reorder(tree, order = "cladewise")
  
  is_tip <- tree$edge[, 2] <= n
  ord_nodes <- tree$edge[is_tip, 2]
  ord_labels <- tree$tip.label[ord_nodes]
  
  X_all <- traits[tree$tip.label, , drop = FALSE]
  
  if (scale_traits) {
    X_all <- scale(X_all)
  }
  
  X <- X_all[ord_labels, , drop = FALSE]
  
  sum(sqrt(rowSums(diff(X)^2)))
}

cost_from_order_labels <- function(order_labels, phy, d) {
  ids <- match(order_labels, phy$tip.label)
  
  if (length(ids) <= 1) {
    return(0)
  }
  
  sum(d[cbind(ids[-length(ids)], ids[-1])])
}

# =====================================================================
# Brute force for very small trees
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
# Original OLO implementation
# =====================================================================

olo_phylo <- function(phy, d) {
  
  n_tip <- length(phy$tip.label)
  if (!is.matrix(d)) d <- as.matrix(d)
  stopifnot(nrow(d) == n_tip, ncol(d) == n_tip)
  
  children <- get_children_list(phy)
  memo <- vector("list", n_tip + phy$Nnode)
  
  dp <- function(node) {
    if (!is.null(memo[[node]])) return(memo[[node]])
    kids <- children[[node]]
    
    if (is.null(kids)) {
      res <- list(
        leaves = node,
        table = matrix(0, 1, 1, dimnames = list(node, node))
      )
      memo[[node]] <<- res
      return(res)
    }
    
    if (length(kids) != 2) {
      stop(
        "Node ", node, " has ", length(kids), " children; the tree ",
        "must be fully bifurcating. Resolve polytomies first with ",
        "ape::multi2di()."
      )
    }
    
    left  <- dp(kids[1])
    right <- dp(kids[2])
    
    Lleaves <- left$leaves
    Rleaves <- right$leaves
    
    LT <- left$table
    RT <- right$table
    
    L_full <- matrix(
      Inf,
      length(Lleaves),
      length(Lleaves),
      dimnames = list(Lleaves, Lleaves)
    )
    
    if (length(Lleaves) == 1) {
      L_full[1, 1] <- 0
    } else {
      L_full[rownames(LT), colnames(LT)] <- LT
      L_full[colnames(LT), rownames(LT)] <- t(LT)
    }
    
    R_full <- matrix(
      Inf,
      length(Rleaves),
      length(Rleaves),
      dimnames = list(Rleaves, Rleaves)
    )
    
    if (length(Rleaves) == 1) {
      R_full[1, 1] <- 0
    } else {
      R_full[rownames(RT), colnames(RT)] <- RT
      R_full[colnames(RT), rownames(RT)] <- t(RT)
    }
    
    tab  <- matrix(
      Inf,
      length(Lleaves),
      length(Rleaves),
      dimnames = list(Lleaves, Rleaves)
    )
    
    ptrB <- matrix(NA_integer_, length(Lleaves), length(Rleaves))
    ptrC <- matrix(NA_integer_, length(Lleaves), length(Rleaves))
    
    for (ui in seq_along(Lleaves)) {
      u <- Lleaves[ui]
      Lrow <- L_full[as.character(u), ]
      
      for (wi in seq_along(Rleaves)) {
        w <- Rleaves[wi]
        Rcol <- R_full[, as.character(w)]
        
        best <- Inf
        bestB <- NA_integer_
        bestC <- NA_integer_
        
        for (bi in seq_along(Lleaves)) {
          Lcost <- Lrow[bi]
          if (!is.finite(Lcost)) next
          
          b <- Lleaves[bi]
          
          total_vec <- Lcost + Rcol + d[as.character(b), as.character(Rleaves)]
          ci <- which.min(total_vec)
          
          if (total_vec[ci] < best) {
            best <- total_vec[ci]
            bestB <- b
            bestC <- Rleaves[ci]
          }
        }
        
        tab[ui, wi] <- best
        ptrB[ui, wi] <- bestB
        ptrC[ui, wi] <- bestC
      }
    }
    
    res <- list(
      leaves = c(Lleaves, Rleaves),
      table = tab,
      Lleaves = Lleaves,
      Rleaves = Rleaves,
      ptrB = ptrB,
      ptrC = ptrC,
      left = left,
      right = right
    )
    
    memo[[node]] <<- res
    res
  }
  
  root <- unique(setdiff(phy$edge[, 1], phy$edge[, 2]))
  stopifnot(length(root) == 1)
  
  top <- dp(root)
  
  get_order <- function(node_res, start, end) {
    if (length(node_res$leaves) == 1) {
      return(node_res$leaves)
    }
    
    Lleaves <- node_res$Lleaves
    Rleaves <- node_res$Rleaves
    
    if (start %in% Lleaves && end %in% Rleaves) {
      ui <- match(start, Lleaves)
      wi <- match(end, Rleaves)
      
      b <- node_res$ptrB[ui, wi]
      c_ <- node_res$ptrC[ui, wi]
      
      c(
        get_order(node_res$left, start, b),
        get_order(node_res$right, c_, end)
      )
    } else if (start %in% Rleaves && end %in% Lleaves) {
      rev(get_order(node_res, end, start))
    } else {
      stop("invalid boundary pair")
    }
  }
  
  min_cost <- min(top$table)
  
  idx <- which(top$table == min_cost, arr.ind = TRUE)[1, ]
  
  best_order <- get_order(
    top,
    top$Lleaves[idx[1]],
    top$Rleaves[idx[2]]
  )
  
  list(
    order = as.integer(best_order),
    cost = min_cost
  )
}

# =====================================================================
# Caterpillar detection
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
# Exact O(N^2) caterpillar algorithm
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
    return(
      list(
        order = ord,
        cost = 0,
        term_counts = NA_integer_,
        stats = list(
          method = "caterpillar",
          candidates_generated = NA_real_,
          max_raw_terms = NA_real_,
          max_surviving_terms = NA_real_,
          total_surviving_terms = NA_real_,
          dominance_skips = NA_real_
        )
      )
    )
  }
  
  if (N == 2) {
    return(
      list(
        order = ord,
        cost = d[ord[1], ord[2]],
        term_counts = NA_integer_,
        stats = list(
          method = "caterpillar",
          candidates_generated = NA_real_,
          max_raw_terms = NA_real_,
          max_surviving_terms = NA_real_,
          total_surviving_terms = NA_real_,
          dominance_skips = NA_real_
        )
      )
    )
  }
  
  F <- matrix(Inf, nrow = N, ncol = N)
  prev_switch <- matrix(NA_integer_, nrow = N, ncol = N)
  
  F[2, 1] <- d[ord[2], ord[1]]
  
  for (k in 2:(N - 1)) {
    q <- ord[k + 1]
    
    if (k >= 3) {
      F[k + 1, 1:(k - 1)] <- F[k, 1:(k - 1)] + d[ord[k], q]
    }
    
    vals <- F[k, 1:(k - 1)] + d[ord[1:(k - 1)], q]
    best_j <- which.min(vals)
    
    F[k + 1, k] <- vals[best_j]
    prev_switch[k + 1, k] <- best_j
  }
  
  j_end <- which.min(F[N, 1:(N - 1)])
  cost <- F[N, j_end]
  
  target_endpoint <- numeric(N)
  
  k <- N
  j <- j_end
  
  while (k > 2) {
    if (j == k - 1) {
      pj <- prev_switch[k, j]
      target_endpoint[k] <- ord[pj]
      j <- pj
    } else {
      target_endpoint[k] <- ord[k - 1]
    }
    
    k <- k - 1
  }
  
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
  
  list(
    order = as.integer(deque),
    cost = cost,
    term_counts = NA_integer_,
    stats = list(
      method = "caterpillar",
      candidates_generated = NA_real_,
      max_raw_terms = NA_real_,
      max_surviving_terms = NA_real_,
      total_surviving_terms = NA_real_,
      dominance_skips = NA_real_
    )
  )
}

# =====================================================================
# Experimental Term-Envelope implementation with diagnostics
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
  
  distm <- d
  
  n_tip <- length(phy$tip.label)
  N <- n_tip
  
  stopifnot(nrow(distm) == N, ncol(distm) == N)
  
  children <- get_children_list(phy)
  memo <- vector("list", n_tip + phy$Nnode)
  
  root <- unique(setdiff(phy$edge[, 1], phy$edge[, 2]))
  stopifnot(length(root) == 1)
  
  # -----------------------------------------------------------------
  # Diagnostic statistics
  # -----------------------------------------------------------------
  
  stats <- new.env()
  stats$candidates_generated <- 0
  stats$max_raw_terms <- 0
  stats$max_surviving_terms <- 0
  stats$dominance_skips <- 0
  
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
  
  collapse_terms <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(empty_terms())
    
    key <- as.numeric(df$p) * (N + 1) + as.numeric(df$q)
    
    ord <- order(key, df$alpha)
    df <- df[ord, , drop = FALSE]
    key <- key[ord]
    
    keep <- !duplicated(key)
    df <- df[keep, , drop = FALSE]
    
    rownames(df) <- NULL
    df
  }
  
  dominance_prune <- function(df) {
    m <- nrow(df)
    if (m <= 1) return(df)
    
    if (m > max_dominance_terms) {
      stats$dominance_skips <- stats$dominance_skips + 1
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
  
  add_candidates <- function(acc, cand) {
    if (nrow(cand) == 0) return(acc)
    
    stats$candidates_generated <- stats$candidates_generated + nrow(cand)
    
    acc <- rbind(acc, cand)
    
    if (nrow(acc) > stats$max_raw_terms) {
      stats$max_raw_terms <- nrow(acc)
    }
    
    needs_prune <- nrow(acc) > prune_interval
    
    if (
      !is.null(hard_term_limit) &&
      is.finite(hard_term_limit) &&
      nrow(acc) > hard_term_limit
    ) {
      needs_prune <- TRUE
    }
    
    if (needs_prune) {
      acc <- collapse_terms(acc)
      acc <- dominance_prune(acc)
    }
    
    if (nrow(acc) > stats$max_surviving_terms) {
      stats$max_surviving_terms <- nrow(acc)
    }
    
    if (
      !is.null(hard_term_limit) &&
      is.finite(hard_term_limit) &&
      nrow(acc) > hard_term_limit
    ) {
      stop(
        "Term envelope exceeded hard_term_limit = ", hard_term_limit,
        ". Profile compression failed or hard_term_limit is too low."
      )
    }
    
    acc
  }
  
  compose_terms <- function(L, R) {
    acc <- empty_terms()
    
    nL <- nrow(L)
    nR <- nrow(R)
    
    if (nL == 0 || nR == 0) return(empty_terms())
    
    # Orientation 1: L then R
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
    
    if (nrow(acc) > stats$max_surviving_terms) {
      stats$max_surviving_terms <- nrow(acc)
    }
    
    acc
  }
  
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
    
    if (nrow(terms) > stats$max_surviving_terms) {
      stats$max_surviving_terms <- nrow(terms)
    }
    
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
  
  stats$total_surviving_terms <- sum(term_counts)
  
  list(
    order = as.integer(best_order),
    cost  = best_cost,
    term_counts = term_counts,
    stats = list(
      method = "term-envelope",
      candidates_generated = stats$candidates_generated,
      max_raw_terms = stats$max_raw_terms,
      max_surviving_terms = stats$max_surviving_terms,
      total_surviving_terms = stats$total_surviving_terms,
      dominance_skips = stats$dominance_skips
    )
  )
}

# =====================================================================
# Fast core wrapper
# =====================================================================

olo_fast_core <- function(phy, d, verbose = FALSE, ...) {
  if (is_caterpillar(phy)) {
    res <- olo_caterpillar_exact(phy, d, verbose = verbose)
    res$method <- "caterpillar-exact"
    return(res)
  }
  
  res <- olo_term_envelope(phy, d, verbose = verbose, ...)
  res$method <- "term-envelope"
  res
}

# =====================================================================
# Optional full reordered tree wrapper for the experimental algorithm
# =====================================================================

reorder_tree_olo_fast <- function(
    phy,
    d,
    verbose = TRUE,
    prune_interval = 20000,
    max_dominance_terms = 1500,
    hard_term_limit = 1e6
) {
  res <- olo_fast_core(
    phy,
    d,
    verbose = verbose,
    prune_interval = prune_interval,
    max_dominance_terms = max_dominance_terms,
    hard_term_limit = hard_term_limit
  )
  
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
    cost  = res$cost,
    stats = res$stats
  )
}

# =====================================================================
# Diagnostic runner for one N
# =====================================================================

diagnose_one <- function(
    nTips,
    seed = 42,
    run_original = TRUE,
    run_fast = TRUE,
    verify_brute = TRUE,
    verbose_fast = FALSE,
    hard_term_limit = 2e6
) {
  set.seed(seed)
  
  my_tree <- rtree(
    nTips,
    tip.label = paste0("Tip", seq_len(nTips))
  )
  
  traits <- matrix(
    runif(nTips * 10, 0, 5),
    nrow = nTips,
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
  
  brute_cost <- NA_real_
  brute_ok <- NA
  
  if (verify_brute && nTips <= 10) {
    bf <- brute_olo(my_tree, d)
    brute_cost <- bf$cost
  }
  
  # -----------------------------------------------------------------
  # Original OLO
  # -----------------------------------------------------------------
  
  orig_time_sec <- NA_real_
  orig_cost <- NA_real_
  orig_direct_cost <- NA_real_
  orig_labels <- NULL
  orig_error <- NA_character_
  
  if (run_original) {
    orig_try <- tryCatch(
      {
        t_orig <- system.time(
          {
            orig_core <- olo_phylo(my_tree, d)
          }
        )
        
        list(core = orig_core, time = t_orig)
      },
      error = function(e) e
    )
    
    if (inherits(orig_try, "error")) {
      orig_error <- conditionMessage(orig_try)
    } else {
      orig_core <- orig_try$core
      orig_time_sec <- as.numeric(orig_try$time["elapsed"])
      
      orig_labels <- my_tree$tip.label[orig_core$order]
      orig_cost <- orig_core$cost
      orig_direct_cost <- cost_from_order_labels(orig_labels, my_tree, d)
    }
  }
  
  # -----------------------------------------------------------------
  # Experimental fast OLO
  # -----------------------------------------------------------------
  
  fast_time_sec <- NA_real_
  fast_cost <- NA_real_
  fast_direct_cost <- NA_real_
  fast_labels <- NULL
  fast_error <- NA_character_
  
  method_fast <- NA_character_
  candidates_generated <- NA_real_
  max_raw_terms <- NA_real_
  max_surviving_terms <- NA_real_
  total_surviving_terms <- NA_real_
  dominance_skips <- NA_real_
  
  if (run_fast) {
    fast_try <- tryCatch(
      {
        t_fast <- system.time(
          {
            fast_core <- olo_fast_core(
              my_tree,
              d,
              verbose = verbose_fast,
              hard_term_limit = hard_term_limit
            )
          }
        )
        
        list(core = fast_core, time = t_fast)
      },
      error = function(e) e
    )
    
    if (inherits(fast_try, "error")) {
      fast_error <- conditionMessage(fast_try)
    } else {
      fast_core <- fast_try$core
      fast_time_sec <- as.numeric(fast_try$time["elapsed"])
      
      method_fast <- fast_core$method
      
      fast_labels <- my_tree$tip.label[fast_core$order]
      fast_cost <- fast_core$cost
      fast_direct_cost <- cost_from_order_labels(fast_labels, my_tree, d)
      
      if (!is.null(fast_core$stats)) {
        candidates_generated <- fast_core$stats$candidates_generated
        max_raw_terms <- fast_core$stats$max_raw_terms
        max_surviving_terms <- fast_core$stats$max_surviving_terms
        total_surviving_terms <- fast_core$stats$total_surviving_terms
        dominance_skips <- fast_core$stats$dominance_skips
      }
    }
  }
  
  # -----------------------------------------------------------------
  # Comparisons
  # -----------------------------------------------------------------
  
  cost_diff <- NA_real_
  costs_equal <- NA
  
  if (!is.na(orig_cost) && !is.na(fast_cost)) {
    cost_diff <- fast_cost - orig_cost
    costs_equal <- isTRUE(abs(cost_diff) < 1e-8)
  }
  
  orders_identical <- NA
  
  if (!is.null(orig_labels) && !is.null(fast_labels)) {
    orders_identical <- identical(orig_labels, fast_labels)
  }
  
  if (verify_brute && nTips <= 10) {
    brute_ok <- TRUE
    
    if (!is.na(orig_cost)) {
      brute_ok <- brute_ok && isTRUE(abs(orig_cost - brute_cost) < 1e-8)
    }
    
    if (!is.na(fast_cost)) {
      brute_ok <- brute_ok && isTRUE(abs(fast_cost - brute_cost) < 1e-8)
    }
  }
  
  data.frame(
    N = nTips,
    method_fast = method_fast,
    
    original_time_sec = orig_time_sec,
    fast_time_sec = fast_time_sec,
    
    speedup_original_over_fast = ifelse(
      !is.na(orig_time_sec) && !is.na(fast_time_sec) && fast_time_sec > 0,
      orig_time_sec / fast_time_sec,
      NA_real_
    ),
    
    original_cost = orig_cost,
    fast_cost = fast_cost,
    cost_diff_fast_minus_orig = cost_diff,
    costs_equal = costs_equal,
    orders_identical = orders_identical,
    
    original_direct_cost = orig_direct_cost,
    fast_direct_cost = fast_direct_cost,
    
    brute_cost = brute_cost,
    brute_ok = brute_ok,
    
    candidates_generated = candidates_generated,
    max_raw_terms = max_raw_terms,
    max_surviving_terms = max_surviving_terms,
    total_surviving_terms = total_surviving_terms,
    dominance_skips = dominance_skips,
    
    original_error = ifelse(is.na(orig_error), "", orig_error),
    fast_error = ifelse(is.na(fast_error), "", fast_error),
    
    stringsAsFactors = FALSE
  )
}

# =====================================================================
# Run diagnostic sweep
# =====================================================================

# Quick sizes. Increase carefully.
QUICK_SIZES <- c(8, 12, 16, 20)

# Optional larger sizes. Original OLO may become very slow.
RUN_LARGE <- FALSE
LARGE_SIZES <- c(30, 40, 50)

# Optional N=500 test.
# Set RUN_N500 to TRUE if you want to run it.
# Set RUN_N500_ORIGINAL to TRUE only if you are willing to wait for original OLO.
RUN_N500 <- FALSE
RUN_N500_ORIGINAL <- FALSE

results <- list()

for (n in QUICK_SIZES) {
  cat("\n\n==================== Diagnosing N =", n, "====================\n")
  
  results[[length(results) + 1]] <- diagnose_one(
    nTips = n,
    seed = 42,
    run_original = TRUE,
    run_fast = TRUE,
    verify_brute = (n <= 10),
    verbose_fast = FALSE,
    hard_term_limit = 2e6
  )
  
  print(tail(results, 1)[[1]])
}

if (RUN_LARGE) {
  for (n in LARGE_SIZES) {
    cat("\n\n==================== Diagnosing large N =", n, "====================\n")
    
    results[[length(results) + 1]] <- diagnose_one(
      nTips = n,
      seed = 42,
      run_original = TRUE,
      run_fast = TRUE,
      verify_brute = FALSE,
      verbose_fast = FALSE,
      hard_term_limit = 5e6
    )
    
    print(tail(results, 1)[[1]])
  }
}

if (RUN_N500) {
  cat("\n\n==================== Diagnosing N = 500 ====================\n")
  
  results[[length(results) + 1]] <- diagnose_one(
    nTips = 500,
    seed = 42,
    run_original = RUN_N500_ORIGINAL,
    run_fast = TRUE,
    verify_brute = FALSE,
    verbose_fast = FALSE,
    hard_term_limit = 10e6
  )
  
  print(tail(results, 1)[[1]])
}

results_df <- do.call(rbind, results)
rownames(results_df) <- NULL

cat("\n\n==================== FINAL DIAGNOSTIC TABLE ====================\n")
print(results_df)