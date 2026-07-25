library(ape)
library(reticulate)

set.seed(42)

# -----------------------------------------------------------------------------
# 1. Setup Python & Hugging Face Model via reticulate
# -----------------------------------------------------------------------------

# IMPORTANT: reticulate::py_install() installs into whatever Python
# environment reticulate is currently configured to use. If you get
# "ModuleNotFoundError: No module named 'transformers'" even after
# installing it manually, it almost always means reticulate is pointed at
# a *different* Python than the one you installed into. Check with:
#   reticulate::py_config()
#
# Pin reticulate to the existing project venv BEFORE anything else runs.
# Use the venv's python.exe directly (most reliable on Windows) - note the
# path below is the venv root, not the \Scripts subfolder.
reticulate::use_python(
  "E:/Dars/Phylogeny Comparative methods/Tree reordering/.venv/Scripts/python.exe",
  required = TRUE
)

required_pkgs <- c("transformers", "torch", "accelerate")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, py_module_available, logical(1))]
if (length(missing_pkgs) > 0) {
  message("Installing missing Python packages: ", paste(missing_pkgs, collapse = ", "))
  py_install(missing_pkgs, pip = TRUE)
}

transformers <- import("transformers")
torch <- import("torch")

cat("Loading Hugging Face model for tree optimization...\n")

model_id  <- "Qwen/Qwen2.5-Coder-0.5B-Instruct"# Fast, lightweight instruction model
tokenizer <- transformers$AutoTokenizer$from_pretrained(model_id)

# float16 generally requires a GPU; fall back to float32 on CPU-only machines
use_cuda <- torch$cuda$is_available()
dtype    <- if (use_cuda) torch$float16 else torch$float32

model <- transformers$AutoModelForCausalLM$from_pretrained(
  model_id,
  torch_dtype = dtype,
  device_map = "auto"
)

# Qwen models don't define a pad token by default; set one so generate()
# doesn't warn or error on batches
if (is.null(tokenizer$pad_token)) {
  tokenizer$pad_token <- tokenizer$eos_token
}

# Do the "keep only the newly generated tokens" step in native Python.
# Slicing a torch tensor with R's `[` mixes up 1-based vs 0-based indexing
# and doesn't build a real Python slice object, so it's done here instead.
py_run_string("
def extract_new_tokens(outputs, input_ids):
    return outputs[0, input_ids.shape[1]:]
")

# Function to query LLM for strategic decision
query_llm_optimizer <- function(system_prompt, user_prompt) {
  messages <- list(
    list("role" = "system", "content" = system_prompt),
    list("role" = "user",   "content" = user_prompt)
  )
  prompt_text <- tokenizer$apply_chat_template(
    messages,
    tokenize = FALSE,
    add_generation_prompt = TRUE
  )
  inputs <- tokenizer(prompt_text, return_tensors = "pt")$to(model$device)
  
  outputs <- model$generate(
    inputs$input_ids,
    max_new_tokens = as.integer(120), # Increased token space for chain-of-thought reasoning
    temperature    = 0.4,
    do_sample      = TRUE,
    pad_token_id   = tokenizer$pad_token_id
  )
  
  gen_tokens <- py$extract_new_tokens(outputs, inputs$input_ids)
  response   <- tokenizer$decode(gen_tokens, skip_special_tokens = TRUE)
  trimws(response)
}


# -----------------------------------------------------------------------------
# 2. Synthetic Tree & Trait Data Setup
# -----------------------------------------------------------------------------
my_tree <- rtree(100, tip.label = paste0("Tip", 1:100))

traits <- matrix(
  runif(100 * 6, 0, 5),
  nrow = 100,
  dimnames = list(my_tree$tip.label, paste0("Trait", 1:6))
)

# -----------------------------------------------------------------------------
# 3. Helper Functions & Sanitization
# -----------------------------------------------------------------------------
sanitize_tree <- function(tree) {
  ape::read.tree(text = ape::write.tree(tree))
}

get_plot_order <- function(tree) {
  n <- length(tree$tip.label)
  idx <- tree$edge[tree$edge[, 2] <= n, 2]
  tree$tip.label[idx]
}

trait_cost <- function(tree, traits) {
  ord <- get_plot_order(tree)
  X <- traits[ord, , drop = FALSE]
  sum(sqrt(rowSums((X[-1, , drop = FALSE] - X[-nrow(X), , drop = FALSE])^2)))
}

get_clade_sizes <- function(tree) {
  Ntip <- length(tree$tip.label)
  Nnode <- tree$Nnode
  sizes <- integer(Ntip + Nnode)
  sizes[1:Ntip] <- 1L
  pp <- ape::prop.part(tree)
  for (i in seq_along(pp)) {
    sizes[Ntip + i] <- length(pp[[i]])
  }
  return(sizes)
}

# TABU-LOCKED Greedy Local Descent:
# Skips tabu nodes to prevent instantly undoing macro strategy perturbations!
greedy_local_descent <- function(tree, traits, max_passes = 5, tabu_list = integer(0)) {
  current <- sanitize_tree(tree)
  current_cost <- trait_cost(current, traits)
  
  Ntip <- length(current$tip.label)
  nodes <- (Ntip + 1):(Ntip + current$Nnode)
  # Exclude nodes locked by Tabu
  nodes <- nodes[!nodes %in% tabu_list]
  
  improved <- TRUE
  pass <- 0
  
  while (improved && pass < max_passes) {
    improved <- FALSE
    pass <- pass + 1
    
    for (nd in nodes) {
      cand <- ape::rotate(current, node = as.integer(nd))
      cand_cost <- trait_cost(cand, traits)
      
      if (cand_cost < current_cost - 1e-4) {
        current <- cand
        current_cost <- cand_cost
        improved <- TRUE
      }
    }
    current <- sanitize_tree(current)
  }
  
  return(list(tree = current, cost = current_cost))
}

find_tree_bottlenecks <- function(tree, traits, top_k = 6) {
  ord <- get_plot_order(tree)
  X <- traits[ord, , drop = FALSE]
  dists <- sqrt(rowSums((X[-1, , drop = FALSE] - X[-nrow(X), , drop = FALSE])^2))
  worst_idx <- order(dists, decreasing = TRUE)[1:top_k]
  
  mrca_matrix <- ape::mrca(tree)
  clade_sizes <- get_clade_sizes(tree)
  
  bottlenecks <- list()
  for (i in seq_along(worst_idx)) {
    idx <- worst_idx[i]
    tip1 <- ord[idx]
    tip2 <- ord[idx + 1]
    mrca_node <- as.integer(mrca_matrix[tip1, tip2])
    
    bottlenecks[[i]] <- list(
      dist = round(dists[idx], 3),
      mrca_node = mrca_node,
      clade_size = clade_sizes[mrca_node]
    )
  }
  return(bottlenecks)
}

# -----------------------------------------------------------------------------
# 4. Strategy Execution Engine (R handles the physics)
# -----------------------------------------------------------------------------
# Trait-Centroid Alignment Operator
trait_centroid_align <- function(tree, traits, target_node) {
  # Rotates target_node if flipping aligns clade trait means closer to neighbors
  cand <- ape::rotate(tree, node = as.integer(target_node))
  if (trait_cost(cand, traits) < trait_cost(tree, traits)) {
    return(cand)
  }
  return(tree)
}

apply_llm_strategy <- function(tree, traits, strategy_choice, bottlenecks, tabu_list, temp = 1.0) {
  Ntip <- length(tree$tip.label)
  clade_sizes <- get_clade_sizes(tree)
  all_internal <- (Ntip + 1):(Ntip + tree$Nnode)
  
  valid_bottlenecks <- sapply(bottlenecks, function(b) b$mrca_node)
  valid_bottlenecks <- valid_bottlenecks[!valid_bottlenecks %in% tabu_list]
  if (length(valid_bottlenecks) == 0) valid_bottlenecks <- sapply(bottlenecks, function(b) b$mrca_node)
  
  nodes_to_rotate <- integer(0)
  mod_tree <- tree
  
  if (strategy_choice == 1) {
    # STRATEGY 1: Single Worst Bottleneck
    nodes_to_rotate <- valid_bottlenecks[1]
    mod_tree <- ape::rotate(mod_tree, node = as.integer(nodes_to_rotate))
    
  } else if (strategy_choice == 2) {
    # STRATEGY 2: Trait-Centroid Alignment on Top 3 Bottleneck Clades
    nodes_to_rotate <- head(valid_bottlenecks, 3)
    for (nd in nodes_to_rotate) {
      mod_tree <- trait_centroid_align(mod_tree, traits, nd)
    }
    
  } else if (strategy_choice == 3) {
    # STRATEGY 3: Deep Clade Structural Reversal
    deep_nodes <- all_internal[clade_sizes[all_internal] >= 15]
    deep_nodes <- deep_nodes[!deep_nodes %in% tabu_list]
    if (length(deep_nodes) == 0) deep_nodes <- all_internal
    nodes_to_rotate <- sample(deep_nodes, min(2, length(deep_nodes)))
    for (nd in nodes_to_rotate) {
      mod_tree <- ape::rotate(mod_tree, node = as.integer(nd))
    }
    
  } else {
    # STRATEGY 4: Multi-Node Annealed Stochastic Kick (5 random rotations)
    rand_nodes <- all_internal[!all_internal %in% tabu_list]
    if (length(rand_nodes) == 0) rand_nodes <- all_internal
    k_nodes <- min(5, length(rand_nodes))
    nodes_to_rotate <- sample(rand_nodes, k_nodes)
    for (nd in nodes_to_rotate) {
      mod_tree <- ape::rotate(mod_tree, node = as.integer(nd))
    }
  }
  
  mod_tree <- sanitize_tree(mod_tree)
  return(list(tree = mod_tree, rotated_nodes = nodes_to_rotate))
}

# -----------------------------------------------------------------------------
# 5. Memory & Stagnation-Aware Strategy Optimizer
# -----------------------------------------------------------------------------
llm_strategy_optimizer <- function(tree, traits, macro_iterations = 100, tabu_size = 12) {
  cat("--- Step 1: Initial Local Hill-Climbing ---\n")
  init_run <- greedy_local_descent(tree, traits)
  current <- init_run$tree
  current_cost <- init_run$cost
  
  best <- current
  best_cost <- current_cost
  tabu_list <- integer(0)
  
  stagnation_count <- 0
  last_chosen_strat <- 0
  strat_repeat_count <- 0
  strat_counts <- c(0, 0, 0, 0)
  
  cat(sprintf("Initial Raw Cost: %.4f | Post-Greedy Descent: %.4f\n\n", 
              trait_cost(tree, traits), current_cost))
  
  system_prompt <- "You are an adaptive optimization controller for tree reordering (Optimal Leaf Ordering).
Your goal is to select strategies that escape local minima and reach global minimal trait distance (~391).

Strategies:
1 = Target Worst Trait Bottleneck (Local exploitation)
2 = Trait-Centroid Clade Alignment (Align sister clade trait medians)
3 = Deep Ancestral Clade Reversal (Flip major sub-clades >=15 tips)
4 = Multi-Node Stochastic Kick (Aggressive 5-node jump)

Analyze past stagnation and choose a strategy. Output format: 'STRATEGY: X'"
  
  for (step in 1:macro_iterations) {
    bottlenecks <- find_tree_bottlenecks(current, traits, top_k = 6)
    
    # Construct stagnation alert if stuck
    stag_warning <- ""
    if (stagnation_count > 3) {
      stag_warning <- sprintf(
        "\nWARNING: Search stuck for %d iterations! Strategy %d failed repeatedly. YOU MUST SWITCH STRATEGY!",
        stagnation_count, last_chosen_strat
      )
    }
    
    user_prompt <- sprintf(
      "Iteration: %d | Current Cost: %.4f | Best Cost: %.4f | Stagnation: %d steps\nTop Bottleneck Costs: [%s]\nStrategy Usage Counts: [S1:%d, S2:%d, S3:%d, S4:%d]%s\n\nSelect Next Strategy (1-4):",
      step, current_cost, best_cost, stagnation_count,
      paste(sapply(bottlenecks, function(b) b$dist), collapse=", "),
      strat_counts[1], strat_counts[2], strat_counts[3], strat_counts[4],
      stag_warning
    )
    
    # Query LLM
    llm_raw <- query_llm_optimizer(system_prompt, user_prompt)
    
    # Parse choice
    strat_match <- regmatches(llm_raw, regexpr("STRATEGY:\\s*([1-4])", llm_raw, ignore.case = TRUE))
    if (length(strat_match) > 0) {
      strat_choice <- as.numeric(sub("STRATEGY:\\s*", "", strat_match, ignore.case = TRUE))
    } else {
      digits <- as.numeric(regmatches(llm_raw, regexpr("[1-4]", llm_raw)))
      strat_choice <- if (length(digits) > 0) digits[1] else sample(1:4, 1)
    }
    
    # Anti-Repetition Hard Guard: Override if model gets stuck picking same stagnant strategy
    if (strat_choice == last_chosen_strat && stagnation_count > 3) {
      strat_repeat_count <- strat_repeat_count + 1
      if (strat_repeat_count >= 2) {
        other_strats <- setdiff(1:4, strat_choice)
        strat_choice <- sample(other_strats, 1)
        strat_repeat_count <- 0
      }
    } else {
      strat_repeat_count <- 0
    }
    
    last_chosen_strat <- strat_choice
    strat_counts[strat_choice] <- strat_counts[strat_choice] + 1
    
    # Execute Strategy
    strat_res <- apply_llm_strategy(current, traits, strat_choice, bottlenecks, tabu_list)
    
    # Local Descent WITH TABU LOCK (Prevents undoing the strategy move!)
    post_descent <- greedy_local_descent(strat_res$tree, traits, tabu_list = tabu_list)
    cand_cost <- post_descent$cost
    
    # Update Tabu Queue
    tabu_list <- c(strat_res$rotated_nodes, tabu_list)
    if (length(tabu_list) > tabu_size) {
      tabu_list <- tabu_list[1:tabu_size]
    }
    
    # Evaluate Move
    diff <- cand_cost - current_cost
    
    if (cand_cost < current_cost - 1e-4) {
      current <- post_descent$tree
      current_cost <- cand_cost
      
      if (current_cost < best_cost - 1e-4) {
        best <- current
        best_cost <- cand_cost
        stagnation_count <- 0
        cat(sprintf("Macro %02d | Strategy %d | NEW BEST COST: %.4f\n", step, strat_choice, best_cost))
      } else {
        cat(sprintf("Macro %02d | Strategy %d | IMPROVED (Cost: %.4f)\n", step, strat_choice, current_cost))
      }
    } else {
      stagnation_count <- stagnation_count + 1
      # Simulated Annealing Basin Jump Acceptance Probability
      sa_prob <- exp(-diff / (0.5 * (1 + 0.1 * stagnation_count)))
      
      if (runif(1) < sa_prob) {
        current <- post_descent$tree
        current_cost <- cand_cost
        cat(sprintf("Macro %02d | Strategy %d | ACCEPTED BASIN JUMP (Cost: %.4f, Δ: +%.4f)\n", 
                    step, strat_choice, cand_cost, diff))
      } else {
        cat(sprintf("Macro %02d | Strategy %d | REJECTED (Cost Δ: +%.4f)\n", step, strat_choice, diff))
      }
    }
  }
  
  return(list(tree = best, cost = best_cost))
}

# -----------------------------------------------------------------------------
# 6. Run Benchmark
# -----------------------------------------------------------------------------
res_llm <- llm_strategy_optimizer(my_tree, traits, macro_iterations = 100)

cat("\n======================================\n")
cat("Original Cost    :", trait_cost(my_tree, traits), "\n")
cat("Optimized Cost   :", res_llm$cost, "\n")
cat("======================================\n")