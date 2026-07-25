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
    max_new_tokens = as.integer(60),
    temperature    = 0.2,
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

# Fast Local Hill Climber
greedy_local_descent <- function(tree, traits, max_passes = 5) {
  current <- sanitize_tree(tree)
  current_cost <- trait_cost(current, traits)
  
  Ntip <- length(current$tip.label)
  nodes <- (Ntip + 1):(Ntip + current$Nnode)
  
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
apply_llm_strategy <- function(tree, traits, strategy_choice, bottlenecks, tabu_list) {
  Ntip <- length(tree$tip.label)
  clade_sizes <- get_clade_sizes(tree)
  all_internal <- (Ntip + 1):(Ntip + tree$Nnode)
  
  # Remove tabu nodes from candidate pools
  valid_bottlenecks <- sapply(bottlenecks, function(b) b$mrca_node)
  valid_bottlenecks <- valid_bottlenecks[!valid_bottlenecks %in% tabu_list]
  if (length(valid_bottlenecks) == 0) valid_bottlenecks <- sapply(bottlenecks, function(b) b$mrca_node)
  
  nodes_to_rotate <- integer(0)
  
  if (strategy_choice == 1) {
    # STRATEGY 1: Single Worst Bottleneck
    nodes_to_rotate <- valid_bottlenecks[1]
    
  } else if (strategy_choice == 2) {
    # STRATEGY 2: Deep Ancestral Node (Clade > 20 tips)
    deep_nodes <- all_internal[clade_sizes[all_internal] >= 20]
    deep_nodes <- deep_nodes[!deep_nodes %in% tabu_list]
    if (length(deep_nodes) == 0) deep_nodes <- all_internal
    nodes_to_rotate <- sample(deep_nodes, 1)
    
  } else if (strategy_choice == 3) {
    # STRATEGY 3: Multi-Node Macro Shakeup (Top 2-3 Bottlenecks)
    nodes_to_rotate <- head(valid_bottlenecks, 3)
    
  } else {
    # STRATEGY 4: Stochastic Exploration
    rand_nodes <- all_internal[!all_internal %in% tabu_list]
    nodes_to_rotate <- sample(rand_nodes, 1)
  }
  
  # Perform rotation
  mod_tree <- tree
  for (nd in nodes_to_rotate) {
    mod_tree <- ape::rotate(mod_tree, node = as.integer(nd))
  }
  mod_tree <- sanitize_tree(mod_tree)
  
  return(list(tree = mod_tree, rotated_nodes = nodes_to_rotate))
}

# -----------------------------------------------------------------------------
# 5. Strategy-Driven LLM Optimizer
# -----------------------------------------------------------------------------
llm_strategy_optimizer <- function(tree, traits, macro_iterations = 20, tabu_size = 8) {
  cat("--- Step 1: Initial Local Hill-Climbing ---\n")
  init_run <- greedy_local_descent(tree, traits)
  current <- init_run$tree
  current_cost <- init_run$cost
  
  best <- current
  best_cost <- current_cost
  tabu_list <- integer(0)
  
  cat(sprintf("Initial Raw Cost: %.4f | Post-Greedy Descent: %.4f\n\n", 
              trait_cost(tree, traits), current_cost))
  
  system_prompt <- "You are a meta-heuristic search controller. Choose the next strategy to jump out of local minima:
1 = Target Worst Trait Bottleneck
2 = Deep Ancestral Node Perturbation (Flip major clade)
3 = Multi-Node Macro Shakeup (Flip top 3 bottlenecks)
4 = Random Stochastic Walk

Respond ONLY with the strategy number (1, 2, 3, or 4)."
  
  for (step in 1:macro_iterations) {
    bottlenecks <- find_tree_bottlenecks(current, traits, top_k = 6)
    
    user_prompt <- sprintf(
      "Current Cost: %.4f | Stagnation Step: %d\nTop Bottleneck MRCA Nodes: [%s]\nTabu List: [%s]\n\nSelect Strategy (1-4):",
      current_cost, step, 
      paste(sapply(bottlenecks, function(b) b$mrca_node), collapse=", "),
      paste(tabu_list, collapse=", ")
    )
    
    # Query LLM
    llm_raw <- query_llm_optimizer(system_prompt, user_prompt)
    
    # Parse choice 1-4 with regex
    strat_choice <- as.numeric(regmatches(llm_raw, regexpr("[1-4]", llm_raw)))
    if (length(strat_choice) == 0) strat_choice <- sample(1:4, 1) # Fallback if parsing fails
    
    # Execute Strategy
    strat_res <- apply_llm_strategy(current, traits, strat_choice, bottlenecks, tabu_list)
    
    # Fast Local Descent on Perturbed State
    post_descent <- greedy_local_descent(strat_res$tree, traits)
    cand_cost <- post_descent$cost
    
    # Update Tabu
    tabu_list <- unique(c(strat_res$rotated_nodes, tabu_list))[1:min(tabu_size, length(tabu_list))]
    
    if (cand_cost < current_cost - 1e-4) {
      current <- post_descent$tree
      current_cost <- cand_cost
      
      if (current_cost < best_cost) {
        best <- current
        best_cost <- cand_cost
      }
      cat(sprintf("Macro %02d | Strategy %d | Rotated [%s] -> NEW BEST COST: %.4f\n", 
                  step, strat_choice, paste(strat_res$rotated_nodes, collapse=","), current_cost))
    } else {
      # Probabilistic Basin Jump Acceptance
      if (runif(1) < 0.20) {
        current <- post_descent$tree
        current_cost <- cand_cost
        cat(sprintf("Macro %02d | Strategy %d | Accepted Basin Jump (Cost: %.4f)\n", step, strat_choice, cand_cost))
      } else {
        cat(sprintf("Macro %02d | Strategy %d | Rejected Move (Cost: %.4f vs Current: %.4f)\n", step, strat_choice, cand_cost, current_cost))
      }
    }
  }
  
  return(list(tree = best, cost = best_cost))
}

# -----------------------------------------------------------------------------
# 6. Run Benchmark
# -----------------------------------------------------------------------------
res_llm <- llm_strategy_optimizer(my_tree, traits, macro_iterations = 20)

cat("\n======================================\n")
cat("Original Cost    :", trait_cost(my_tree, traits), "\n")
cat("Optimized Cost   :", res_llm$cost, "\n")
cat("======================================\n")