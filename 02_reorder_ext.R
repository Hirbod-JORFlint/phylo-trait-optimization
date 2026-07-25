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

model_id  <- "Qwen/Qwen2.5-Coder-1.5B-Instruct" # Fast, lightweight instruction model
tokenizer <- transformers$AutoTokenizer$from_pretrained(model_id)

# float16 generally requires a GPU; fall back to float32 on CPU-only machines
use_cuda <- torch$cuda$is_available()
dtype    <- if (use_cuda) torch$float16 else torch$float32

model <- transformers$AutoModelForCausalLM$from_pretrained(
  model_id,
  torch_dtype = dtype,
  device_map  = "auto"
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
# 3. Tree Trait Evaluation & Bottleneck Analysis
# -----------------------------------------------------------------------------
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

# Find adjacent tip transitions with highest continuous trait distance
find_tree_bottlenecks <- function(tree, traits, top_k = 5) {
  ord <- get_plot_order(tree)
  X <- traits[ord, , drop = FALSE]
  
  # Step-wise Euclidean distances
  dists <- sqrt(rowSums((X[-1, , drop = FALSE] - X[-nrow(X), , drop = FALSE])^2))
  
  # Order by worst distance jumps
  worst_idx <- order(dists, decreasing = TRUE)[1:top_k]
  
  mrca_matrix <- ape::mrca(tree)
  bottlenecks <- list()
  
  for (i in seq_along(worst_idx)) {
    idx <- worst_idx[i]
    tip1 <- ord[idx]
    tip2 <- ord[idx + 1]
    mrca_node <- mrca_matrix[tip1, tip2]
    
    bottlenecks[[i]] <- list(
      tip1 = tip1,
      tip2 = tip2,
      dist = round(dists[idx], 3),
      mrca_node = mrca_node
    )
  }
  return(bottlenecks)
}

# -----------------------------------------------------------------------------
# 4. Enhanced Candidate Lookahead & Feature Extraction
# -----------------------------------------------------------------------------
evaluate_candidate_moves <- function(tree, traits, candidate_nodes, current_cost, tabu_list) {
  evals <- list()
  for (i in seq_along(candidate_nodes)) {
    nd <- candidate_nodes[i]
    cand_tree <- ape::rotate(tree, node = nd)
    cand_cost <- trait_cost(cand_tree, traits)
    delta <- cand_cost - current_cost
    
    # Check clade size for context
    clade_size <- length(ape::extract.clade(tree, node = nd)$tip.label)
    is_tabu <- nd %in% tabu_list
    
    evals[[i]] <- list(
      node = nd,
      cost = cand_cost,
      delta = round(delta, 4),
      clade_size = clade_size,
      tabu = is_tabu
    )
  }
  return(evals)
}

# -----------------------------------------------------------------------------
# 5. Upgraded LLM Optimizer Loop with Tabu Search & Strategy Selection
# -----------------------------------------------------------------------------
llm_optimize_tree_v2 <- function(tree, traits, iterations = 50, tabu_size = 8) {
  current <- tree
  current_cost <- trait_cost(current, traits)
  best <- current
  best_cost <- current_cost
  
  tabu_list <- integer(0)
  stagnation_counter <- 0
  
  cat("--- Starting High-Efficiency LLM Tabu Optimizer ---\n")
  cat("Initial Cost:", round(current_cost, 4), "\n\n")
  
  system_prompt <- "You are an expert combinatorial optimizer for tree leaf ordering.
Your task is to analyze candidate node rotations, their cost deltas (negative = improvement), clade sizes, and tabu status.
Select the best node to rotate, prioritizing negative deltas and non-tabu nodes.
If all deltas are positive, pick a non-tabu ancestral node with a large clade size to escape the local minimum.
Respond ONLY with a JSON object: {\"selected_node\": <integer_node_id>, \"reasoning\": \"<short_reason>\"}"
  
  for (step in 1:iterations) {
    # 1. Extract bottleneck nodes (MRCA of worst trait jumps)
    bottlenecks <- find_tree_bottlenecks(current, traits, top_k = 6)
    candidate_nodes <- unique(sapply(bottlenecks, function(x) x$mrca_node))
    
    # 2. Compute exact lookahead delta cost for each candidate
    candidate_evals <- evaluate_candidate_moves(current, traits, candidate_nodes, current_cost, tabu_list)
    
    # 3. Format candidate options for LLM prompt
    candidates_formatted <- paste(sapply(candidate_evals, function(e) {
      status <- if (e$tabu) " [TABU - RECENTLY VISITED]" else ""
      sprintf("- Node %d: Delta = %.4f | Clade Size = %d tips%s", 
              e$node, e$delta, e$clade_size, status)
    }), collapse = "\n")
    
    user_prompt <- sprintf(
      "Current Cost: %.4f | Stagnation Level: %d\nRecently Tried Nodes (Tabu): [%s]\nCandidate Moves:\n%s\n\nWhich Node ID should be rotated?",
      current_cost, stagnation_counter, paste(tabu_list, collapse = ", "), candidates_formatted
    )
    
    # 4. Query LLM Policy
    llm_raw_response <- query_llm_optimizer(system_prompt, user_prompt)
    
    # Parse selected node ID
    selected_node <- as.numeric(regmatches(llm_raw_response, regexpr("[0-9]+", llm_raw_response)))
    
    # Fallback policy: Pick best non-tabu candidate if parsing fails
    valid_candidates <- candidate_nodes[!candidate_nodes %in% tabu_list]
    if (length(valid_candidates) == 0) valid_candidates <- candidate_nodes # Reset tabu if all blocked
    
    if (length(selected_node) == 0 || !(selected_node %in% candidate_nodes)) {
      # Fallback to greedy non-tabu node
      deltas <- sapply(candidate_evals, function(e) ifelse(e$tabu, Inf, e$delta))
      selected_node <- candidate_evals[[which.min(deltas)]]$node
    }
    
    # 5. Apply selected rotation
    candidate <- ape::rotate(current, node = selected_node)
    candidate_cost <- trait_cost(candidate, traits)
    delta <- candidate_cost - current_cost
    
    # 6. Update Search State & Tabu Memory
    tabu_list <- c(selected_node, tabu_list)
    if (length(tabu_list) > tabu_size) tabu_list <- tabu_list[1:tabu_size]
    
    if (candidate_cost < current_cost) {
      current <- candidate
      current_cost <- candidate_cost
      stagnation_counter <- 0
      
      if (current_cost < best_cost) {
        best <- current
        best_cost <- current_cost
      }
      cat(sprintf("Step %02d | Rotated Node %d | COST IMPROVED: %.4f (Delta: %.4f)\n", 
                  step, selected_node, current_cost, delta))
    } else {
      stagnation_counter <- stagnation_counter + 1
      
      # Simulated annealing acceptance probability for uphill moves
      temp <- max(1.0, 10.0 * (0.95^step))
      accept_prob <- exp(-delta / temp)
      
      if (stagnation_counter < 5 && runif(1) < accept_prob) {
        current <- candidate
        current_cost <- candidate_cost
        cat(sprintf("Step %02d | Rotated Node %d | Accepted Uphill Move: %.4f (Delta: +%.4f)\n", 
                    step, selected_node, current_cost, delta))
      } else {
        cat(sprintf("Step %02d | Rotated Node %d | Rejected (Delta: +%.4f)\n", 
                    step, selected_node, delta))
      }
    }
  }
  
  return(list(tree = best, cost = best_cost))
}

# -----------------------------------------------------------------------------
# 3. Run Optimization Benchmark
# -----------------------------------------------------------------------------
res_llm <- llm_optimize_tree_v2(my_tree, traits, iterations = 50, tabu_size = 8)

cat("\n======================================\n")
cat("Original Tree Cost :", trait_cost(my_tree, traits), "\n")
cat("LLM Tabu Opt Cost  :", res_llm$cost, "\n")
cat("======================================\n")