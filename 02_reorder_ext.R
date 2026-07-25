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
# 3. Helper Functions & Tree Sanitization
# -----------------------------------------------------------------------------

# Prevents ape C-level matrix type coercion errors during repeated rotations
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

# Pre-calculate clade sizes without calling heavy ape::extract.clade
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

# -----------------------------------------------------------------------------
# 4. Fast Native R Greedy Hill-Climber (Milliseconds)
# -----------------------------------------------------------------------------
greedy_local_descent <- function(tree, traits, max_passes = 3) {
  current <- sanitize_tree(tree)
  current_cost <- trait_cost(current, traits)
  
  Ntip <- length(current$tip.label)
  nodes <- (Ntip + 1):(Ntip + current$Nnode)
  
  improved <- TRUE
  pass <- 0
  
  while (improved && pass < max_passes) {
    improved <- FALSE
    pass <- pass + 1
    
    # Evaluate all internal nodes sequentially
    for (nd in nodes) {
      cand <- ape::rotate(current, node = as.integer(nd))
      cand_cost <- trait_cost(cand, traits)
      
      if (cand_cost < current_cost - 1e-5) {
        current <- cand
        current_cost <- cand_cost
        improved <- TRUE
      }
    }
    current <- sanitize_tree(current)
  }
  
  return(list(tree = current, cost = current_cost))
}

# -----------------------------------------------------------------------------
# 5. Find Trait Bottlenecks for LLM Context
# -----------------------------------------------------------------------------
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
      tip1 = tip1,
      tip2 = tip2,
      dist = round(dists[idx], 3),
      mrca_node = mrca_node,
      clade_size = clade_sizes[mrca_node]
    )
  }
  return(bottlenecks)
}

# -----------------------------------------------------------------------------
#  6. LLM Meta-Optimizer Loop (Local Descent + Strategic Basin Hopping)
# -----------------------------------------------------------------------------
llm_basin_hopping_optimizer <- function(tree, traits, macro_iterations = 15, tabu_size = 5) {
  # Step 0: Immediate initial local descent
  cat("--- Running Initial Local Hill-Climbing Descent ---\n")
  init_descent <- greedy_local_descent(tree, traits)
  current <- init_descent$tree
  current_cost <- init_descent$cost
  
  best <- current
  best_cost <- current_cost
  
  tabu_list <- integer(0)
  
  cat(sprintf("Initial Raw Cost: %.4f | Post-Local Descent Cost: %.4f\n\n", 
              trait_cost(tree, traits), current_cost))
  
  system_prompt <- "You are an expert meta-heuristic optimization agent for phylogenetic optimal leaf ordering.
The search is currently trapped in a local minimum.
Your job is to select TWO distinct MRCA node IDs from the bottleneck list to rotate SIMULTANEOUSLY as a macro-perturbation step to jump out of the local minimum.
Prefer nodes with large clade sizes or non-tabu nodes.
Respond ONLY with JSON: {\"primary_node\": <int>, \"secondary_node\": <int>, \"reasoning\": \"<short>\"}"
  
  for (step in 1:macro_iterations) {
    # 1. Identify bottlenecks
    bottlenecks <- find_tree_bottlenecks(current, traits, top_k = 6)
    cand_nodes <- unique(sapply(bottlenecks, function(x) x$mrca_node))
    
    # 2. Format context for LLM
    bottleneck_text <- paste(sapply(bottlenecks, function(b) {
      status <- if (b$mrca_node %in% tabu_list) " [TABU]" else ""
      sprintf("- MRCA Node %d: Trait Jump = %.3f | Clade Size = %d tips%s", 
              b$mrca_node, b$dist, b$clade_size, status)
    }), collapse = "\n")
    
    user_prompt <- sprintf(
      "Current Local Min Cost: %.4f\nRecently Perturbed Nodes (Tabu): [%s]\n\nBottleneck Candidates:\n%s\n\nSelect TWO node IDs to rotate together to jump out of this basin.",
      current_cost, paste(tabu_list, collapse = ", "), bottleneck_text
    )
    
    # 3. Query LLM Policy
    llm_raw_response <- query_llm_optimizer(system_prompt, user_prompt)
    
    # Extract numbers
    parsed_nodes <- as.integer(regmatches(llm_raw_response, gregexpr("[0-9]+", llm_raw_response))[[1]])
    parsed_nodes <- parsed_nodes[parsed_nodes %in% cand_nodes]
    
    # Fallback if parsing fails
    if (length(parsed_nodes) < 2) {
      valid_nodes <- cand_nodes[!cand_nodes %in% tabu_list]
      if (length(valid_nodes) < 2) valid_nodes <- cand_nodes
      parsed_nodes <- sample(valid_nodes, min(2, length(valid_nodes)))
    }
    
    n1 <- parsed_nodes[1]
    n2 <- parsed_nodes[2]
    
    # 4. Apply LLM Macro-Perturbation
    perturbed_tree <- current
    perturbed_tree <- ape::rotate(perturbed_tree, node = n1)
    if (!is.na(n2) && n2 != n1) {
      perturbed_tree <- ape::rotate(perturbed_tree, node = n2)
    }
    perturbed_tree <- sanitize_tree(perturbed_tree)
    
    # 5. Run Fast Local Hill-Climbing on the Perturbed Basin
    post_descent <- greedy_local_descent(perturbed_tree, traits)
    candidate_tree <- post_descent$tree
    candidate_cost <- post_descent$cost
    
    # 6. Update Tabu & Accept/Reject Logic
    tabu_list <- unique(c(n1, n2, tabu_list))[1:min(tabu_size, length(tabu_list))]
    
    if (candidate_cost < current_cost - 1e-4) {
      current <- candidate_tree
      current_cost <- candidate_cost
      
      if (current_cost < best_cost) {
        best <- current
        best_cost <- current_cost
      }
      cat(sprintf("Macro Step %02d | LLM Perturbed Nodes [%d, %d] -> NEW BEST COST: %.4f\n", 
                  step, n1, n2, current_cost))
    } else {
      # Probabilistic acceptance of basin jump
      delta <- candidate_cost - current_cost
      if (runif(1) < exp(-delta / 5.0)) {
        current <- candidate_tree
        current_cost <- candidate_cost
        cat(sprintf("Macro Step %02d | LLM Perturbed Nodes [%d, %d] -> Accepted Basin Jump (Cost: %.4f)\n", 
                    step, n1, n2, current_cost))
      } else {
        cat(sprintf("Macro Step %02d | LLM Perturbed Nodes [%d, %d] -> Rejected (Cost: %.4f vs Current: %.4f)\n", 
                    step, n1, n2, candidate_cost, current_cost))
      }
    }
  }
  
  return(list(tree = best, cost = best_cost))
}

# -----------------------------------------------------------------------------
# 7. Run Execution Benchmark
# -----------------------------------------------------------------------------
res_llm <- llm_basin_hopping_optimizer(my_tree, traits, macro_iterations = 15, tabu_size = 6)

cat("\n======================================\n")
cat("Original Tree Cost   :", trait_cost(my_tree, traits), "\n")
cat("LLM Hybrid Opt Cost  :", res_llm$cost, "\n")
cat("======================================\n")