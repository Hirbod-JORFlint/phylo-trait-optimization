library(ape)
library(reticulate)

# Set seed for reproducibility
set.seed(42)

# -----------------------------------------------------------------------------
# 1. Setup Python & Hugging Face Model via reticulate
# -----------------------------------------------------------------------------
# Ensure python environment has transformers, torch, and accelerate
# py_install(c("transformers", "torch", "accelerate"))

transformers <- import("transformers")
torch <- import("torch")

cat("Loading Hugging Face model for tree optimization...\n")
model_id <- "Qwen/Qwen2.5-Coder-1.5B-Instruct" # Fast, lightweight instruction model

tokenizer <- transformers$AutoTokenizer$from_pretrained(model_id)
model <- transformers$AutoModelForCausalLM$from_pretrained(
  model_id,
  torch_dtype = torch$float16,
  device_map = "auto"
)

# Function to query LLM for strategic decision
query_llm_optimizer <- function(system_prompt, user_prompt) {
  messages <- list(
    list("role" = "system", "content" = system_prompt),
    list("role" = "user", "content" = user_prompt)
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
    temperature = 0.2,
    do_sample = TRUE
  )
  
  # Extract generated answer
  gen_tokens <- outputs[0, (inputs$input_ids$shape[1]:outputs$shape[2])]
  response <- tokenizer$decode(gen_tokens, skip_special_tokens = TRUE)
  return(trimws(response))
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
# 4. LLM-Guided Optimization Loop
# -----------------------------------------------------------------------------
llm_optimize_tree <- function(tree, traits, iterations = 50) {
  current <- tree
  current_cost <- trait_cost(current, traits)
  best <- current
  best_cost <- current_cost
  
  cat("--- Starting LLM-Guided Optimal Leaf Ordering ---\n")
  cat("Initial Cost:", current_cost, "\n\n")
  
  system_prompt <- "You are an expert discrete optimization agent for phylogenetic trees.
Given a list of top bottleneck nodes that cause high trait dissimilarity in leaf ordering, select the single best node ID to rotate.
Respond ONLY with a valid JSON object: {\"selected_node\": <integer_node_id>, \"reasoning\": \"<short_explanation>\"}"
  
  for (step in 1:iterations) {
    # 1. Extract structural bottlenecks
    bottlenecks <- find_tree_bottlenecks(current, traits, top_k = 5)
    
    # 2. Format context for the LLM
    candidate_nodes <- sapply(bottlenecks, function(x) x$mrca_node)
    
    prompt_text <- paste0(
      "Current Total Cost: ", round(current_cost, 2), "\n",
      "Top 5 Trait Jumps between adjacent leaves:\n",
      paste(sapply(bottlenecks, function(b) {
        paste0("- JUMP: ", b$tip1, " <-> ", b$tip2, " (Distance: ", b$dist, ") -> Controlled by MRCA Node: ", b$mrca_node)
      }), collapse = "\n"), "\n\n",
      "Which MRCA Node ID from [", paste(candidate_nodes, collapse = ", "), "] should we rotate to smooth out leaf traits?"
    )
    
    # 3. Query LLM via reticulate
    llm_raw_response <- query_llm_optimizer(system_prompt, prompt_text)
    
    # Parse node ID using regex
    selected_node <- as.numeric(regmatches(llm_raw_response, regexpr("[0-9]+", llm_raw_response)))
    
    # Fallback if parsing fails or node is invalid
    if (length(selected_node) == 0 || !(selected_node %in% candidate_nodes)) {
      selected_node <- sample(candidate_nodes, 1)
    }
    
    # 4. Perform targeted rotation
    candidate <- ape::rotate(current, node = selected_node)
    candidate_cost <- trait_cost(candidate, traits)
    
    # Accept if improved
    if (candidate_cost < current_cost) {
      current <- candidate
      current_cost <- candidate_cost
      
      if (current_cost < best_cost) {
        best <- current
        best_cost <- current_cost
      }
      cat(sprintf("Step %02d | LLM Rotated Node %d | COST IMPROVED: %.3f\n", step, selected_node, current_cost))
    } else {
      # Probabilistic acceptance (Simulated Annealing step)
      if (runif(1) < 0.15) {
        current <- candidate
        current_cost <- candidate_cost
        cat(sprintf("Step %02d | LLM Rotated Node %d | Accepted Uphill Move: %.3f\n", step, selected_node, current_cost))
      } else {
        cat(sprintf("Step %02d | LLM Rotated Node %d | Rejected Move (Cost: %.3f)\n", step, selected_node, candidate_cost))
      }
    }
  }
  
  return(list(tree = best, cost = best_cost))
}

# -----------------------------------------------------------------------------
# 5. Run & Compare
# -----------------------------------------------------------------------------
res_llm <- llm_optimize_tree(my_tree, traits, iterations = 30)

cat("\n======================================\n")
cat("Original Tree Cost :", trait_cost(my_tree, traits), "\n")
cat("LLM Optimized Cost :", res_llm$cost, "\n")
cat("======================================\n")