# install.packages("R.utils")
library(ape)
library(phytools)
library(RMaCzek)
library(seriation)


set.seed(42)

my_tree <- rtree(
  100,
  tip.label=paste0("Tip",1:100)
)

traits <- matrix(
  runif(100*6,0,5),
  nrow=100,
  dimnames=list(my_tree$tip.label,
                paste0("Trait",1:6))
)



# Cost function


# get_plot_order <- function(tree){
# 
#     oldpar <- par(no.readonly = TRUE)
# 
#     on.exit(par(oldpar))
# 
#     plot(
#         tree,
#         plot=FALSE
#     )
# 
#     pp <- get(
#         "last_plot.phylo",
#         envir=.PlotPhyloEnv
#     )
# 
#     n <- length(tree$tip.label)
# 
#     tree$tip.label[
#         order(pp$yy[1:n])
#     ]
# }
# 
# trait_cost <- function(tree, traits){
# 
#     ord <- get_plot_order(tree)
# 
#     X <- traits[ord,,drop=FALSE]
# 
#     sum(
#         sqrt(
#             rowSums(
#                 (X[-1,]-X[-nrow(X),])^2
#             )
#         )
#     )
# }

# Fast tree tip order (no plotting)


get_plot_order <- function(tree){
  
  n <- length(tree$tip.label)
  
  # children of each internal node
  children <- split(
    tree$edge[,2],
    tree$edge[,1]
  )
  
  root <- setdiff(
    tree$edge[,1],
    tree$edge[,2]
  )[1]
  
  
  # recursive depth-first traversal
  traverse <- function(node){
    
    # tip
    if(node <= n){
      return(node)
    }
    
    unlist(
      lapply(
        children[[as.character(node)]],
        traverse
      )
    )
  }
  
  
  idx <- traverse(root)
  
  tree$tip.label[idx]
}


trait_cost <- function(tree, traits){
  
  ord <- get_plot_order(tree)
  
  X <- traits[ord,,drop=FALSE]
  
  sum(
    sqrt(
      rowSums(
        (X[-1,]-X[-nrow(X),])^2
      )
    )
  )
  
}


# Random rotation of one internal node


random_rotate <- function(tree){
  
  Ntip <- length(tree$tip.label)
  
  internal_nodes <- (Ntip+1):(Ntip+tree$Nnode)
  
  node <- sample(internal_nodes,1)
  
  tree2 <- tree
  
  children <- which(tree2$edge[,1]==node)
  
  if(length(children)==2){
    
    tree2$edge[children,2] <-
      rev(tree2$edge[children,2])
    
  }
  
  tree2
  
}




# Simulated annealing optimizer


anneal_tree_order <- function(
    tree,
    traits,
    max_iter=100000,
    temperature=1,
    cooling=0.99995,
    patience=5000
){
  
  current_tree <- tree
  
  current_cost <- trait_cost(
    current_tree,
    traits
  )
  
  cat(
    "Starting cost:",
    current_cost,
    "\n"
  )
  
  best_tree <- current_tree
  
  best_cost <- current_cost
  
  
  no_improve <- 0
  
  
  for(iter in 1:max_iter){
    
    
    candidate <- random_rotate(
      current_tree
    )
    
    
    candidate_cost <- trait_cost(
      candidate,
      traits
    )
    
    
    delta <- candidate_cost-current_cost
    
    
    
    # Metropolis acceptance rule
    
    
    accept <- FALSE
    
    if(delta < 0){
      
      accept <- TRUE
      
    } else {
      
      probability <- exp(
        -delta/temperature
      )
      
      if(runif(1)<probability){
        accept <- TRUE
      }
      
    }
    
    
    if(accept){
      
      current_tree <- candidate
      
      current_cost <- candidate_cost
      
    }
    
    
    
    # Save best
    
    
    if(current_cost < best_cost){
      
      best_tree <- current_tree
      
      best_cost <- current_cost
      
      no_improve <- 0
      
    } else {
      
      no_improve <- no_improve + 1
      
    }
    
    
    
    # Cooling
    
    
    temperature <- temperature * cooling
    
    
    
    
    # Early stopping
    
    
    if(no_improve >= patience){
      
      message(
        "Stopped early at iteration ",
        iter
      )
      
      break
      
    }
    
    
    if(temperature < 1e-8){
      
      message(
        "Temperature exhausted"
      )
      
      break
      
    }
    
    
    if(iter %% 1000 == 0){
      
      message(
        "Iteration ",
        iter,
        " Best cost = ",
        best_cost
      )
      
    }
    
  }
  
  
  list(
    tree=best_tree,
    cost=best_cost,
    iterations=iter
  )
  
}




# Run optimization


result <- anneal_tree_order(
  my_tree,
  traits,
  max_iter=200000,
  temperature=10,
  cooling=0.99997,
  patience=10000
)



# Results


optimized_tree <- result$tree

init_cost <- trait_cost(my_tree,traits)
cat(
  "Initial cost:",
  init_cost,
  "\n"
)


cat(
  "Final cost:",
  result$cost,
  "\n"
)

cat(
  "Iterations:",
  result$iterations,
  "\n"
)



# Plot example


par(mfrow=c(1,2))

plot(my_tree,
     show.tip.label=FALSE,
     main="Original")

plot(optimized_tree,
     show.tip.label=FALSE,
     main="Optimized")


# use RMaCzek now

ord <- get_plot_order(optimized_tree)

cz <- czek_matrix(
  traits,
  order = NULL  # tell RMaCzek not to re-seriate; I already picked the order from anealing
)

par(mfrow = c(1, 3))
plot(my_tree, show.tip.label = FALSE, main = "Original")
plot(optimized_tree, show.tip.label = FALSE, main = "Optimized")
plot(cz, main = "Czekanowski Diagram")

# helper function
trait_cost_from_order <- function(ord, traits){
  
  X <- traits[ord, , drop = FALSE]
  
  sum(
    sqrt(
      rowSums(
        (X[-1, ] - X[-nrow(X), ])^2
      )
    )
  )
}

cz_auto <- czek_matrix(traits, order = 'OLO')

# check what czek_matrix actually returns to get the order
str(cz_auto)

# I assume the second code is what needed to calculate the cost

ord_auto <- rownames(cz_auto)
ord_auto <- names(attr(cz_auto, "order"))


head(ord_auto)

cost_auto <- trait_cost_from_order(ord_auto, traits)
cost_Mine <- result$cost

cat("RMaCzek seriation cost:", cost_auto, "\n")
cat("My annealing cost:   ", cost_Mine, "\n")

library(R.utils)

# 1. Gather all distance-based seriation methods
methods_to_test <- seriation::list_seriation_methods("dist")

# Remove Enumerate (completely impractical for large datasets)
methods_to_test <- setdiff(methods_to_test, "Enumerate")

# Add RMaCzek's native genetic algorithm
methods_to_test <- c(methods_to_test, "ga")

# 2. Initialize results table
benchmark_results <- data.frame(
  Method = character(),
  Cost = numeric(),
  Status = character(),
  stringsAsFactors = FALSE
)

# 3. Benchmark each method with a 30-second timeout
for (m in methods_to_test) {
  
  message("Evaluating method: ", m, "...")
  
  res <- tryCatch({
    
    withTimeout({
      
      # Run seriation
      cz_res <- czek_matrix(
        traits,
        order = m,
        scale_data = FALSE
      )
      
      # Extract ordering
      ord_idx <- attr(cz_res, "order")
      
      # Compute cost
      current_cost <- trait_cost_from_order(ord_idx, traits)
      
      data.frame(
        Method = m,
        Cost = current_cost,
        Status = "Success",
        stringsAsFactors = FALSE
      )
      
    }, timeout = 30, onTimeout = "error")
    
  }, TimeoutException = function(e) {
    
    data.frame(
      Method = m,
      Cost = NA,
      Status = "Timeout (>30s)",
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    
    data.frame(
      Method = m,
      Cost = NA,
      Status = paste("Error:", conditionMessage(e)),
      stringsAsFactors = FALSE
    )
    
  })
  
  benchmark_results <- rbind(benchmark_results, res)
}

# 4. Append simulated annealing result
annealing_res <- data.frame(
  Method = "My_Simulated_Annealing",
  Cost = result$cost,
  Status = "Success",
  stringsAsFactors = FALSE
)

benchmark_results <- rbind(benchmark_results, annealing_res)

# 5. Sort leaderboard
benchmark_results <- benchmark_results[
  order(benchmark_results$Cost, na.last = TRUE),
]

print(benchmark_results, row.names = FALSE)