import random
import math
import numpy as np
import pandas as pd
import copy
import time
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from scipy.cluster.hierarchy import leaves_list, optimal_leaf_ordering
from scipy.spatial.distance import pdist, squareform

# Set seeds
np.random.seed(42)
random.seed(42)

n_tips = 500
tip_labels = [f"Tip{i}" for i in range(1, n_tips + 1)]

# Helper function to generate a random binary tree topology
def make_random_clade(labels):
    if len(labels) == 1:
        return Clade(name=labels[0], branch_length=random.uniform(0.1, 1.0))

    # Split labels randomly into two groups
    split = random.randint(1, len(labels) - 1)
    left_labels = labels[:split]
    right_labels = labels[split:]

    left_child = make_random_clade(left_labels)
    right_child = make_random_clade(right_labels)

    return Clade(
        clades=[left_child, right_child], branch_length=random.uniform(0.1, 1.0)
    )


# Build Biopython Tree
root_clade = make_random_clade(tip_labels)
my_tree = Tree(root=root_clade, rooted=True)

# Generate Trait Matrix
trait_labels = [f"Trait{i}" for i in range(1, 11)]

traits = pd.DataFrame(
    np.random.uniform(0.0, 5.0, size=(n_tips, 10)),
    index=tip_labels,
    columns=trait_labels,
)

print(traits.head())

# Plot with Matplotlib
fig, ax = plt.subplots(figsize=(10, 12))

Phylo.draw(
    my_tree,
    axes=ax,
    do_show=False,  # Prevents showing the plot immediately so you can customize
    label_func=lambda clade: "",  # Hides individual leaf labels for large trees
)

ax.set_title("500-Tip Random Phylogeny", fontsize=14)
ax.set_xlabel("Branch Length")
ax.set_ylabel("")  # Clean up Y axis labels

# Only display label if it ends with "0" (e.g., Tip10, Tip20)
#label_func = lambda c: (
#    c.name if c.is_terminal() and c.name.endswith("0") else ""
#)

# Replace plt.show() with:
fig.savefig("phylogeny_500tips.png", dpi=300, bbox_inches="tight")
# Or vector graphics:
fig.savefig("phylogeny_500tips.svg", format="svg")

# Color the first main branch red
my_tree.root.clades[0].color = "red"

#plt.tight_layout()
#plt.show()


# 2. Set initial vertical view limits (e.g., show tips 0 to 30)
view_window = 30
ax.set_ylim(0, view_window)


# 3. Add scroll event handler
def on_scroll(event):
    if event.inaxes != ax:
        return

    cur_ymin, cur_ymax = ax.get_ylim()
    step = (cur_ymax - cur_ymin) * 0.2  # Scroll distance per step

    if event.button == "up":  # Scroll up
        new_ymin = min(n_tips - view_window, cur_ymin + step)
        new_ymax = min(n_tips, cur_ymax + step)
    elif event.button == "down":  # Scroll down
        new_ymin = max(0, cur_ymin - step)
        new_ymax = max(view_window, cur_ymax - step)
    else:
        return

    ax.set_ylim(new_ymin, new_ymax)
    fig.canvas.draw_idle()


# Connect mouse wheel event to plot
fig.canvas.mpl_connect("scroll_event", on_scroll)

ax.set_title(
    "Interactive Scrollable Tree (Use Mouse Wheel to Scroll Up/Down)",
    fontsize=12,
)
plt.tight_layout()
#plt.show()

# ---------------------------------------------------------
# 1. Shared Helper Functions
# ---------------------------------------------------------


def get_plot_order(tree):
  """Extracts tip labels in their exact visual/traversal order (left-to-right)."""
  return [leaf.name for leaf in tree.get_terminals()]


def trait_cost(tree, traits_df):
  """Calculates total Euclidean distance between adjacent leaves in tree layout."""
  ord_labels = get_plot_order(tree)
  X = traits_df.loc[ord_labels].values
  diffs = np.diff(X, axis=0)
  return float(np.sum(np.linalg.norm(diffs, axis=1)))


def tree_to_linkage(tree, tip_labels):
  """Converts a Bio.Phylo Tree into a SciPy linkage matrix Z (N-1, 4)."""
  label_to_id = {name: i for i, name in enumerate(tip_labels)}
  Z_rows = []
  next_id = len(tip_labels)

  def _traverse(clade):
    nonlocal next_id
    if clade.is_terminal():
      return label_to_id[clade.name], 1, 0.0

    children_info = [_traverse(c) for c in clade.clades]

    curr = children_info
    while len(curr) > 1:
      c1_id, c1_count, c1_h = curr.pop(0)
      c2_id, c2_count, c2_h = curr.pop(0)

      new_count = c1_count + c2_count
      new_h = max(c1_h, c2_h) + 1.0

      Z_rows.append([c1_id, c2_id, new_h, new_count])
      curr.insert(0, (next_id, new_count, new_h))
      next_id += 1

    return curr[0]

  _traverse(tree.root)
  return np.array(Z_rows, dtype=float)


def reorder_tree_clades(tree, optimal_tip_order):
  """Recursively rotates internal clades in-place to match optimal linear leaf ordering."""
  pos_map = {name: i for i, name in enumerate(optimal_tip_order)}

  def _sort_clade(clade):
    if clade.is_terminal():
      return
    for child in clade.clades:
      _sort_clade(child)
    clade.clades.sort(
        key=lambda child: min(
            pos_map[term.name] for term in child.get_terminals()
        )
    )

  _sort_clade(tree.root)
  
# ---------------------------------------------------------
# 3. Solve Optimal Leaf Ordering
# ---------------------------------------------------------

def olo_tree(tree, traits_df):
  """Solves Optimal Leaf Ordering (OLO) globally using SciPy's Dynamic Programming.

  Returns: (optimized_tree, final_cost, elapsed_time_sec)
  """
  start_time = time.time()

  # Extract current tip order
  tip_labels = [t.name for t in tree.get_terminals()]

  # 1. Convert Bio.Phylo tree to SciPy linkage matrix
  Z = tree_to_linkage(tree, tip_labels)

  # 2. Extract traits matrix corresponding to initial tip order
  X = traits_df.loc[tip_labels].values

  # 3. Solve OLO via Dynamic Programming (Fast Exact O(N^3) solver)
  Z_opt = optimal_leaf_ordering(Z, X, metric="euclidean")

  # 4. Extract optimal leaf index sequence and reorder tree clades
  optimal_indices = leaves_list(Z_opt)
  optimal_tip_order = [tip_labels[i] for i in optimal_indices]
  reorder_tree_clades(tree, optimal_tip_order)

  # 5. Evaluate final cost and execution time
  final_cost = trait_cost(tree, traits_df)
  elapsed_time = time.time() - start_time

  return tree, final_cost, elapsed_time


## simulated annealing in python
def anneal_tree(
    tree,
    traits_df,
    max_iter=200000,
    temperature=10.0,
    cooling=0.99995,  # Fixed: changed 0.099995 to 0.99995
    patience=10000,
):
  internal_nodes = [
      clade for clade in tree.find_clades() if not clade.is_terminal()
  ]

  current_cost = trait_cost(tree, traits_df)
  best_cost = current_cost
  no_change = 0

  print(f"Starting cost: {current_cost:.4f}")

  for i in range(1, max_iter + 1):
    target_node = random.choice(internal_nodes)
    target_node.clades.reverse()

    candidate_cost = trait_cost(tree, traits_df)
    delta = candidate_cost - current_cost

    # Safe acceptance check: prevent division by zero or float underflow
    if delta < 0:
      accept = True
    elif temperature < 1e-12:
      accept = False
    else:
      try:
        accept = random.random() < math.exp(-delta / temperature)
      except (OverflowError, ZeroDivisionError):
        accept = False

    # Apply move decision
    if accept:
      current_cost = candidate_cost
      if current_cost < best_cost:
        best_cost = current_cost
        no_change = 0
      else:
        no_change += 1
    else:
      target_node.clades.reverse()  # Undo rotation
      no_change += 1

    temperature *= cooling

    if no_change > patience:
      print(f"Early stopping at iteration {i} (patience reached).")
      return tree, best_cost, i

  return tree, best_cost, max_iter

# ---------------------------------------------------------
# Comparative Execution
# ---------------------------------------------------------

# 1. Standardize trait matrix (Recommended for multi-dimensional data)
traits_scaled = (traits - traits.mean()) / traits.std()

# Prepare deep copies of the tree so both start from the exact same layout
tree_for_sa = copy.deepcopy(my_tree)
tree_for_olo = copy.deepcopy(my_tree)

# Compute starting cost
init_cost = trait_cost(my_tree, traits_scaled)

print("=" * 55)
print(f"INITIAL UNOPTIMIZED COST: {init_cost:.4f}")
print("=" * 55)

# ---------------------------------------------------------
# Method 1: Dynamic Programming OLO
# ---------------------------------------------------------
olo_result_tree, olo_cost, olo_time = olo_tree(tree_for_olo, traits_scaled)

print("\n--- [Method 1: Dynamic Programming (OLO)] ---")
print(f"Final Cost:     {olo_cost:.4f}")
print(f"Cost Reduction: {(1 - olo_cost / init_cost) * 100:.2f}%")
print(f"Execution Time: {olo_time:.3f} seconds")

# ---------------------------------------------------------
# Method 2: Simulated Annealing
# ---------------------------------------------------------
print("\n--- [Method 2: Simulated Annealing] ---")
sa_start = time.time()
sa_result_tree, sa_cost, sa_iters = anneal_tree(
    tree_for_sa,
    traits_scaled,
    max_iter=200000,
    temperature=10.0,
    cooling=0.099995,
    patience=10000,
)
sa_time = time.time() - sa_start

print(f"Final Cost:     {sa_cost:.4f}")
print(f"Cost Reduction: {(1 - sa_cost / init_cost) * 100:.2f}%")
print(f"Execution Time: {sa_time:.3f} seconds ({sa_iters} iterations)")

# ---------------------------------------------------------
# Summary Comparison
# ---------------------------------------------------------
gap = ((sa_cost - olo_cost) / olo_cost) * 100
print("\n" + "=" * 55)
print("BENCHMARK SUMMARY:")
print(f"• Exact Global Optimum (OLO Cost): {olo_cost:.4f}")
print(f"• Heuristic Estimate  (SA Cost):  {sa_cost:.4f}")
if gap > 0:
  print(f"• Simulated Annealing was {gap:.2f}% worse than exact OLO.")
else:
  print("• Simulated Annealing matched the exact global OLO solution.")
print(f"• Speedup: OLO was {sa_time / olo_time:.1f}x faster than SA.")
print("=" * 55)
