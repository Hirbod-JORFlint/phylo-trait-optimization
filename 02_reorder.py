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

## Bar_joseph dynamic programming
def bar_joseph_olo(tree, traits_df):

    start = time.time()

    leaves = [x.name for x in tree.get_terminals()]
    leaf_id = {x:i for i,x in enumerate(leaves)}

    X = traits_df.loc[leaves].values

    D = np.linalg.norm(
        X[:,None,:] - X[None,:,:],
        axis=2
    )


    class NodeDP:
        def __init__(self):
            self.leaves = []
            self.cost = {}
            self.choice = {}


    def solve(clade):

        node = NodeDP()


        # --------------------------
        # Terminal node
        # --------------------------
        if clade.is_terminal():

            i = leaf_id[clade.name]

            node.leaves = [i]
            node.cost[(i,i)] = 0

            return node



        A = solve(clade.clades[0])
        B = solve(clade.clades[1])

        node.leaves = A.leaves + B.leaves



        # --------------------------
        # Orientation function
        # --------------------------
        def merge(left,right,orientation):

            cost = {}
            choice = {}


            for i in left.leaves:
                for j in right.leaves:

                    best = np.inf
                    bx = None
                    by = None


                    for x in left.leaves:

                        if (i,x) not in left.cost:
                            continue

                        for y in right.leaves:

                            if (y,j) not in right.cost:
                                continue

                            value = (
                                left.cost[(i,x)]
                                +
                                D[x,y]
                                +
                                right.cost[(y,j)]
                            )


                            if value < best:
                                best=value
                                bx=x
                                by=y


                    cost[(i,j)] = best
                    choice[(i,j)] = (
                        orientation,
                        bx,
                        by
                    )

            return cost, choice



        # left -> right
        LR_cost, LR_choice = merge(
            A,B,"LR"
        )


        # right -> left
        RL_cost, RL_choice = merge(
            B,A,"RL"
        )


        # --------------------------
        # Select better orientation
        # --------------------------

        for i,j in LR_cost:

            if LR_cost[(i,j)] <= RL_cost.get((j,i),np.inf):

                node.cost[(i,j)] = LR_cost[(i,j)]

                node.choice[(i,j)] = LR_choice[(i,j)]

            else:

                node.cost[(j,i)] = RL_cost[(j,i)]

                node.choice[(j,i)] = RL_choice[(j,i)]


        return node



    root_dp = solve(tree.root)


    # global optimum

    best = np.inf
    best_pair=None

    for pair,value in root_dp.cost.items():

        if value < best:
            best=value
            best_pair=pair



    # --------------------------
    # Recover ordering
    # --------------------------

    def recover(clade, start_leaf, end_leaf):

        if clade.is_terminal():
            return [clade.name]


        A,B = clade.clades


        # Try both possible child orientations

        order1 = recover(
            A,
            start_leaf,
            end_leaf
        ) + recover(
            B,
            start_leaf,
            end_leaf
        )


        order2 = recover(
            B,
            start_leaf,
            end_leaf
        ) + recover(
            A,
            start_leaf,
            end_leaf
        )


        c1 = trait_cost_from_order(
            order1,
            traits_df
        )

        c2 = trait_cost_from_order(
            order2,
            traits_df
        )


        if c1 <= c2:
            clade.clades=[A,B]
            return order1

        else:
            clade.clades=[B,A]
            return order2



    def trait_cost_from_order(order,traits):

        Z = traits.loc[order].values

        return np.sum(
            np.linalg.norm(
                np.diff(Z,axis=0),
                axis=1
            )
        )


    order = recover(
        tree.root,
        best_pair[0],
        best_pair[1]
    )


    reorder_tree_clades(
        tree,
        order
    )


    elapsed=time.time()-start


    return tree,best,elapsed


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
# Method 1: Bar-Joseph Dynamic Programming OLO
# ---------------------------------------------------------

olo_result_tree, olo_cost, olo_time = bar_joseph_olo(
    tree_for_olo,
    traits_scaled
)

print("\n--- [Method 1: Bar-Joseph Dynamic Programming OLO] ---")
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
    cooling=0.99995,
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
print(f"• Speedup: OLO was {sa_time / olo_time:.1f}x faster than SA.")
print("=" * 55)
