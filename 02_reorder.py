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

class OLOState:

    def __init__(self, leaves):

        self.leaves = leaves

        # dictionary:
        # (first leaf, last leaf) -> cost
        self.cost = {}

        # traceback:
        # (first,last) -> (orientation,left endpoints,right endpoints)
        self.trace = {}



def bar_joseph_olo(tree, traits_df):

    start_time=time.time()


    # ----------------------------
    # Leaf indexing
    # ----------------------------

    terminals = tree.get_terminals()

    names=[x.name for x in terminals]

    leaf_id={
        name:i for i,name in enumerate(names)
    }


    X=traits_df.loc[names].values


    # pairwise distances
    D=np.linalg.norm(
        X[:,None,:]-X[None,:,:],
        axis=2
    )



    # ----------------------------
    # Bottom-up DP
    # ----------------------------

    def DP(clade):


        # terminal
        if clade.is_terminal():

            idx=leaf_id[clade.name]

            state=OLOState([idx])

            state.cost[(idx,idx)] = 0

            return state



        left=DP(clade.clades[0])
        right=DP(clade.clades[1])


        state=OLOState(
            left.leaves + right.leaves
        )


        # ---------------------------------
        # left subtree before right subtree
        # ---------------------------------

        for a in left.leaves:

            for b in right.leaves:


                best=np.inf
                best_x=None
                best_y=None


                for x in left.leaves:

                    c1=left.cost.get(
                        (a,x),
                        np.inf
                    )

                    if np.isinf(c1):
                        continue


                    for y in right.leaves:

                        c2=right.cost.get(
                            (y,b),
                            np.inf
                        )


                        if np.isinf(c2):
                            continue


                        value=(
                            c1
                            +
                            D[x,y]
                            +
                            c2
                        )


                        if value < best:

                            best=value
                            best_x=x
                            best_y=y



                state.cost[(a,b)] = best

                state.trace[(a,b)] = (
                    "LR",
                    best_x,
                    best_y
                )



        # ---------------------------------
        # right subtree before left subtree
        # ---------------------------------

        for a in right.leaves:

            for b in left.leaves:


                best=np.inf
                best_x=None
                best_y=None


                for x in right.leaves:

                    c1=right.cost.get(
                        (a,x),
                        np.inf
                    )


                    if np.isinf(c1):
                        continue



                    for y in left.leaves:


                        c2=left.cost.get(
                            (y,b),
                            np.inf
                        )


                        if np.isinf(c2):
                            continue



                        value=(
                            c1
                            +
                            D[x,y]
                            +
                            c2
                        )



                        if value < best:

                            best=value
                            best_x=x
                            best_y=y



                state.cost[(a,b)] = best

                state.trace[(a,b)] = (
                    "RL",
                    best_x,
                    best_y
                )



        return state



    root_state=DP(tree.root)



    # ----------------------------
    # Choose global optimum
    # ----------------------------

    optimum=np.inf
    end_pair=None


    for pair,value in root_state.cost.items():

        if value < optimum:

            optimum=value
            end_pair=pair



    # ----------------------------
    # Traceback
    # ----------------------------

    def traceback(clade,a,b,state):


        if clade.is_terminal():

            return [
                clade.name
            ]


        left=DP_cache[id(clade.clades[0])]
        right=DP_cache[id(clade.clades[1])]


        direction,x,y = state.trace[(a,b)]


        if direction=="LR":

            left_order=traceback(
                clade.clades[0],
                a,
                x,
                left
            )

            right_order=traceback(
                clade.clades[1],
                y,
                b,
                right
            )


            return left_order+right_order


        else:

            right_order=traceback(
                clade.clades[1],
                a,
                x,
                right
            )

            left_order=traceback(
                clade.clades[0],
                y,
                b,
                left
            )


            return right_order+left_order



    # store DP states
    DP_cache={}


    def fill_cache(clade):

        state=DP(clade)

        DP_cache[id(clade)] = state

        if not clade.is_terminal():

            for child in clade.clades:
                fill_cache(child)


    fill_cache(tree.root)



    root_state=DP_cache[id(tree.root)]


    order=traceback(
        tree.root,
        end_pair[0],
        end_pair[1],
        root_state
    )


    reorder_tree_clades(
        tree,
        order
    )


    return (
        tree,
        optimum,
        time.time()-start_time
    )


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
