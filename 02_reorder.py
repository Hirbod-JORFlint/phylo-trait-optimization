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
    labels = list(labels)

    if len(labels) == 1:
        return Clade(
            name=labels[0],
            branch_length=random.uniform(0.1, 1.0)
        )

    random.shuffle(labels)

    split = random.randint(1, len(labels)-1)

    left = make_random_clade(labels[:split])
    right = make_random_clade(labels[split:])

    return Clade(
        clades=[left, right],
        branch_length=random.uniform(0.1, 1.0)
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

## Bar-Joseph dynamic programming (exact OLO)
#
# M_v(u, w) = minimum adjacent-sum cost of ordering every leaf in subtree v
#             so that leaf u is at one end and leaf w is at the other.
#
# Leaf (base case):  M_v(u, u) = 0.
#
# Internal node v with children L, R  (u ∈ L, w ∈ R):
#   M_v(u, w) = min_{b ∈ L, c ∈ R} [ M_L(u, b) + D(b, c) + M_R(c, w) ]
#
# R-before-L layouts are recovered by reversal symmetry:
#   the reversed ordering has the same adjacent-sum cost because D is symmetric,
#   so M_v(w, u) = M_v(u, w).  We only store the L-before-R rectangular table.
#
# Complexity: O(n³) total across the tree via factored minimisation:
#   Pass 1 –  g(b, w) = min_{c ∈ R} [ D(b, c) + M_R(c, w) ]   O(|L|·|R|²)
#   Pass 2 –  M_v(u, w) = min_{b ∈ L} [ M_L(u, b) + g(b, w) ]  O(|L|²·|R|)


class OLOState:
    """DP state for one subtree.

    Attributes
    ----------
    leaves : list[int]
        Global leaf indices in this subtree, in an arbitrary canonical order.
    left_leaves / right_leaves : list[int] | None
        Leaf indices of the left / right child (None for a tip).
    cost : dict[(int,int) -> float]
        M_v(u, w) for u in left_leaves, w in right_leaves (or u == w for a tip).
    trace : dict[(int,int) -> (int,int)]
        Optimal seam pair (b, c) that achieves M_v(u, w).
    left / right : OLOState | None
        Child states.
    """

    def __init__(self, leaves, left_leaves=None, right_leaves=None):
        self.leaves = list(leaves)
        self.left_leaves = (
            list(left_leaves) if left_leaves is not None else None
        )
        self.right_leaves = (
            list(right_leaves) if right_leaves is not None else None
        )
        self.cost = {}
        self.trace = {}
        self.left = None
        self.right = None


def bar_joseph_olo(tree, traits_df):
    start_time = time.time()

    # ------------------------------------------------------------------
    # Leaf indexing  (global IDs 0 … n-1 matched to trait-matrix rows)
    # ------------------------------------------------------------------
    terminals = tree.get_terminals()
    names = [t.name for t in terminals]
    leaf_id = {name: i for i, name in enumerate(names)}

    X = traits_df.loc[names].values
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)   # (n, n)

    # ------------------------------------------------------------------
    # Bottom-up DP
    # ------------------------------------------------------------------
    def DP(clade):
        if clade.is_terminal():
            idx = leaf_id[clade.name]
            state = OLOState([idx])
            state.cost[(idx, idx)] = 0.0          # M_tip(idx, idx) = 0
            return state

        if len(clade.clades) != 2:
            raise ValueError(
                "Tree must be fully bifurcating; resolve polytomies first."
            )

        left  = DP(clade.clades[0])
        right = DP(clade.clades[1])
        Lleaves = left.leaves      # global indices
        Rleaves = right.leaves
        nL = len(Lleaves)
        nR = len(Rleaves)

        # Build full symmetric M matrices for both children so that
        #   M_L[i, j] = M_L(Lleaves[i], Lleaves[j])
        # Entry (i, j) is finite only when i ≠ j (or the subtree is a tip).
        L_mat = _symmetric_cost_matrix(left)
        R_mat = _symmetric_cost_matrix(right)

        state = OLOState(Lleaves + Rleaves, Lleaves, Rleaves)
        state.left  = left
        state.right = right

        # ----------------------------------------------------------
        # Pass 1  –  g[b, w] = min_{c∈R} [ D(b,c) + M_R(c,w) ]
        #           + argmin index  g_arg[b, w]
        # ----------------------------------------------------------
        g     = np.full((nL, nR), np.inf)
        g_arg = np.zeros((nL, nR), dtype=int)

        for bi in range(nL):
            b_global = Lleaves[bi]
            D_b_R = D[b_global, Rleaves]          # D(b, c) for every c ∈ R

            for wi in range(nR):
                # R_mat[:, wi] = M_R(c, w) for every c ∈ R
                seam_costs = D_b_R + R_mat[:, wi]  # D(b,c) + M_R(c,w)
                ci = int(np.argmin(seam_costs))
                g[bi, wi]     = seam_costs[ci]
                g_arg[bi, wi] = ci

        # ----------------------------------------------------------
        # Pass 2  –  M_v(u,w) = min_{b∈L} [ M_L(u,b) + g(b,w) ]
        #           + seam pair (b*, c*) stored in state.trace
        # ----------------------------------------------------------
        for ui in range(nL):
            u = Lleaves[ui]
            M_L_u = L_mat[ui, :]          # M_L(u, b) for every b ∈ L

            for wi in range(nR):
                w = Rleaves[wi]
                totals = M_L_u + g[:, wi]  # M_L(u,b) + g(b,w) for every b ∈ L
                bi = int(np.argmin(totals))
                state.cost[(u, w)]  = float(totals[bi])
                state.trace[(u, w)] = (Lleaves[bi], Rleaves[g_arg[bi, wi]])

        return state

    root_state = DP(tree.root)

    # ------------------------------------------------------------------
    # Choose the globally optimal end-pair at the root
    # ------------------------------------------------------------------
    optimum  = np.inf
    end_pair = None
    for pair, value in root_state.cost.items():
        if value < optimum:
            optimum  = value
            end_pair = pair

    # ------------------------------------------------------------------
    # Traceback  –  recover the full leaf ordering from stored seam pairs
    # ------------------------------------------------------------------
    def traceback(state, a, b):
        """Return leaf-name list for subtree *state* with a at left end, b at right."""
        if len(state.leaves) == 1:
            return [names[state.leaves[0]]]

        left_set  = set(state.left_leaves)
        right_set = set(state.right_leaves)

        if a in left_set and b in right_set:
            x, y = state.trace[(a, b)]        # seam: x ∈ L, y ∈ R
            return (traceback(state.left,  a, x)
                    + traceback(state.right, y, b))

        if a in right_set and b in left_set:
            # M is symmetric; reversing the (b,a) layout gives the (a,b) layout.
            return traceback(state, b, a)[::-1]

        raise ValueError(f"Invalid boundary pair ({a}, {b})")

    order = traceback(root_state, end_pair[0], end_pair[1])
    reorder_tree_clades(tree, order)

    return tree, optimum, time.time() - start_time


def _symmetric_cost_matrix(state):
    """Build the full |leaves|×|leaves| M matrix for a subtree.

    Only cross-child entries (u∈L, w∈R) are stored in state.cost.
    The matrix is completed by the symmetry M(w,u) = M(u,w).
    Diagonal entries are 0 for tips, inf otherwise (u cannot equal w when
    the subtree has >1 leaf).
    """
    leaves = state.leaves
    n = len(leaves)
    M = np.full((n, n), np.inf)
    if n == 1:
        M[0, 0] = 0.0
        return M
    pos = {leaf: i for i, leaf in enumerate(leaves)}
    for (u, w), value in state.cost.items():
        i, j = pos[u], pos[w]
        M[i, j] = value
        M[j, i] = value
    return M


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
