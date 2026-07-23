import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

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
label_func = lambda c: (
    c.name if c.is_terminal() and c.name.endswith("0") else ""
)

# Replace plt.show() with:
fig.savefig("phylogeny_500tips.png", dpi=300, bbox_inches="tight")
# Or vector graphics:
fig.savefig("phylogeny_500tips.svg", format="svg")

# Color the first main branch red
my_tree.root.clades[0].color = "red"

plt.tight_layout()
plt.show()
