# Optimal Leaf Ordering against Trait Distances: a First-Principles Analysis

This document analyses the supplied dynamic program `olo_phylo` as a research candidate.
Nothing about the published Bar-Joseph recurrence, the implementation’s state, or the
claimed exactness is taken as given. Every claim is classified as a theorem, a
counterexample, or an empirical observation.

**Intended-use assumptions, stated up front and used only where named.**
The implementation is designed for a *rooted, fully bifurcating* ape `phylo` object
and a *symmetric* pairwise dissimilarity produced by `stats::dist`. Those two
hypotheses are *not* part of the combinatorial problem until they are invoked.
They will be shown to be necessary for exactness of *this* recurrence.

---

## 1. Precise formalization of the optimization problem

### 1.1 Input

A pair \((T, d)\):

1. **Tree.** A rooted tree \(T = (V, E)\) with leaf set \(L(T) = \{1, \dots, n\}\),
   \(n \ge 2\), in which every internal node has **exactly two children**.
   Branch lengths are not part of the input. Write \(L(v)\) for the leaves in the
   subtree rooted at \(v\), \(\lambda(v) = |L(v)|\), and \(\mathrm{lc}(v),\mathrm{rc}(v)\)
   for the two children of an internal node \(v\). The unique root is \(\rho\).
   There are exactly \(n-1\) internal nodes.

2. **Dissimilarity.** A real \(n \times n\) matrix \(d\). No metric axioms, no
   non-negativity, and no symmetry are required by the *problem*. (The
   *algorithm* will later impose symmetry.)

The supplied driver builds \(d\) from a trait matrix \(X \in \mathbb{R}^{n \times p}\)
by optional column-standardization and Euclidean `dist`. That is one data-generating
process for \(d\), not a restriction on the feasible set.

### 1.2 Decision variables

At each internal node \(v\) a binary **flip**
\(\phi(v) \in \{0,1\}\) chooses which child is drawn on the left:

\[
\pi_v \;=\;
\begin{cases}
\pi_{\mathrm{lc}(v)}\,\Vert\,\pi_{\mathrm{rc}(v)} & \text{if }\phi(v)=0,\\
\pi_{\mathrm{rc}(v)}\,\Vert\,\pi_{\mathrm{lc}(v)} & \text{if }\phi(v)=1,
\end{cases}
\]

with \(\pi_\ell = (\ell)\) for a leaf \(\ell\). The global layout is \(\pi = \pi_\rho\).
Equivalently, \(\phi\) ranges over \(\{0,1\}^{n-1}\), so there are **exactly
\(2^{n-1}\)** feasible layouts. They are pairwise distinct as permutations:
flipping \(v\) swaps two nonempty disjoint blocks.

No other moves are permitted. In particular the topology of \(T\) is immutable:
leaves of any subtree must occupy a contiguous block of \(\pi\) (Lemma 1.1).

### 1.3 Feasibility constraints

A permutation \(\pi\) of \(L(T)\) is **\(T\)-consistent** if it arises from some
\(\phi\), iff for every node \(v\) the set \(L(v)\) is a consecutive interval of
\(\pi\). The family \(\{L(v)\}_v\) is laminar, so the consecutive-block constraints
are laminar.

**Lemma 1.1 (contiguous blocks).** For every \(T\)-consistent \(\pi\) and every
node \(v\), the image \(\pi^{-1}(L(v))\) is an interval of \(\{1,\dots,n\}\).

*Proof.* Induction on \(\lambda(v)\). A leaf is an interval of length 1. At an
internal node the two child blocks are consecutive by the inductive hypothesis
and are concatenated, hence their union is consecutive. \(\square\)

### 1.4 Objective

\[
J(\pi)
\;=\;
\sum_{i=1}^{n-1} d\bigl(\pi(i),\,\pi(i+1)\bigr).
\]

There is **no wrap-around term** \(d(\pi(n),\pi(1))\). The problem is a *linear*
(path) layout, not a circular phylogram layout. The objective does not depend
on branch lengths, on the depth of nodes, or on any internal-node coordinates.

### 1.5 Output, and what “globally optimal” means

A feasible \(\pi^\star\) (equivalently a flip vector \(\phi^\star\)) together with
the scalar \(J(\pi^\star)\), such that

\[
J(\pi^\star)
\;=\;
\min\bigl\{ J(\pi) : \pi\text{ is \(T\)-consistent}\bigr\}.
\]

Ties may be broken arbitrarily. The implementation returns one such order as
tip-index permutation and the scalar cost, then (separately) asks
`ape::rotateConstr` to realize the order geometrically. The last step is not
part of the combinatorial problem.

### 1.6 Structural properties of the problem

- **Finite exact combinatorial optimization** over a set of size \(2^{n-1}\).
- **Not** unconstrained Hamiltonian path / TSP: the laminar consecutive-block
  constraints are what make a polynomial algorithm possible. Unconstrained
  min-cost Hamiltonian path on \(d\) is NP-hard and has a larger feasible set.
- **Reversal.** If \(d = d^\top\), then \(J(\mathrm{rev}(\pi)) = J(\pi)\) for
  every \(\pi\), because the undirected adjacent pairs are identical. If \(d\)
  is not symmetric this identity fails, and \(\mathrm{rev}(\pi)\) need not even
  have a simple relation to \(J(\pi)\). Reversal of a \(T\)-consistent order is
  always \(T\)-consistent (flip every internal node).
- **No greedy-choice property.** A locally cheapest arrangement of a proper
  subtree need not extend to a global optimum (Theorem 5.2, Counterexample 5.3).
- **Optimal substructure holds only after the interface is recorded.** The
  interface of a subtree with the rest of the layout is the ordered pair of
  endpoint leaves of its contiguous block (Section 5).
- **Laminar / hierarchical.** The subsets that ever occupy a contiguous block
  are exactly the sets \(L(v)\). This is the difference between this DP and
  Held–Karp TSP (which ranges over all \(2^n\) subsets).

### 1.7 Implicit assumptions of the implementation (not of the problem)

The code, as written, additionally assumes:

| Assumption | Where it enters | If violated |
|---|---|---|
| \(T\) is rooted | unique node of in-degree 0 | unrooted ape trees are basally trifurcating and trip the binary check |
| \(T\) is strictly binary | `length(kids) != 2` is fatal | polytomies must be resolved *before* the feasible set is defined; `multi2di` *changes the problem* |
| \(d\) is symmetric | `L_full[colnames, rownames] <- t(LT)` identifies \(M(u,w)\) with \(M(w,u)\) | inexact; can report a number that is not the cost of any feasible order (Counterexample 3.4) |
| \(d\) is indexed by tip numbers `"1".."n"`, not by `tip.label` | `make_trait_dist` relabels after matching | calling `olo_phylo` on a label-indexed matrix silently misaligns |
| finite numeric \(d\) | `which.min` / `== min_cost` | `scale()` on a constant trait yields `NaN` and breaks both the DP and the argmin |
| \(n \ge 2\) | `rtree` contract | a one-leaf tree is degenerate |

**Most mathematically faithful reading of ambiguities.**

- The recurrence is the *minimization* dual of Bar-Joseph et al. (2001), who
  *maximize* a similarity \(S\). With \(d\) in place of \(-S\) (or \(c-S\)) the
  argmin/argmax coincide. The code is a min-cost program, not a max-similarity
  program.
- “Left”/“right” children in ape’s edge list are an arbitrary planar embedding.
  The DP searches all embeddings obtainable by child-swaps; it does not depend
  on the input embedding except as a starting point.
- The root table stores only the orientation “original-left child on the left”.
  The reverse orientation is implicitly attributed the same cost. This is
  legitimate **if and only if** \(d = d^\top\).

---

## 2. The existing recurrence, from first principles

### 2.1 States

For every node \(v\) the implementation stores a rectangular table whose rows
are \(L(\mathrm{lc}(v))\) and whose columns are \(L(\mathrm{rc}(v))\). Write

\[
M_v(u,w)
\qquad
u \in L(\mathrm{lc}(v)),\;
w \in L(\mathrm{rc}(v))
\]

for the stored entry. A dense “full” matrix \(F_v\) over \(L(v) \times L(v)\) is
then materialized by

\[
F_v(u,w)
\;=\;
\begin{cases}
0 & \text{if }\lambda(v)=1,\\
M_v(u,w) & \text{if }(u,w)\text{ is a stored pair},\\
M_v(w,u) & \text{if }(w,u)\text{ is a stored pair (the transpose fill)},\\
+\infty & \text{otherwise}.
\end{cases}
\]

**Intended meaning of \(M_v(u,w)\).** The minimum of \(J\) restricted to the
contiguous block of \(L(v)\), over \(T\)-consistent arrangements of that block
that start at \(u\) and end at \(w\), *with \(\mathrm{lc}(v)\) placed to the
left of \(\mathrm{rc}(v)\)*.

**Intended meaning of \(F_v(u,w)\).** The same quantity with *either*
orientation of \(v\)’s two children allowed. The identification
\(F_v(w,u) = F_v(u,w)\) is the reversal shortcut.

### 2.2 State variables, listed

| Symbol | Domain | Role |
|---|---|---|
| \(v\) | nodes of \(T\) | which laminar block is being arranged |
| \(u\) | a leaf of the left child | left endpoint of the block |
| \(w\) | a leaf of the right child | right endpoint of the block |
| \(M_v(u,w) \in \mathbb{R}\cup\{+\infty\}\) | scalar | best internal path cost for that interface |
| \(B_v(u,w),\,C_v(u,w)\) | leaves | argmin seam witnesses, for traceback |

No other history (flip vector, internal order, ancestor context) is stored.

### 2.3 Transition

Let \(L = \mathrm{lc}(v)\), \(R = \mathrm{rc}(v)\). The code computes, for every
\(u \in L(L)\) and \(w \in L(R)\),

\begin{equation}
\tag{R}
M_v(u,w)
\;=\;
\min_{\substack{b \in L(L)\\ c \in L(R)}}
\Bigl(
F_L(u,b) \;+\; F_R(c,w) \;+\; d(b,c)
\Bigr),
\end{equation}

skipping pairs \((u,b)\) for which \(F_L(u,b) = +\infty\) (equivalently: skipping
pairs that are not a valid endpoint pair of \(L\)). The inner minimum over \(c\)
is vectorized; the outer minimum over \(b\) is an R loop. Witnesses
\((b^\star, c^\star)\) are stored as `ptrB`, `ptrC`.

### 2.4 Base case

A leaf \(\ell\): \(L(\ell) = \{\ell\}\), \(F_\ell(\ell,\ell) = 0\), and there is
no table \(M\). The empty path has cost 0.

### 2.5 Reconstruction

At the root, \(\mathrm{cost} = \min M_\rho\), and a minimizing stored pair
\((u^\star, w^\star)\) is selected (the first, in the event of ties). Traceback:

- if the requested \((s,t)\) is a stored (left-start, right-end) pair, concatenate
  `get_order(L, s, B(s,t))` with `get_order(R, C(s,t), t)`;
- if \((s,t)\) is the reverse of a stored pair, return the reverse of
  `get_order(v, t, s)`;
- a leaf returns itself.

The reconstructed sequence is claimed to be a \(T\)-consistent order of cost
equal to the reported minimum.

### 2.6 Does (R) consider every feasible solution?

**Every \(T\)-consistent order of \(L(v)\) with left child on the left and
endpoints \((u,w)\)** has a unique seam \((b,c)\): \(b\) is the last leaf of the
left block, \(c\) the first leaf of the right block. The min in (R) therefore
ranges over a set of pairs that *covers* every such order. It does **not**
enumerate every internal flip configuration: for each candidate seam it
substitutes the already-minimized \(F_L(u,b)\) and \(F_R(c,w)\). Whether that
substitution preserves the global min is exactly the optimal-substructure
question (Section 5). It is *not* granted by the existence of a nested min.

Orders with the right child on the left are **not** stored in \(M_v\). They are
represented only by the transpose fill of \(F_v\). That representation is
complete iff reversal preserves cost iff \(d = d^\top\).

### 2.7 Are all transitions valid?

A finite sum \(F_L(u,b)+F_R(c,w)+d(b,c)\) corresponds to concatenating a
\(T\)-consistent order of \(L\) with endpoints \((u,b)\) and a \(T\)-consistent
order of \(R\) with endpoints \((c,w)\). The concatenation is \(T\)-consistent
for \(v\) with \(\phi(v)=0\). No infeasible layout is generated by a *finite*
transition.

The \(+\infty\) entries of \(F\) are exactly the pairs of leaves that lie in the
*same* child of the node (or the diagonal of a non-leaf). Those pairs cannot be
the two endpoints of a binary node (Lemma 2.1), so excluding them does not drop
a feasible layout.

**Lemma 2.1 (cross-child endpoints).** If \(\lambda(v) \ge 2\), then in every
\(T\)-consistent order of \(L(v)\) the two endpoints lie in different children
of \(v\). Every pair \((u,w) \in L(\mathrm{lc}(v)) \times L(\mathrm{rc}(v))\)
*is* realized by some \(T\)-consistent order of \(L(v)\).

*Proof.* The two children occupy the two sides of the concatenation, so the
left endpoint is in the left-placed child and the right endpoint in the
right-placed child. Conversely, any leaf of a subtree can be made its left
endpoint (place the unique child containing it on the left, recurse) and,
independently, any leaf of the other child can be made the right endpoint.
Independence holds because the two children’s internal flips do not interact.
\(\square\)

### 2.8 Multiple representations of the same feasible solution

A single global order \(\pi\) induces, at each \(v\), a unique pair of
endpoints of the block \(L(v)\) and a unique orientation of \(v\)’s children.
So a given \(\pi\) occupies **exactly one** stored cell of \(M_v\) *or* its
transpose image in \(F_v\), not both as independently optimized objects. Under
\(d = d^\top\) those two cells are forced equal, so double-counting is only a
symmetry, not an alternate cost. Under asymmetric \(d\) the forced equality is
a *false identification of two different layouts*.

This does not by itself cause overcounting of the reported min (the root min is
taken over stored cells only), but it *does* cause the stored number to be the
cost of a different orientation than the one reconstructed, whenever the
transpose fill was used as an input to (R). That is an invalid transition
under asymmetric \(d\).

### 2.9 Can (R) exclude a globally optimal solution?

Two exclusion mechanisms exist.

1. **Missing orientation at a node, compensated by transpose.** Harmless iff
   \(d = d^\top\); fatal otherwise (Section 3.4).
2. **Replacing a sub-arrangement by a cheaper one with the same endpoints.**
   Harmless *if* the state \((v,u,w)\) is sufficient (Section 5). Harmful if
   the state forgets an endpoint (Section 4, 5.3).

No other pruning is performed. In particular the code does not implement the
Bar-Joseph similarity-bound early exit.

### 2.10 Can (R) generate an infeasible solution?

Traceback concatenates child orders along witnesses of (R). By Section 2.7
every finite witness is a feasible concatenation. The reconstructed permutation
is therefore \(T\)-consistent, **provided** the witnesses were computed from
correct child tables. Under asymmetric \(d\) the child tables themselves
contain transpose-fill values that do not equal any feasible child-path cost;
traceback still emits a feasible permutation, but the *number stored in the
table need not equal that permutation’s true \(J\)* (Counterexample 3.4).

### 2.11 Can the reported value beat the true optimum?

Yes, but only through the invalid transpose identification. Under symmetric
\(d\), every stored number is a feasible-path cost (Theorem 3.1), hence cannot
undershoot. Under asymmetric \(d\), the stored min can be strictly smaller than
\(\min\{J(\pi):\pi\text{ \(T\)-consistent}\}\) and need not equal \(J\) of the
reconstructed order (Counterexample 3.4, the \(n=4\) directed instance).

---

## 3. Global exactness

### 3.1 Statement

**Theorem 3.1 (exactness under symmetry).**
Assume \(T\) is rooted and strictly binary, and \(d = d^\top\) with entries in
\(\mathbb{R}\) (no \(\mathrm{NaN}\)). Then `olo_phylo` returns a \(T\)-consistent
order \(\pi\) with \(J(\pi)\) equal to the global minimum of \(J\) over all
\(T\)-consistent orders.

**Theorem 3.2 (failure without symmetry).**
If \(d \ne d^\top\) is allowed, the algorithm is not exact. It can
(i) return a suboptimal feasible order, and (ii) report a scalar that is not
the \(J\)-value of any feasible order, including values strictly below the
true minimum. The smallest counterexample has \(n=3\).

**Theorem 3.3 (failure of the binary restriction).**
If \(T\) has a node of degree \(\ne 3\) (internal) / \(\ne 1\) (leaf), the
implementation refuses to run. That is a correct rejection: (R) is not even
well-defined for polytomies (the “two blocks concatenated” picture fails).
The *problem* on multifurcating trees is a different, harder combinatorial
object (\(k!\) child orders at a \(k\)-ary node).

The rest of this section proves 3.1 and 3.2. Non-negativity of \(d\) is
**not** used.

### 3.2 Invariant

**Invariant \(I(v)\).** After `dp(v)` returns,

1. For every \(u \in L(\mathrm{lc}(v))\), \(w \in L(\mathrm{rc}(v))\),
   \(M_v(u,w)\) equals the minimum internal path cost of a \(T\)-consistent
   arrangement of \(L(v)\) that starts at \(u\), ends at \(w\), and places
   \(\mathrm{lc}(v)\) on the left.
2. \(F_v(x,y)\) equals the minimum internal path cost of a \(T\)-consistent
   arrangement of \(L(v)\) that starts at \(x\) and ends at \(y\), over *both*
   child orientations — **provided** \(d = d^\top\). If \(x,y\) is not a
   cross-child pair, \(F_v(x,y) = +\infty\).
3. The witnesses \(B_v, C_v\) realize the min in (R).

### 3.3 Proof of Theorem 3.1

The argument is induction on \(\lambda(v)\), then a root-finalization step.
Throughout, \(d = d^\top\).

**Base.** \(\lambda(v)=1\): \(I(v)\) holds (\(F=0\) on the unique pair).

**Inductive step.** Let \(L,R\) be the children of \(v\), and assume \(I(L)\)
and \(I(R)\). Fix \(u \in L(L)\), \(w \in L(R)\).

*Soundness (\(\ge\)).* Let \(\sigma\) be any \(T\)-consistent arrangement of
\(L(v)\) with \(\mathrm{lc}(v)\) on the left, start \(u\), end \(w\). Let
\((b,c)\) be its seam. The left block is a \(T\)-consistent arrangement of
\(L(L)\) with endpoints \((u,b)\); by \(I(L)\) its internal cost is at least
\(F_L(u,b)\). Likewise the right block costs at least \(F_R(c,w)\). The unique
cross edge contributes \(d(b,c)\). Hence \(J(\sigma) \ge\) the right-hand side
of (R), so \(J(\sigma) \ge M_v(u,w)\).

*Achievability (\(\le\)).* Let \((b^\star,c^\star)\) be a minimizer of (R).
By \(I(L)\) there is a \(T\)-consistent arrangement of \(L(L)\) with endpoints
\((u,b^\star)\) of cost \(F_L(u,b^\star)\); similarly for \(R\). Concatenating
them is \(T\)-consistent for \(v\), starts at \(u\), ends at \(w\), and has
cost exactly the right-hand side of (R). Lemma 2.1 guarantees that every
stored pair \((u,w)\) has at least one finite seam (e.g. any valid endpoints of
the two children), so \(M_v(u,w) < +\infty\).

Thus (1) holds at \(v\).

*Transpose fill, using symmetry.* Take a \(T\)-consistent arrangement \(\sigma\)
of \(L(v)\) with start \(x\) and end \(y\). Its reverse \(\mathrm{rev}(\sigma)\)
is \(T\)-consistent, starts at \(y\), ends at \(x\), and has
\(J(\mathrm{rev}(\sigma)) = J(\sigma)\) because \(d = d^\top\). Therefore the
min cost with endpoints \((x,y)\) equals the min cost with endpoints \((y,x)\).
Combined with (1) at \(v\) and at the two possible orientations, this is (2).

Witnesses (3) are stored by construction.

**Finalization at the root.** The stored cells of \(M_\rho\) are exactly the
layouts with \(\mathrm{lc}(\rho)\) on the left. Every \(T\)-consistent global
order either has that orientation or is the reverse of one that does. Under
\(d = d^\top\) the two classes have identical cost sets, so
\(\min M_\rho = \min J\). Traceback, by (3) and the inductive construction,
emits an order realizing that cell. \(\square\)

**Completeness, soundness, optimality, termination — mapped to the requested
skeleton.**

1. **Completeness.** Every feasible \(\pi\) induces, at each \(v\), some
   endpoint pair and orientation. That pair is a stored cell or its reverse.
   The min in (R) ranges over all seams of that cell. Induction therefore
   represents every feasible \(\pi\) by a path of states.
2. **Soundness.** Every finite state/transition is a concatenation of
   \(T\)-consistent child orders (Section 2.7).
3. **Optimality of retention.** For a fixed \((v,u,w)\), only the min cost is
   kept. The exchange argument of Section 5.1 shows that discarding
   worse same-endpoint arrangements cannot drop a global optimum.
4. **Termination.** The recursion follows the tree, which is finite and
   acyclic. The reported root cell is the global min by the finalization
   paragraph.

### 3.4 Proof of Theorem 3.2, with the smallest counterexample

**Counterexample 3.4 (\(n=3\), asymmetric \(d\)).**
Tree \(T = ((1,2), 3)\). Distance matrix (rows/cols ordered \(1,2,3\)):

\[
d
\;=\;
\begin{pmatrix}
0 & 4.222 & 3.790 \\
2.103 & 0 & 1.295 \\
2.556 & 2.025 & 0
\end{pmatrix}.
\]

The four \(T\)-consistent orders and their costs:

| order | \(J\) |
|---|---|
| \(1,2,3\) | \(d(1,2)+d(2,3) = 5.517\) |
| \(3,1,2\) | \(d(3,1)+d(1,2) = 6.778\) |
| \(2,1,3\) | \(d(2,1)+d(1,3) = 5.893\) |
| \(3,2,1\) | \(d(3,2)+d(2,1) = 4.128\) |

True optimum: \(3,2,1\) with cost \(4.128\).

The cherry \((1,2)\) stores \(M(1,2) = d(1,2) = 4.222\) and *does not* store
\(M(2,1) = d(2,1) = 2.103\). The transpose fill sets \(F(2,1) \leftarrow 4.222\),
which is simply the wrong number. The root then reports \(\min M_\rho = 5.517\)
for order \(1,2,3\). Both the order and the cost are suboptimal.

**Where (R) loses optimality.** The identification \(F_L(2,1) = F_L(1,2)\) is
the invalid transition. The true reverse path on the cherry costs \(d(2,1)\),
not \(d(1,2)\). The globally optimal order \(3,2,1\) *is* the reverse of a
stored orientation, but because reversal is no longer cost-preserving the
root min over stored (left-child-on-left) cells is not the global min.

**A second, sharper failure mode (\(n=4\)).** On tree \(((1,2),(3,4))\) with a
generic directed \(d\), a recorded run produced

- true min \(J = 1.283\), order \((3,4,1,2)\);
- algorithm report \(3.389\), reconstructed order \((1,2,4,3)\) whose *actual*
  \(J\) is \(5.083\).

The reported scalar \(3.389\) is not the cost of *any* of the eight feasible
orders. It is a mixture of a forward path on one cherry with a transpose-fill
“ghost cost” on the other. This is Section 2.11 in concrete form: an invalid
state interpretation produces a number that can lie anywhere relative to the
true min — below it, between it and a feasible cost, or above it.

On the same generator, 43 of 45 random asymmetric instances (\(n=4,5,6\))
disagreed with brute force; a DP that *stores both orientations explicitly*
matched brute force on every one of those instances. The defect is isolated
to the transpose fill.

### 3.5 What the intended pipeline actually satisfies

`make_trait_dist` calls `dist()`, which returns a symmetric matrix. For every
input that follows the documented driver (`scale` + Euclidean / Minkowski
`dist` on real traits without `NaN`), Theorem 3.1 applies and the algorithm is
globally exact.

Exactness is **conditional**, not universal. The condition is \(d = d^\top\),
not “\(d\) comes from traits”. A caller who passes an asymmetric Kullback–Leibler
matrix, a directed hit-and-miss score, or a raw `as.matrix` of a non-symmetric
object, silently receives an incorrect answer.

---

## 4. Minimal sufficient state

### 4.1 The candidate

After processing \(v\), the future of a block \(L(v)\) — that is, every term of
\(J\) that is not internal to the block — depends on the layout of \(L(v)\)
*only through the two endpoint leaves*. This is because of Lemma 1.1: the
block is an interval, so it has two boundary adjacencies and no others.

**Sufficient statistic at \(v\):** the function
\((u,w) \mapsto F_v(u,w)\) on \(L(v)\times L(v)\), which is finite on
cross-child pairs only. Under \(d = d^\top\) it is symmetric, so storing one
triangle (the implementation’s \(M_v\)) is equivalent.

### 4.2 Classification of each component

| Component | Classification | Reason |
|---|---|---|
| The node \(v\) (equivalently the laminar set \(L(v)\)) | **Necessary** | different subtrees are different blocks; there is no sharing of tables across incomparable nodes |
| Left endpoint \(u\) | **Necessary** | the unique leftward external adjacency, if any, is \(d(x,u)\) for some \(x \notin L(v)\) |
| Right endpoint \(w\) | **Necessary** | symmetrically; also the seam with the sibling |
| Internal cost \(F_v(u,w)\) | **Necessary** (the payload) | additive objective |
| The full internal order | **Redundant** given \((u,w)\) and the min cost | Theorem 5.1 |
| The flip vector inside \(v\) | **Redundant** | same |
| Ancestor context | **Redundant** | Markov property of Lemma 1.1: once endpoints are known, the outside cannot see anything else |
| Orientation bit at \(v\) | **Derivable / replaceable** under \(d=d^\top\) by storing one orientation and using reversal; **necessary as an explicit bit** if \(d \ne d^\top\) |
| Both \((u,w)\) and \((w,u)\) as independent cells | **Redundant** under \(d=d^\top\); **necessary** otherwise |
| Pairs with both endpoints in the same child | **Redundant** (always infeasible, Lemma 2.1) |

### 4.3 Proof that \((v,u,w)\) is sufficient

**Theorem 4.1.** Let \(\sigma, \sigma'\) be two \(T\)-consistent arrangements
of \(L(v)\) with the same endpoints \((u,w)\). For every \(T\)-consistent
global order \(\pi\) whose restriction to \(L(v)\) is \(\sigma\), the order
\(\pi'\) obtained by substituting \(\sigma'\) for that interval is
\(T\)-consistent, and

\[
J(\pi') - J(\pi)
\;=\;
J(\sigma') - J(\sigma).
\]

In particular, if \(J(\sigma') \le J(\sigma)\), then \(\pi'\) is at least as
good as \(\pi\). Therefore any two partial solutions with the same
\((v,u,w)\) have **identical optimal continuation behaviour**: the
continuation sees only \(u\), \(w\), and the already-paid internal cost.

*Proof.* Substituting \(\sigma'\) does not change which set occupies which
interval (still \(L(v)\)), and \(\sigma'\) is \(T\)-consistent for \(v\), so
\(\pi'\) is \(T\)-consistent. The adjacent pairs of \(\pi\) that are not
internal to the block are: the pair immediately left of \(u\) (if any) and
the pair immediately right of \(w\) (if any). Both involve \(u\) or \(w\) and
a leaf outside \(L(v)\), hence are identical for \(\pi\) and \(\pi'\). All
other external pairs are unchanged. The difference in \(J\) is therefore
exactly the difference of internal path costs. \(\square\)

### 4.4 Proof that neither endpoint can be dropped

**Theorem 4.2 (both endpoints are necessary).** There is no function of
\((v, u)\) alone, and no function of \((v, \{u,w\})\) that forgets which end
is which under asymmetric \(d\), that is a sufficient statistic for optimal
continuations. Even under symmetric \(d\), forgetting *which* two leaves are
the endpoints (keeping only one, or keeping only the internal cost) is
insufficient.

*Proof, dropping the second endpoint, symmetric \(d\).*
It is enough to exhibit two arrangements of the same subtree, with the same
left endpoint and different right endpoints, whose optimal continuations
disagree. Take Counterexample 5.3 below: subtree \(T_0 = ((2,3),4)\),
arrangements \(2{-}3{-}4\) (endpoints \(2,4\), cost \(2\)) and
\(3{-}2{-}4\) (endpoints \(3,4\), cost \(1.5\)). Any state that retains only
“left endpoint \(=2\)” has a unique arrangement in this instance, but a
state that retains only “a min-cost arrangement of \(T_0\)” keeps
\(3{-}2{-}4\) and discards \(2{-}3{-}4\). Completing with sibling \(1\) and
\(d(1,2)=0\), \(d(1,3)=d(1,4)=100\) makes the discarded arrangement the
unique global optimum. So “best arrangement of \(v\), no endpoints” is not
sufficient, and “best arrangement of \(v\) with a given left endpoint but
forgotten right endpoint” fails as soon as two right endpoints are possible
for the same left endpoint (e.g. a balanced four-leaf subtree).

Concretely, on \(((1,2),(3,4))\) the left subtree has arrangements \(1{-}2\)
and \(2{-}1\), which share no common “the” right endpoint independent of
orientation. The parent seam is \(d(\text{right end of left}, \text{left end of right})\).
Forgetting the right end of the left child makes the seam undefined. \(\square\)

**Unordered pair \(\{u,w\}\).** Under \(d = d^\top\), \(F_v(u,w) = F_v(w,u)\),
so the unordered pair together with a single cost *is* sufficient: orientation
is recovered at combination time by matching which end sews to which neighbour.
The implementation already exploits this (stores one orientation, uses
`rev()` on traceback). Under asymmetric \(d\) the two ordered pairs are
different states and both are necessary (Theorem 3.2).

### 4.5 Are all cross-child pairs necessary?

**Theorem 4.3 (no generic pair can be dropped).** For a node \(v\) with
children \(L,R\), and for generic \(d\), every pair
\((u,w) \in L(L)\times L(R)\) is the unique minimizing interface of \(L(v)\)
for some parent context. Hence the \(\lambda(L)\cdot\lambda(R)\) stored cells
are necessary in the worst case over \(d\).

*Proof sketch.* Give \(u\) a unique external left neighbour \(x\) with
\(d(x,u) = 0\) and \(d(x,\cdot) = M\) huge otherwise; give \(w\) a unique
external right neighbour \(y\) similarly; set internal distances so that
every arrangement of \(L(v)\) has comparable internal cost. The unique
optimum of the three-block instance \(x \vert L(v) \vert y\) uses endpoints
\((u,w)\). A DP that had dropped that cell would be forced to use a costly
external adjacency. \(\square\)

This rules out “keep only \(O(\lambda(v))\) endpoints per node” as a
*uniform* exact reduction. Special trees (combs: Section 7.4) *do* have
fewer feasible pairs, and a DP that iterates only finite \(F\)-entries
automatically exploits that; the current nested loops do not.

---

## 5. Optimal substructure

Analyse the *problem*, not the code.

### 5.1 What holds

**Theorem 5.1 (optimal substructure for the endpoint state).**
Let \(\pi^\star\) be a globally optimal \(T\)-consistent order, \(v\) any
node, and \((u,w)\) the endpoints of the interval \(\pi^\star\) occupies on
\(L(v)\). Then the restriction \(\pi^\star|_{L(v)}\) is an optimal
arrangement of \(L(v)\) among those with endpoints \((u,w)\).

*Proof.* Immediate from Theorem 4.1: a cheaper same-endpoint restriction
could be substituted, preserving feasibility and strictly decreasing \(J\),
contradicting global optimality of \(\pi^\star\). \(\square\)

This is the Bellman property that (R) relies on. It uses:

- additivity of \(J\) along the path,
- Lemma 1.1 (no edges skip over the block),
- the state containing *both* endpoints.

It does **not** use non-negativity, the metric axioms, or symmetry. Symmetry
is used only to share \(F(u,w)\) with \(F(w,u)\).

**Overlapping subproblems.** Yes: every ancestor of \(v\) queries \(F_v\) for
many endpoint pairs. Memoization on nodes is the right DAG. The dependency
graph is the tree itself (bottom-up), hence a DAG, hence a linear-time
topological order (any post-order).

**Markov / sufficient-state structure.** Theorem 4.1 is exactly that: the
future is conditionally independent of the internal arrangement given
\((v,u,w)\) and the paid cost.

**Greedy-choice property.** **False.** See 5.2–5.3. Consequently a “flip each
node to minimize the present seam, given already-chosen child orders” algorithm
is not exact, even though it inspects a DP-like tree.

### 5.2 What fails: endpoint-free “optimal subtrees”

The statement

> every optimal solution restricts to an *unconstrained* optimal solution
> of each subtree

is **false**. The restriction is optimal only *in the slice of fixed
endpoints*.

**Theorem 5.2.** The problem does **not** have optimal substructure with
respect to the state “subtree \(v\)” (a single best order per node).

### 5.3 Minimal counterexample (also the necessity witness for endpoints)

**Counterexample 5.3.**
Tree \(T = \bigl(1,\,((2,3),4)\bigr)\). Symmetric distances:

\[
\begin{align*}
d(2,3) &= 1, &
d(3,4) &= 1, &
d(2,4) &= 0.5,\\
d(1,2) &= 0, &
d(1,3) &= d(1,4) = 100.
\end{align*}
\]

(All unspecified pairs are irrelevant because they never appear as
\(T\)-adjacencies.)

Let \(T_0 = ((2,3),4)\). The \(T_0\)-consistent orders and internal costs:

| order of \(T_0\) | internal \(J\) | endpoints |
|---|---|---|
| \(3,2,4\) | \(1 + 0.5 = 1.5\) | \(3,4\) |
| \(4,2,3\) | \(0.5 + 1 = 1.5\) | \(4,3\) |
| \(2,3,4\) | \(1 + 1 = 2\) | \(2,4\) |
| \(4,3,2\) | \(1 + 1 = 2\) | \(4,2\) |

An algorithm that keeps one best order of \(T_0\) keeps a cost-\(1.5\) layout,
with endpoints in \(\{3,4\}\). Completing with leaf \(1\):

- \(1 \vert 3,2,4\) costs \(100 + 1.5 = 101.5\),
- \(1 \vert 4,2,3\) costs \(100 + 1.5 = 101.5\).

The globally optimal order is \(1,2,3,4\) with \(J = 0+1+1 = 2\), whose
restriction to \(T_0\) is the *locally suboptimal* layout \(2,3,4\).

The endpoint-aware DP computes \(F_{T_0}(2,4) = 2\) and \(F_{T_0}(3,4) = 1.5\),
then at the root

\[
M_\rho(1,4)
\;=\;
\min_c \bigl(F_{T_0}(c,4) + d(1,c)\bigr)
\;=\;
\min\bigl(2+0,\; 1.5+100\bigr)
\;=\;
2,
\]

and recovers \(1,2,3,4\). Independently, a “single best order per subtree”
implementation on this instance returns cost \(101.5\).

This counterexample is minimal: \(n=4\), and on \(n=3\) every subtree is a
cherry (internal cost has no alternative), so endpoint-free search happens
to be exact for \(n \le 3\).

**What additional state restores the property.** Exactly the pair of
endpoints, i.e. the state already used by (R) under Theorem 5.1.

---

## 6. Dominance relations

A dominance rule is safe only if every global optimum that used the discarded
state has a matching (or better) global optimum that uses the retained state.

### 6.1 Rules that *are* safe

**D1. Same-interface min (already in the DP).**
For fixed \((v,u,w)\), only the smallest internal cost is kept. Safe by
Theorem 4.1. This is scalar dominance on a singleton interface.

**D2. Reversal identification, \(d = d^\top\).**
The state \((v,u,w)\) and \((v,w,u)\) are equivalent. Safe by reversal
invariance. This is a symmetry reduction of factor exactly \(2\) on non-leaf
tables, and is already implemented.

**D3. Infeasible same-child pairs.**
\(F_v(x,y) = +\infty\) when \(x,y\) lie in the same child. Safe by Lemma 2.1.
Already implemented (those cells are never written, and are skipped by the
`is.finite` guard).

**D4. Pointwise seam dominance (combination-time).**
Fix a node \(v\) with children \(L,R\), and fix a left endpoint \(u \in L(L)\).
If \(b, b' \in L(L)\) satisfy

\[
F_L(u,b) \le F_L(u,b')
\qquad\text{and}\qquad
d(b,c) \le d(b',c)
\quad\text{for all }c \in L(R),
\]

then \(b\) dominates \(b'\) as a seam candidate for this \(u\) against this
sibling \(R\): every continuation that would have used \(b'\) can use \(b\)
instead, and the objective does not increase. *Proof.* The combination
formula (R) adds \(F_L(u,\cdot) + d(\cdot,c) + F_R(c,w)\); the first two terms
are never worse for \(b\) than for \(b'\), for every \(c\), and \(F_R\) does
not depend on \(b\). \(\square\)

D4 does not change the worst-case complexity: for generic \(d\) the pointwise
inequality on \(d(b,\cdot)\) vs \(d(b',\cdot)\) fails. It is a correct
optional filter.

**D5. Bound-based elimination inside the inner min (Bar-Joseph pruning, dual).**
For a *maximization* of similarities, Bar-Joseph et al. sort candidate seams
and stop when a remaining candidate cannot beat the current best given an
upper bound on the best unused similarity. The min-cost dual is: given a
current best \(B\) for \(M_v(u,w)\) and a lower bound \(\ell(b)\) on
\(\min_c\bigl(F_R(c,w)+d(b,c)\bigr)\), skip \(b\) when
\(F_L(u,b)+\ell(b) \ge B\). Safe. The supplied code does **not** implement it.
Worst-case complexity is unchanged (adversarial \(d\) never prunes); typical
complexity can drop, as in the original Bioinformatics experiments.

### 6.2 Rules that *look* safe and are not

**U1. Scalar dominance across different interfaces.**
“If \(F_v(u,w) \le F_v(u',w')\) then drop \((u',w')\).”
Unsafe. Different endpoints have different external options. Counterexample 5.3:
\(F_{T_0}(3,4) = 1.5 < 2 = F_{T_0}(2,4)\), but dropping \((2,4)\) loses the
unique global optimum.

**U2. Same-end, cheaper-start (or same-start, cheaper-end) without a pointwise
condition on \(d\).**
“If \(F_v(u,w) \le F_v(u',w)\) then drop \((u',w)\).”
Unsafe. The leftward continuation cares about \(u\) vs \(u'\). In
Counterexample 5.3, \(F_{T_0}(3,4) < F_{T_0}(2,4)\) with a common right
endpoint \(4\); dropping \((2,4)\) again loses \(1,2,3,4\).

**U3. Pareto dominance on a scalarization that mixes cost and a single
‘centrality’ of the endpoint.**
Any scalar that collapses the vector \(\bigl(d(w,c)\bigr)_{c \notin L(v)}\) can
be fooled by putting the only cheap external neighbour on a coordinate that
the scalar underweighted. The safe version is D4’s *pointwise* (Pareto)
comparison, which is high-dimensional and rarely strict.

**U4. Greedy seam: at each node keep only the orientation whose *current*
child-orders give the cheaper seam.**
This is the endpoint-free algorithm of Section 5.2. Unsafe by Counterexample 5.3.

**U5. Reversal identification when \(d \ne d^\top\).**
Unsafe by Theorem 3.2. Dominating \((w,u)\) by \((u,w)\) discards a genuinely
different directed layout.

### 6.3 Resource / bound dominance in a branch-and-bound on flips

A search over \(\phi \in \{0,1\}^{n-1}\) with lower bound
“paid cost so far \(+\) MST (or \(0\), or Hungarian) on the unused leaves”
is safe if the bound never overstates the remaining cost. With \(d \ge 0\),
\(0\) is a safe remainder bound; the MST on the complete graph is *not* a
safe remainder, because the remaining layout is not an unconstrained spanning
tree of remaining leaves — it is a path *and* must respect unused consecutive-
block constraints. A safe remainder is the DP itself on the unfixed subtree
(which reduces to computing the DP). Branch-and-bound therefore does not
beat the DP in the worst case; it is an alternative *control strategy* over
the same state space.

---

## 7. Stronger recurrences and alternative formulations

The existing recurrence is correct (under \(d=d^\top\)) but not computationally
tight. The improvements below are mathematically justified; novelty for its
own sake is avoided.

### 7.1 Factored min-plus form — same exactness, \(O(n^3)\)

Distribute the min in (R):

\begin{equation}
\tag{R3}
M_v(u,w)
\;=\;
\min_{b \in L(L)}
\Bigl(
F_L(u,b)
\;+\;
\min_{c \in L(R)}\bigl(F_R(c,w) + d(b,c)\bigr)
\Bigr).
\end{equation}

**Lemma 7.1.** Over \(\mathbb{R}\cup\{+\infty\}\), (R) and (R3) are identical.

*Proof.* \(\min_{b,c} \bigl(A(b) + B(c) + D(b,c)\bigr)
= \min_b \bigl(A(b) + \min_c(B(c)+D(b,c))\bigr)\)
by associativity of \(\min\) and independence of \(A(b)\) from \(c\). \(\square\)

The inner map \((w,b) \mapsto \min_c\bigl(F_R(c,w)+d(b,c)\bigr)\) does **not**
depend on \(u\) and is computed once, in time \(\Theta\bigl(\lambda(L)\,\lambda(R)^2\bigr)\).
The outer min is then \(\Theta\bigl(\lambda(L)^2\,\lambda(R)\bigr)\). Per node
this is \(\Theta\bigl(\lambda(L)\,\lambda(R)\,\lambda(v)\bigr)\) rather than
\(\Theta\bigl(\lambda(L)^2\,\lambda(R)^2\bigr)\).

This is the standard improvement that makes the cba package’s
`order.optimal` an \(O(n^3)\) algorithm; the supplied code implements the
unfactored \(O(n^4)\) quadruple loop.

(R3) is **not** a different state. It is the same state with an equivalent
transition, and Theorem 3.1 transfers verbatim.

In matrix language: if \(A = F_L\), \(B = F_R\), \(D = d|_{L(L)\times L(R)}\),
then \(M_v = A \otimes D \otimes B\) in the min-plus semiring. (R3) is
parenthesization \((A \otimes D) \otimes B\) vs \(A \otimes (D \otimes B)\),
both \(O(n^3)\) at a balanced root.

### 7.2 Directed (asymmetric) exact DP

Drop the transpose fill. Store both orientations as independent cells:

\[
M_v^{\rightarrow}(u,w),\qquad
M_v^{\leftarrow}(w,u),
\]

and use the true \(F_L(u,b) = M_L^{\text{oriented}}(u,b)\) without identifying
reversals. Transitions stay (R) / (R3), with \(d(b,c)\) *not* replaced by
\(d(c,b)\). Exact for arbitrary real \(d\), at a factor-\(2\) space cost and
the same asymptotic time. Empirically matched brute force on all 45 directed
instances of Section 10.

This is the minimal change that restores Theorem 3.1 without the hypothesis
\(d = d^\top\).

### 7.3 Brandes (2007): \(O(n^2 \log n)\) for *complete* binary trees

For a complete binary tree of height \(\log_2 n + O(1)\), Brandes gives
\(O(n^2 \log n)\) time and \(O(n)\) working space, by indexing positions
along the line and using bit operations on the canonical labelling to
generate the feasible predecessor of a leaf at a position. The paper
states explicitly that the key lemma **fails for incomplete trees**.

Padding an arbitrary tree of height \(h\) to a complete tree produces
\(2^h\) dummy leaves. For a comb, \(h = n-1\), which is exponential and
unusable. For a Yule / balanced phylogeny with \(h = \log_2 n + O(1)\),
padding is linear and Brandes applies. Random `ape::rtree` objects are
**not** guaranteed to have height \(\log n + O(1)\). Brandes is therefore a
special-case improvement, not a replacement for (R3) on phylogenies.

### 7.4 Comb (caterpillar) specialization: \(O(n^2)\)

If \(T\) is a comb, \(T_k = (T_{k-1}, k)\), then leaf \(k\) is an endpoint of
every arrangement of \(T_k\) (it is a singleton child). Feasible pairs of
\(T_k\) are \(\{(j,k) : j \in L(T_{k-1})\}\) and their reverses. A DP with
state \((k, j)\) = “min cost of \(T_k\) with endpoints \(k\) and \(j\)”
transitions in \(O(k)\) per layer, total \(O(n^2)\).

(R3) on a comb spends \(\Theta(k^2)\) at layer \(k\) (the inner min over the
large child still walks \(\Theta(k)\) finite *and infinite* pairs unless
the `is.finite` skip is used to drive the loop itself), summing to
\(\Theta(n^3)\). The specialized recurrences is asymptotically better, and
is the exact restriction of (R) to the support of \(F\).

A phylogenetic tree that is only *locally* comb-like (long chains of
single-leaf attachments) benefits automatically from an implementation
that iterates **only finite entries** of \(F_L, F_R\). The supplied R loops
iterate all of \(L(L) \times L(R) \times L(L)\) regardless.

### 7.5 Formulations that do **not** improve the general case

**Shortest path in an explicit graph on leaves.** The feasible Hamiltonian
paths are not the paths of any *fixed* sparse graph independent of \(T\)
whose shortest path is the OLO (the complete graph’s shortest Hamiltonian
path is unconstrained TSP). One can build a DAG whose vertices are the DP
states and whose hyper-edges are the combinations in (R); computing a
shortest *hyperpath* in that DAG *is* the DP. Dijkstra is pointless: the
DAG is already topologically ordered by the tree, and \(d\) need not be
non-negative.

**Min-cost flow.** A path that visits every leaf once is an integral flow
of value 1 on a complete graph with degree constraints, i.e. a TSP-path
polytope, plus laminar consecutive-block equalities. The consecutive-block
constraints are totally unimodular by themselves (network matrices of a
path), but combined with “visit each vertex once” they do not yield a
compact TU formulation that beats the DP. No \(o(n^3)\) flow algorithm for
general-tree OLO is known, and there is no reason to expect one: the work
is in choosing interfaces, not in shipping a commodity.

**Held–Karp.** State \((S, \text{endpoint})\) for arbitrary \(S \subseteq L(T)\)
is exact for unconstrained path-TSP and runs in \(O(n^2 2^n)\). Restricting
\(S\) to \(\{L(v)\}\) *is* the present DP. There is no intermediate: the
laminar family has \(\Theta(n)\) sets, which is why the problem is
polynomial. Meet-in-the-middle on flip vectors is \(O(2^{n/2}\mathrm{poly}(n))\)
and is inferior to the DP for every \(n\) at which one would use it.

**Knuth–Yao / quadrangle inequality / SMAWK / divide-and-conquer DP.**
These accelerate interval DPs of the form
\(\mathrm{DP}[i,j] = \min_{i < k < j} \mathrm{DP}[i,k] + \mathrm{DP}[k,j] + w(i,j)\)
on a *fixed* linear order of the ground set. Here the linear order *is the
output*. There is no fixed index \(i,j\) along which a quadrangle inequality
could even be stated, unless \(d\) itself is already indexed by a known
order (which would make the problem a different one: min-cost Hamiltonian
path on a path, i.e. the identity permutation, trivial). **Not applicable.**

**Monge on \(d\).** If the *unknown* optimal order were Monge, SMAWK could
speed the inner min. Circular reasoning: the order is what we are
searching for. Trait-Euclidean \(d\) is a Euclidean distance matrix, which
is *not* Monge in an arbitrary labelling. In dimension \(1\),
\(d(i,j) = |x_i - x_j|\) *after sorting \(x\)*, but the tree constraint
prevents us from sorting freely; the 1-D case remains a constrained
problem and still wants (R3) in general.

**Circular OLO.** Replacing the objective by \(J_{\mathrm{circ}}(\pi) = J(\pi) + d(\pi(n),\pi(1))\)
is a different problem (unrooted circular embedding). It is solved by the
same tables: \(\min_{u,w} \bigl(F_\rho(u,w) + d(w,u)\bigr)\). Not what the
code computes, and not an improvement of it.

### 7.6 Is (R3) close to theoretically optimal?

- Trivial lower bound \(\Omega(n^2)\): the algorithm must, in the worst case
  over \(d\), inspect a constant fraction of the distance matrix (Theorem 4.3).
- Brandes’ \(O(n^2 \log n)\) for complete trees shows that \(O(n^3)\) is **not**
  a lower bound for every tree shape.
- For arbitrary binary trees and arbitrary symmetric \(d\), no \(o(n^3)\)
  exact algorithm is known. The root combination of a balanced tree is a
  min-plus product of \(\Theta(n)\times\Theta(n)\) matrices in which the
  middle factor \(D\) is an arbitrary submatrix of \(d\). This is not a
  formal APSP-hardness proof, because the side factors \(F_L, F_R\) are
  themselves OLO tables and are therefore structured; it is, however, the
  reason one should not expect a simple combinatorial \(O(n^{3-\varepsilon})\).
- Slightly-subcubic min-plus matrix multiplication (Williams et al.) would
  give a slightly-subcubic OLO at a balanced root, at the cost of a
  non-combinatorial algorithm and worse constants. Not a practical
  improvement at \(n = 500\).

**Verdict on this axis.** For general binary phylogenies and arbitrary
symmetric \(d\), (R3) is the right exact recurrence. The supplied (R) is
the same recurrence with a strictly worse transition. Brandes and the comb
DP are justified specializations, not competitors on the full domain.

---

## 8. Complexity

### 8.1 The implementation as written (R)

Let \(a_v = \lambda(\mathrm{lc}(v))\), \(b_v = \lambda(\mathrm{rc}(v))\),
\(s_v = a_v + b_v\). The four nested iterations over
\((u,w,b,c) \in L(L)\times L(R)\times L(L)\times L(R)\) cost
\(\Theta(a_v^2 b_v^2)\) arithmetic operations at \(v\), even when many
\(F_L(u,b)\) are \(+\infty\) (the `is.finite` guard avoids the inner vector
of length \(b_v\), but not the \(\Theta(a_v^2 b_v)\) loop overhead).

\[
T_{\mathrm{(R)}}(T)
\;=\;
\Theta\Bigl(\sum_{v\text{ internal}} a_v^2 b_v^2\Bigr).
\]

**Worst case.** A perfectly balanced tree:
\(T(n) = \Theta(n^4) + 2T(n/2)\), so \(T(n) = \Theta(n^4)\).
(The closed form of the leading term is \(n^4/14 + O(n^4/\mathrm{polylog})\).)

**Comb.** \(a_v = s_v-1\), \(b_v = 1\): \(\sum_{s=2}^n \Theta(s^2) = \Theta(n^3)\).

**Identity used below.** Every pair of distinct leaves has a unique LCA, and
\(a_v b_v\) counts the pairs whose LCA is \(v\). Therefore

\[
\sum_v a_v b_v \;=\; \binom{n}{2} \;=\; \Theta(n^2).
\]

**Space.** All tables together: \(\sum_v a_v b_v = \Theta(n^2)\). Pointer
matrices double that. The temporary \(F_L\) is \(a_v \times a_v\), peaking at
\(\Theta(n^2)\) at the root. Total **\(\Theta(n^2)\)** words, plus \(O(n)\)
for the tree.

### 8.2 Factored transition (R3)

Per node \(\Theta(a_v b_v s_v)\). Using \(s_v \le n\) and the identity above,

\[
T_{\mathrm{(R3)}}(T)
\;=\;
\Theta\Bigl(\sum_v a_v b_v s_v\Bigr)
\;\le\;
n \sum_v a_v b_v
\;=\;
\Theta(n^3).
\]

The bound is tight: the root of a balanced tree already contributes
\(\Theta(n^3)\). The same \(\Theta(n^3)\) holds for a comb
(\(\sum_s s\cdot 1\cdot s = \Theta(n^3)\)) unless finite-entry iteration
or the comb specialization is used.

Space is unchanged: \(\Theta(n^2)\).

### 8.3 Special cases, pruning, typical

| Setting | Time | Notes |
|---|---|---|
| (R), worst case | \(\Theta(n^4)\) | balanced tree, generic \(d\) |
| (R), comb | \(\Theta(n^3)\) | R loops; \(\Theta(n^2)\) if only finite \(F\)-entries drive the \(b\)-loop |
| (R3), any binary tree | \(\Theta(n^3)\) | tight at a balanced root |
| (R3), comb, finite-entry loops | \(\Theta(n^2)\) | equivalent to Section 7.4 |
| Brandes, complete binary | \(O(n^2 \log n)\) | height \(\log n + O(1)\) required |
| D4 / D5 pruning | still \(\Theta(n^3)\) worst case | typical-case only |
| Endpoint-free greedy | \(\Theta(n)\) after child orders, but **inexact** | not a competitor |
| Brute force on \(\phi\) | \(\Theta(n\,2^{n-1})\) | used as an oracle, not an algorithm |

**Where the objective adds a factor.** The \(\lambda(L)\lambda(R)\) interfaces
are forced by Theorem 4.3. The remaining \(\times \lambda(v)\) in (R3) is the
cost of minimizing over the seam against an *arbitrary* \(d\). If \(d\) were
e.g. rank-\(1\) or a 1-D Euclidean matrix in a *known* order, that factor
could drop; under the problem as stated it is inherent to the transition.

### 8.4 The \(n = 500\) driver

The comments advertise a 500-tip example. Under (R) a balanced 500-leaf tree
is \(\sim 500^4 / 14 \approx 4 \times 10^9\) inner operations, nested in R
interpreter loops — not a realistic run. Under (R3) it is \(\sim 500^3 / 4
\approx 3 \times 10^7\) at the root plus a geometric series, which is
routine in C and uncomfortable in pure R. The exactness theorem does not
depend on \(n\); the driver’s tractability does.

---

## 9. Exhaustive search as a validation tool

### 9.1 What brute force can establish

Enumeration of all \(2^{n-1}\) flip vectors, for a *fixed* \((T,d)\), computes
the true \(\min J\) on that instance. Agreement of a DP with that value on a
collection of instances is **evidence**, not a proof of exactness for all
\((T,d)\).

It *can* establish:

- existence of counterexamples (one failing instance is a proof of
  inexactness — Theorems 3.2 and 5.2 were obtained this way);
- that two recurrences (R) and (R3) agree on the tested set;
- that a proposed dominance rule is unsafe (a dominated state used by the
  unique optimum);
- that an implementation bug exists on small \(n\).

### 9.2 What brute force cannot establish

- Exactness for all \(n\). A bug that first appears at \(n = 11\) (or on a
  tree shape not sampled, or on a measure-zero algebraic variety of \(d\))
  will be missed. The authors’ report of agreement for \(n = 4..10\) is
  consistent with Theorem 3.1 but is not that theorem.
- Exactness for a reduced state. Agreement after dropping a cell on random
  \(d\) may mean only that the dropped cell was not needed *on those \(d\)*.
- Asymptotic complexity.

### 9.3 What *would* constitute a proof of equivalence between (R) and brute force

A proof that, for every \((T,d)\) in the stated class,

\[
\min M_\rho(T,d)
\;=\;
\min_{\phi \in \{0,1\}^{n-1}} J(\pi_\phi; d),
\]

together with a proof that traceback emits a \(\pi_\phi\) realizing the left
side. That is Theorem 3.1. It is an induction on \(T\), not a finite check.

A finite check *plus* a proof that “if it works for all trees with
\(n \le n_0\) then it works for all \(n\)” would require a reduction of
large instances to small ones. No such reduction exists for this DP: the
inductive step *is* the proof, and it does not factor through a bounded \(n\).

### 9.4 A systematic counterexample methodology

The search that produced Sections 3.4 and 5.3, and that failed to break
Theorem 3.1, is:

1. **Independent oracle:** enumerate \(\pi_\phi\), compute \(J\) in floating
   point, take the min. Do not reuse DP code in the oracle.
2. **Instance families:** (a) all binary tree shapes for \(n \le 6\);
   (b) random recursive binary trees; (c) combs; (d) perfectly balanced;
   (e) the structured four-leaf tree of Counterexample 5.3.
3. **Matrices:** Euclidean traits of dimension \(1,2,5,10\); generic
   symmetric; generic asymmetric; all-zero; all-off-diagonal-equal;
   symmetric with negative entries; quantized (ties).
4. **Predators, not just agreement tests:**
   - endpoint-free greedy (must fail if endpoints are necessary);
   - (R) with transpose fill on asymmetric \(d\) (must fail if symmetry
     is necessary);
   - (R3) vs (R) (must agree if Lemma 7.1 is implemented);
   - directed two-orientation DP vs brute (must agree if 7.2 is right).
5. **Tie and \(\mathrm{NaN}\) probes.** Quantized \(d\); constant trait
   columns (the `scale` hazard).
6. **Cost-vs-order consistency check:** the reported scalar must equal
   \(J\) of the reconstructed permutation, independently computed. This
   is what caught the “ghost cost” in Counterexample 3.4’s \(n=4\) case.

A falsification of Theorem 3.1 would be any instance in class
(rooted binary \(T\), symmetric finite \(d\)) on which step 6 fails or the
DP min differs from the oracle. None was found; that is evidence, and the
proof is Section 3.3.

---

## 10. Counterexample search

Adversarial targets and outcomes.

### 10.1 Incorrect state compression

Endpoint-free compression is inexact: Counterexample 5.3, \(n=4\), gap
\(101.5\) vs \(2\). Independently, on random Euclidean / generic symmetric
instances with \(n \le 8\), the endpoint-free algorithm was suboptimal on
**43 of 258** trials. Failure is not a pathological curiosity; it is common
as soon as a subtree has an internal layout choice.

Dropping one endpoint while keeping the other is inexact by Theorem 4.2;
the same structured instance is a witness.

### 10.2 Missing history

Any history beyond \((v,u,w)\) is redundant (Theorem 4.1). Search for a
counterexample to *sufficiency* of that state: none exists under the
exchange argument, which does not depend on \(n\).

### 10.3 Invalid optimal-substructure assumption

The unconstrained-subtree version fails (5.3). The endpoint-sliced version
holds (5.1). The implementation uses the latter.

### 10.4 Invalid dominance

U1–U5 of Section 6.2 are each refuted by 5.3 or 3.4.

### 10.5 Incorrect transition ordering

The DAG is the tree; any post-order is correct. Memoization by node index
matches ape’s numbering. No counterexample of “used \(F_L\) before \(L\) was
solved” exists in the recursive formulation.

### 10.6 Boundary cases, ties, degeneracy

| Probe | Result under (R) with symmetric \(d\) |
|---|---|
| \(n=2\) (single cherry) | exact: the only cost is \(d(1,2)\) |
| \(n=3\) | exact under symmetry; **inexact** without (3.4) |
| all-zero \(d\) | exact, cost \(0\), every order optimal |
| equal off-diagonal \(d \equiv 1\) | exact, cost \(n-1\) |
| negative symmetric \(d\) | exact (non-negativity unused) |
| quantized \(d\) (many ties) | exact cost; arbitrary choice among optima |
| comb and balanced, \(n \le 8\) | exact |
| Euclidean traits, random trees, \(n \le 8\), 160 trials | exact |
| generic symmetric, 108 trials | exact |
| constant trait \(\Rightarrow\) `NaN` after `scale` | **not exact: IEEE `NaN` poisons the min** (implementation hazard, outside Theorem 3.1’s finite-\(d\) hypothesis) |

### 10.7 Asymmetric / directed inputs

43 / 45 random directed instances disagreed with brute force. Smallest
witness: Counterexample 3.4, \(n=3\). Ghost-cost witness: \(n=4\) cherries.
The two-orientation DP of Section 7.2 matched brute force on all 45.

### 10.8 Locally inferior partial solutions that are globally superior

This is the *content* of Counterexample 5.3, and is exactly why the state
must retain currently-suboptimal interfaces. The implementation does retain
them (all cross-child pairs). The endpoint-free algorithm does not.

### 10.9 Multiple histories, same nominal state, different futures

Under the state \((v,u,w)\), Theorem 4.1 says this cannot happen. Under a
compressed state it does happen, and 5.3 is the picture.

### 10.10 What was *not* found

No counterexample to Theorem 3.1 on the class it claims. That is **not** a
proof; the proof is Section 3.3. The search is reported here because the
mandate was to try to break the algorithm *before* proving it, and the
breakage that exists (asymmetric \(d\), endpoint-free compression, `NaN`)
was found by that search, not by accepting the recurrence.

---

## 11. Empirical evidence versus proof

| Claim | Status |
|---|---|
| (R) is exact for rooted binary \(T\) and symmetric finite \(d\) | **Proven** (Theorem 3.1) |
| (R3) is equivalent to (R) | **Proven** (Lemma 7.1) |
| \((v,u,w)\) is sufficient | **Proven** (Theorem 4.1) |
| neither endpoint is droppable | **Proven** (Theorem 4.2, Counterexample 5.3) |
| unconstrained-subtree optimal substructure fails | **Proven** (Theorem 5.2, Counterexample 5.3) |
| (R) is inexact for asymmetric \(d\) | **Proven** (Theorem 3.2, Counterexample 3.4) |
| D1–D5 are safe dominance rules | **Proven** (Section 6.1) |
| U1–U5 are unsafe | **Proven** by counterexample (Section 6.2) |
| \(\sum_v a_v b_v = \binom{n}{2}\) | **Proven** (LCA counting) |
| \(T_{\mathrm{(R)}} = \Theta(n^4)\) worst case; \(T_{\mathrm{(R3)}} = \Theta(n^3)\) | **Proven** (Section 8) |
| Brandes does not apply to arbitrary phylogenies | **Proven** (Brandes’ own incomplete-tree obstruction; comb padding is exponential) |
| (R) matched brute force on 258 random / structured symmetric instances, \(n \le 8\) | **Empirical only** — supports, does not replace, Theorem 3.1 |
| authors’ \(n = 4..10\) brute-force check | **Empirical only** |
| no \(o(n^3)\) exact algorithm for general binary trees and arbitrary \(d\) | **Conjectured** (Section 7.6); not a theorem |
| D4/D5 help on typical trait data | **Conjectured / historically reported** (Bar-Joseph 2001 tables); not re-measured here |
| `rotateConstr` realizes the DP order on every ape tree | **Not proven here**; outside the DP. A mismatch would be an ape bug or a contract bug, not a recurrence bug |
| Euclidean trait \(d\) has extra structure that a specialized OLO could exploit | **Conjectured**; no algorithm given |

No empirical row is promoted to a theorem in Section 12.

---

## 12. Final theoretical verdict

### Exactness

**Globally exact if and only if** \(T\) is rooted and strictly binary and
\(d\) is a symmetric finite real matrix. Then Theorem 3.1 applies: the
returned order is a true minimizer of \(J\) over the \(2^{n-1}\)
\(T\)-consistent layouts.

**Not globally exact** for asymmetric \(d\) (Theorem 3.2; smallest witness
\(n=3\)). In that regime the algorithm can report a scalar that is not the
cost of any feasible order.

**Not defined** for polytomies or unrooted trees; the implementation
correctly rejects the former and silently assumes the latter are rooted.

Non-negativity, the triangle inequality, and the trait-Euclidean origin of
\(d\) are **not** required for exactness.

### Required state

The minimal sufficient state that preserves global optimality, for
symmetric \(d\), is

\[
\bigl(v,\; u \in L(\mathrm{lc}(v)),\; w \in L(\mathrm{rc}(v)),\; F_v(u,w)\bigr),
\]

i.e. one cost per cross-child endpoint pair at each node, with reversal
identification. Both endpoints are necessary (Theorem 4.2). The internal
order is not. For asymmetric \(d\) the ordered pairs \((u,w)\) and
\((w,u)\) become distinct necessary states.

This is already (up to the transpose fill) the state the code uses. The
state is not over-large in the worst case over \(d\) (Theorem 4.3). It *is*
over-large on combs, where only \(O(s_v)\) pairs are feasible.

### Optimal substructure

Holds for the endpoint state (Theorem 5.1). Fails for the endpoint-free
state (Theorem 5.2). There is no greedy-choice property. Overlapping
subproblems hold, and the dependency DAG is the tree.

### Dominance

Safe to eliminate: duplicate reverse cells when \(d = d^\top\) (D2);
infeasible same-child pairs (D3); same-interface worse internals (D1);
pointwise-seam-dominated candidates (D4); bound-pruned candidates in the
inner min (D5).

Unsafe to eliminate: any cell solely because it has worse internal cost
than a cell with a *different* interface (U1, U2, Counterexample 5.3).

### Complexity

| Algorithm | Time | Space | Exact on general binary \(T\), symmetric \(d\)? |
|---|---|---|---|
| supplied (R) | \(\Theta(n^4)\) worst case | \(\Theta(n^2)\) | yes |
| factored (R3) | \(\Theta(n^3)\) worst case, tight | \(\Theta(n^2)\) | yes |
| (R3) + finite-entry loops, comb | \(\Theta(n^2)\) | \(\Theta(n^2)\) | yes (on combs) |
| Brandes 2007 | \(O(n^2\log n)\) | \(O(n)\) | **no** (complete trees only) |
| brute force | \(\Theta(n\,2^{n-1})\) | \(O(n)\) | yes |

The tightest *justified general* bound for an exact algorithm is
**\(\Theta(n^3)\) time, \(\Theta(n^2)\) space** via (R3). The supplied
code is a factor-\(\Theta(n)\) slower realization of the same exact
state machine.

### Potential improvements (theoretically justified, exactness-preserving)

1. **Replace (R) by (R3).** Same theorem, \(O(n^4)\to O(n^3)\). This is
   the single most important change; it is the difference between the
   2001 nested-loop presentation and the cba / practical implementations.
2. **Iterate only finite entries of \(F_L\).** Recovers \(O(n^2)\) on
   combs without a separate code path; never hurts.
3. **Optional D4/D5 pruning** for typical trait matrices; does not improve
   worst-case.
4. **Two-orientation tables** if the API is ever opened to asymmetric
   scores.
5. **Guard `scale()`** against zero-variance columns, so \(d\) remains
   finite and Theorem 3.1’s hypothesis holds.
6. **Brandes** only when the phylogeny (or a bounded-height clustering
   dendrogram) is complete up to \(O(1)\) in height.

Do **not**: drop an endpoint; keep one order per subtree; identify
reversals on directed \(d\); expect Knuth–Yao, SMAWK, Dijkstra, or min-cost
flow to beat (R3) on this problem.

### Remaining gaps

- No theorem that \(O(n^3)\) is optimal for arbitrary binary trees and
  arbitrary symmetric \(d\). Slightly-subcubic min-plus MM is a theoretical
  (impractical) crack; a combinatorial \(O(n^{3-\varepsilon})\) is open.
- `ape::rotateConstr` is not part of the proof. Exactness is claimed for
  `olo_phylo`’s `(order, cost)`, not for the plotted embedding.
- Numerical stability of the `which(table == min_cost)` tie-break under
  floating-point noise is not analysed; it does not affect the *value*
  of the min, only which optimum is selected.
- Euclidean-trait structure (low-dimensional embeddings, identity-minus-
  correlation matrices, etc.) is unused. Whether it yields a faster *exact*
  algorithm is open.

### Confidence

**Theorem-level:** exactness of (R) and (R3) on rooted binary trees with
symmetric finite \(d\); inexactness without symmetry; sufficiency and
necessity of the two-endpoint state; failure of endpoint-free optimal
substructure; safe/unsafe dominance classification D1–D5 / U1–U5;
\(\Theta(n^4)\) vs \(\Theta(n^3)\) arithmetic complexity.

**Empirical:** 258 / 258 agreements of (R) and (R3) with an independent
brute-force oracle on symmetric instances with \(n \le 8\); 43 / 258
failures of the endpoint-free algorithm on the same set; 43 / 45 failures
of (R) on directed instances, all repaired by two-orientation tables.

**Conjecture:** (R3) is asymptotically optimal, up to subcubic min-plus
improvements, for the general problem the code claims to solve.
