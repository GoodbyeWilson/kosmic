# ADR-012: CellTypist broad labels group lymphoid and myeloid subtypes, and the vote is over broad labels

- Status: Accepted by Kirk (2026-10-14).

## Problem

With *Broad labels*, CellTypist's fine labels are collapsed to broad
types by their prefix in the model's names (`vCM1`…`vCM5` →
Ventricular Cardiomyocyte). Majority voting chose each cluster's label
among the fine labels and collapsed the winner afterwards. Two things
followed on the DCM atlas, annotated with the `Healthy_Adult_Heart`
model:

- A type that the model splits into several fine labels lost the vote
  to a type with one, even where its cells were the majority. A T/NK
  cluster was named CD4+ T Cell from `CD4+T_naive` at 20.6% of its
  cells, while CD8 T cells were 40% and CD4 T cells 32%. After
  subclustering, a subcluster was named CD16+ Monocyte at 38% while 46%
  of its cells were LYVE1+ macrophages (three fine labels), and another
  CD4+ T Cell at 39% while 55% were CD8 T cells.
- The broad types kept distinctions that single-nucleus data does not
  support reliably: CD4 against CD8 T cells, and the macrophage subsets
  (`LYVE1+` macrophages against monocyte-derived macrophages). Each of
  these clusters needed relabelling by hand.

## Decision

1. **Broad labels group lymphoid and myeloid subtypes one level up.**
   CD4+ T, CD8+ T, MAIT, γδ T and cycling T/NK → T Cell; NK and ILC →
   NK Cell; LYVE1+ macrophages and MoMP → Macrophage; CD14+ and CD16+
   monocytes → Monocyte. Dendritic, mast, B, neutrophil and every
   non-immune type are unchanged. The mapping, like the existing one,
   is written for the heart model's label names.
2. **With broad labels, the majority vote is over broad types.**
   `majority_vote(clusters, labels, group_of)` takes the most common
   broad type in each cluster; the cluster's fine label
   (`cell_type_fine`) is the most common fine label within that type,
   so the fine and broad labels agree. Without broad labels the vote is
   over fine labels as before.
3. One vote serves both paths. A large dataset is predicted in chunks
   and voted by KOSMIC already; a small one was voted by CellTypist's
   own `majority_voting`, which is the same plurality rule. Both are now
   predicted cell by cell and voted by `majority_vote`.
4. `uns['celltypist_model']` records `majority_voting_level` ('broad'
   or 'fine'), so the annotation's provenance and methods state it.

## Consequences

- CD4/CD8 T cells and the macrophage subsets are no longer separate
  broad labels. Where a study supports them (whole-cell data, for
  instance), they are available as fine labels in `cell_type_fine` and
  can be given to a cluster or subcluster by hand.
- Files annotated before this keep their labels until re-annotated.
- On the DCM atlas the vote over broad types changed only the T/NK
  cluster (1 of 32 clusters); every cluster whose cells agreed kept its
  label.
