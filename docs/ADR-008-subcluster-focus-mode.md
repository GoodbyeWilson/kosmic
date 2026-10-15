# ADR-008: Subclustering a population in a focus mode

- Status: Proposed (2026-10-06). Agreed by Kirk; awaiting Calum.

## Problem

The Cluster step embeds, clusters and annotates a whole study at one
resolution, and one resolution does not suit every population. On an
atlas the major cell types separate well, but the states within them
(cardiomyocyte, fibroblast and endothelial subtypes, for example) are
either merged into one cluster or split in ways that the marker genes do
not explain. The highly variable genes and principal components are
chosen over every cell, so they are dominated by the genes that separate
the major types. Raising the resolution over the whole study splits the
largest populations first, and not necessarily along the lines that
separate their subtypes.

The only route KOSMIC offers today is the Subset step, which extracts a
cell type into a new study ("Subsets are studies" in CODEBASE.md). That
suits analysing a cell type on its own, but not labelling: the refined
labels stay in the derived study and never return to the parent. For
the atlas this means they cannot reach the per-study files through label
propagation.

## Decision

**In the Cluster step, the user can select a population on the embedding
and open it in a focus mode. Focus mode re-embeds, re-clusters and
annotates only that population. On leaving focus mode, the user either
applies the result, which writes subcluster ids and cell-type labels
back onto the full study, or discards it, which leaves the study
unchanged.**

### Selecting a population

1. With the embedding coloured by a categorical obs column (`leiden`,
   `cell_type` or another), clicking a cell selects its group.
   Ctrl-click adds further groups from the same column, so that, for
   example, two clusters that both look like fibroblasts can be
   subclustered together. Unselected groups are greyed out.
2. A button names the selection and its cell count, and opens focus
   mode. Clicking the embedding only selects; it never starts a
   computation.

### Inside focus mode

3. The population's included cells form a working subset
   (`workset.working_subset`). Cells in the population with
   `_role == 'exclude'` are left out, as in the PCA and cluster steps,
   and keep their labels. A population of fewer than 50 cells is not
   opened.
4. The subset gets its own embedding and clustering. Highly variable
   genes, PCA, Harmony (when used), the neighbour graph, UMAP and Leiden
   are all recomputed for it, and marker genes are ranked between the
   subclusters. The steps are the functions the PCA and cluster steps
   already use. Their settings default to the values in the study's
   last `embedding` and `cluster` provenance stages, including the
   Harmony batch key, and the user can change them before running.
   Harmony is skipped when the population's cells come from a single
   batch. The parent's neighbour graph is not reused
   (`sc.tl.leiden(restrict_to=...)`), because the parent's genes and
   components were chosen to separate the major types. Re-running with
   only a new resolution re-runs Leiden on the population's own graph and
   keeps its embedding and labels; any other change starts the
   population afresh.
5. The Annotate controls work on the subset: the cluster table, the
   top-genes panel, manual naming, and the automated methods (CellTypist,
   Reference Atlas, CellMarker 2.0 and PanglaoDB). The automated methods
   stay available until it is established which of them suit subtypes.
6. Nothing in focus mode writes the study's file. Only one population is
   in focus at a time, and a subcluster cannot be opened in focus mode
   while its parent's focus is still open.

### Applying the result

7. Applying writes the following onto the full study and saves its
   processed file:
   - `obs['leiden_sub']`. Cells outside the focus keep their current
     `leiden_sub`, or their `leiden` id when the column is new. Focused
     cells get `<groups>.<subcluster>`, where `<groups>` is the selected
     group, or the selected groups joined with `+` in the column's
     order: `7.0`, `7+12.0`, `Fibroblast.0`.
   - `obs['cell_type']`. Each named subcluster's cells are labelled
     through `set_cell_type`. Cells in a subcluster left unnamed keep
     their current label.
   - A `subcluster` provenance stage with:
     - the column the population was selected from;
     - each selected group and its number of cells;
     - the number of excluded cells left out;
     - the settings the subset was embedded and clustered with (number
       of HVGs and selection method, number of PCs, Harmony and its batch
       key, number of neighbours, resolution);
     - the annotation method with its settings, recorded as the
       `annotate` stage records them, or "manual";
     - each hand relabel, with the resolution it was made at, since a
       relabel made before a re-run at another resolution names a
       subcluster of that earlier clustering;
     - for each subcluster, its id, label and number of cells, and, when
       several groups were selected, how many of its cells came from each.

     Each applied focus appends a stage, so the record lists them in the
     order they were applied. The Methods page renders the stage under
     the title "Subclustering", and the methods document lists every
     `subcluster` stage, where for other stages it reports the last run.
8. The subset's embedding, graph and marker-gene table are discarded.
   The parent embedding is unchanged, and is shown recoloured by the new
   labels.
9. Once applied, a subcluster can be selected on the parent embedding
   (coloured by `leiden_sub`) and opened in focus mode in turn; its
   subclusters get ids such as `7.1.0`. A later focus on cells that
   already have subcluster ids overwrites them.
10. Leaving focus mode, switching study or closing the application with
    an unapplied result asks whether to apply or discard it.

### Later steps on the full study

11. Re-running clustering on the full study removes `leiden_sub`,
    because its ids no longer correspond to the new `leiden`. It keeps
    `cell_type`, which is a result in its own right.
12. Running automated annotation on the full study replaces `cell_type`
    for every cell, including labels applied from focus mode, and so does
    the Annotate step's *Undo Changes*. Both ask first when the study has
    `subcluster` stages.

### Propagation to studies

13. Label propagation already carries `cell_type`, so refined atlas
    labels reach each study as `cell_type_atlas`. It also carries
    `leiden_sub` when the atlas has that column. `leiden_sub` is written
    only as `leiden_sub_atlas`, never into a study's plain `leiden_sub`,
    because a study's own `leiden_sub` must correspond to its own
    `leiden`. Combine carries `leiden_sub_atlas` onto a rebuilt atlas, as
    it does `cell_type_atlas` and `leiden_atlas`.
14. The propagation stage recorded in each study names the atlas's
    `subcluster` stages alongside its `annotate` stage, so a study's
    methods describe how the labels it received were made.

### Code placement

15. The analysis is a Qt-free module, `kosmic/scrna/cluster/subcluster.py`,
    which builds the subset, runs the recomputation and applies the
    result. It is built on `workset.py` and the existing HVG, PCA,
    Harmony, neighbour, Leiden and marker-gene functions, and has tests.
16. Focus mode is a state of the existing Cluster page, not a new page,
    archetype or workspace.
17. The Annotate panel is refactored so that the caller sets the dataset
    it works on and where it may save. Today it always reads the
    workspace's `current_adata`, and every automated annotation worker
    writes the processed file when it finishes. In focus mode the panel
    reads the subset, its workers run without an output path, and its
    Save button is hidden; Apply, in the bar above the embedding, takes
    its place. This refactor is made and tested
    on its own, with the panel's behaviour outside focus mode unchanged,
    before focus mode is built on it.

## Consequences

- An atlas can be labelled in depth within the atlas, and the labels
  reach every study through the existing propagation. Extracting cell
  types into derived studies is no longer needed for labelling; the
  Subset step remains for analysing a cell type as its own study.
- `leiden` and the parent embedding are never changed by focus mode, so
  the major-type clustering stays as it was computed and recorded.
- A study's processed file gains `obs['leiden_sub']`, and its provenance
  gains `subcluster` stages. Code that does not know the column ignores
  it.
- A focus costs the time and memory of the PCA and cluster steps on the
  population's cells, plus a copy of their rows of the expression
  matrix. The largest population of an atlas can be a substantial
  fraction of the atlas. If sketch mode (ADR-005) is adopted, focus mode
  runs on the population's sketched cells and labels the rest by the
  same projection.
- The Annotate panel refactor is the main risk to existing behaviour,
  which is why it is a separate, tested change.
- Propagation changes in two places: `leiden_sub` is requested only when
  the atlas has it (requesting a missing column is an error), and it is
  exempt from filling the plain column.
- CODEBASE.md gains focus mode under the data contracts and the
  `subcluster` stage under Provenance. The Cluster help page
  (`scrna/cluster.md`), the atlas page (`project/atlas.md`) and the
  Cluster screenshots are updated.
