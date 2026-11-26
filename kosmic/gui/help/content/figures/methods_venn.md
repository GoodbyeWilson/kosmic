# Genes called (Venn)

A Venn diagram of the genes that two or three methods of a [Methods
Comparison](../comparison/index.md) call. Each circle is one method,
labelled with its name and the number of genes it calls. The number in each
region is the number of genes called by exactly those methods: in a diagram
of merged pseudobulk, REML and vote counting, the centre holds the genes all
three call, and the region of merged pseudobulk alone holds the genes only
merged pseudobulk calls.

The circles are all the same size and are not drawn to scale, so the
overlap is read from the numbers, not from the areas. The note under the
diagram gives the number of genes compared (the genes every method of the
comparison tested) and how many of them none of the methods shown calls.

A gene is called under the rule chosen in the Methods Comparison workspace,
which the title states: by default an adjusted p-value below 0.05, and for
vote counting its vote pass.

| Control | What it does |
|---|---|
| **First method**, **Second method** | The two methods always drawn. They open on merged pseudobulk and REML. |
| **Third method** | A third method, or **None** for a diagram of two. It opens on vote counting. |

The diagram draws the comparison held in the Methods Comparison workspace,
or else the newest comparison saved in the project. When the project has
more than one saved comparison, for example one per cell type, the **Run**
list above the figure list chooses which one is drawn.

Exporting saves the figure in the comparison's folder,
`methods_comparison/<counts>/<cell type>/`, and writes two files beside it:

- a CSV with one row per gene called by at least one of the methods shown,
  giving the methods that call it (`region`) and one column per method;
- `<figure name>_params.json`, with the comparison's parameters: the call
  rule, the merged pseudobulk design, where each method's result came from,
  and the methods shown.
