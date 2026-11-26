# Genes called (Venn)

An area-proportional Venn diagram of the genes that two to four methods of a
[Methods Comparison](../comparison/index.md) call.

- Each shape's area is proportional to the number of genes its method
  calls, so the diagram shows how many genes each method calls relative to
  the others. A method that calls few genes is drawn small.
- Each overlap's area is proportional to the genes the methods share. A
  method whose calls are all called by another method is drawn inside it.
  In the DCM ventricular cardiomyocytes, for example, vote counting lies
  inside REML, which lies almost entirely inside merged pseudobulk.
- The number in each region is the number of genes called by exactly those
  methods.

Two or three methods are drawn as circles, and four as ellipses. The
layout is the one the Meta-Analysis workspace's Venn diagram uses: the
shapes are fitted so that each region's area approaches its share of the
genes. With four methods the fit is approximate, so read the overlaps from
the numbers. A region too thin to hold its number is drawn without it; the
CSV written on export lists every gene.

The key beside the diagram names each method and the number of genes it
calls. The note under the diagram gives:

- the call rule chosen in the Methods Comparison workspace (by default an
  adjusted p-value below 0.05, and for vote counting its vote pass);
- the number of genes compared (the genes every method of the comparison
  tested), and how many of them none of the methods shown calls.

| Control | What it does |
|---|---|
| **First method**, **Second method** | The two methods always drawn. They open on merged pseudobulk and REML. |
| **Third method**, **Fourth method** | Further methods, or **None**. They open on Stouffer and vote counting, when the comparison includes them. |

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
