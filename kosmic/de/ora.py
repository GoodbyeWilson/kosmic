# Over-representation analysis (ORA).
#
# For each term in a gene-set library, tests whether the overlap between
# the term's genes and the user-supplied study set is larger than expected
# under the hypergeometric distribution (one-sided Fisher's exact test).
#
# Library-agnostic: pass any '{term -> list_of_genes}' mapping and you get
# back a DataFrame with the standard columns. Used by the Enrichr ORA path
# in the DE workspace; reusable for any other gene-set source (custom
# GMTs, KEGG, etc.).
#
# Gene symbols are matched case-insensitively (both sides upper-cased),
# the same convention as kosmic.de.fgsea.run_fgsea. Enrichr and MSigDB
# libraries use human upper-case symbols, so exact matching tests no terms
# on mouse data ('Col4a1' vs 'COL4A1'). Upper-casing is used rather than
# the bundled ortholog table (kosmic.reference.orthologs): it needs no
# species detection, keeps ORA and GSEA on one rule, and the table lists
# several mouse genes for some human genes (e.g. ALDH1A1 -> Aldh1a1 and
# Aldh1a7). Mouse genes whose human ortholog has a different symbol are
# not matched.
from __future__ import annotations

from typing import Callable, Iterable, Mapping, Optional

import pandas as pd

from kosmic.numerical import bh_fdr


def run_ora(study_genes: Iterable[str],
            background_genes: Iterable[str],
            library: Mapping[str, Iterable[str]],
            *,
            min_genes: int = 5,
            progress_callback: Optional[Callable[[int, int], None]] = None,
            ) -> pd.DataFrame:
    """One-sided Fisher's exact test for over-representation.

    Parameters
    ----------
    study_genes
        Significant / "selected" genes (e.g. DE hits at FDR < 0.05).
    background_genes
        Universe to test against -- typically all genes that were tested
        in the upstream DE step. Restricts both 'study_genes' and the
        per-term gene sets so the contingency table denominators are
        consistent.
    library
        Mapping 'term -> list_of_genes'. Genes outside the background
        are dropped per term before testing.

        Genes in all three inputs are matched case-insensitively. If two
        background genes differ only in case they count as one gene.
    min_genes
        Minimum number of background-restricted genes a term must contain
        to be tested. Defaults to 5.
    progress_callback
        Optional 'callback(i_done, n_total)' called every 50 terms.

    Returns
    -------
    pd.DataFrame
        Rows sorted ascending by p-value; columns: 'Term', 'P_value',
        'Fold_Enrichment', 'Gene_Count', 'Term_Size', 'Genes',
        'FDR'. 'Genes' lists the overlapping genes as they are spelled in
        'background_genes', not as in the library. Empty if no term had
        at least one overlap (or no terms passed 'min_genes').
    """
    from scipy.stats import fisher_exact

    # Upper-cased symbol -> the background's own spelling, for 'Genes'.
    # Sorted first so a case collision resolves the same way every run.
    bg_names: dict[str, str] = {}
    for g in sorted({str(g) for g in background_genes}):
        bg_names.setdefault(g.upper(), g)
    bg_set = set(bg_names)
    study_in_bg = {str(g).upper() for g in study_genes} & bg_set
    N = len(bg_set)

    rows: list[dict] = []
    n_terms = len(library)
    for i, (term, genes) in enumerate(library.items()):
        if progress_callback is not None and i % 50 == 0:
            progress_callback(i, n_terms)

        term_genes = {str(g).upper() for g in genes} & bg_set
        if len(term_genes) < min_genes:
            continue

        overlap = study_in_bg & term_genes
        a = len(overlap)
        if a == 0:
            continue
        b = len(study_in_bg) - a
        c = len(term_genes) - a
        d = N - a - b - c

        _, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        fold = ((a / max(len(study_in_bg), 1))
                / (len(term_genes) / max(N, 1)))

        rows.append({
            'Term': term,
            'P_value': pval,
            'Fold_Enrichment': fold,
            'Gene_Count': a,
            'Term_Size': len(term_genes),
            'Genes': sorted(bg_names[g] for g in overlap),
        })

    if progress_callback is not None:
        progress_callback(n_terms, n_terms)

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows).sort_values('P_value').reset_index(drop=True)
    df['FDR'] = bh_fdr(df['P_value'])
    return df
