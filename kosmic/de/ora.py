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
        'FDR'. Empty if no term had at least one overlap (or no terms
        passed 'min_genes').
    """
    from scipy.stats import fisher_exact

    study_set = set(study_genes)
    bg_set = set(background_genes)
    study_in_bg = study_set & bg_set
    N = len(bg_set)

    rows: list[dict] = []
    n_terms = len(library)
    for i, (term, genes) in enumerate(library.items()):
        if progress_callback is not None and i % 50 == 0:
            progress_callback(i, n_terms)

        term_genes = set(genes) & bg_set
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
            'Genes': sorted(overlap),
        })

    if progress_callback is not None:
        progress_callback(n_terms, n_terms)

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows).sort_values('P_value').reset_index(drop=True)
    df['FDR'] = bh_fdr(df['P_value'])
    return df
