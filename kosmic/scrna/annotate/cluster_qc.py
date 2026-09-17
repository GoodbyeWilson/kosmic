# Cluster-level quality flags for the Annotate table.
#
# The QC step filters cells one at a time, and a doublet or a stressed
# nucleus passes every per-cell threshold: a doublet has more genes than
# a real cell and ordinary mitochondrial content; a damaged nucleus at
# 0.7% mitochondrial reads clears a 10% cut-off. They only become visible
# after clustering, as a cluster that matches no reference type well,
# whose markers are mitochondrial genes, whose doublet score is high, or
# which the authors of the deposit could not label. Across the five DCM
# studies these were 2-3% of cells per study and every one of them was
# labelled with confidence 'High' -- the margin rule there was 0.05.
#
# 'flag_clusters' reads those signals, all of which are already in the
# file after clustering and annotation, and returns a short reason string
# per cluster. Nothing is decided here: the Annotate table shows the
# flags so that the person deciding looks at those clusters first.
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from kosmic.scrna.annotate.author_labels import author_breakdown, author_label_series

# Thresholds, from the five DCM studies (Koenig, Reichart-LV, Guo, Youness,
# Chaffin, HeartMap LV reference): every clean cluster scored >= 0.85 with
# a margin >= 0.17; every junk cluster failed at least one rule below.
from kosmic.scrna.annotate.reference import AMBIGUOUS_MARGIN as WEAK_MARGIN, WEAK_SCORE  # noqa: E402
STRESS_MARKERS_MIN = 3          # of the top 10 markers
STRESS_PREFIXES = ('MT-', 'RPL', 'RPS', 'MRPL', 'MRPS')
DOUBLET_SCORE_RATIO = 2.0       # cluster median / study median
DOUBLET_GENES_RATIO = 1.3       # cluster median n_genes / median of its assigned type
AUTHOR_MIXED_FRACTION = 0.60
DOUBLET_COLUMNS = ('doublet_score', 'scrublet_score', 'scrublet_score_z',
                   'cellranger_doublet_scores', 'scrublet')


@dataclass
class ClusterFlags:
    cluster: str
    reasons: list[str] = field(default_factory=list)

    @property
    def flagged(self) -> bool:
        return bool(self.reasons)

    def text(self) -> str:
        return '; '.join(self.reasons)


def doublet_score_column(obs) -> Optional[str]:
    """The first numeric doublet-score column present in obs, or None."""
    for col in DOUBLET_COLUMNS:
        if col in obs.columns and pd.api.types.is_numeric_dtype(obs[col]):
            return col
    return None


def _stress_marker_count(names) -> int:
    return sum(1 for g in names if str(g).upper().startswith(STRESS_PREFIXES))


def flag_clusters(adata, cluster_col: str = 'leiden',
                  type_col: str = 'cell_type',
                  author_col: Optional[str] = 'cell_type_author',
                  details: Optional[dict] = None,
                  top_n_markers: int = 10) -> dict[str, ClusterFlags]:
    """Quality flags per cluster, from what clustering and annotation left in the file.

    Parameters
    ----------
    adata
        Clustered study; ``obs[cluster_col]`` required. Cells with no
        cluster (an excluded arm) are ignored.
    type_col
        The assigned labels; used for the gene-count comparison against
        the cluster's own type. Skipped when absent.
    author_col
        The authors' labels, if the deposit carried any. Skipped when absent.
    details
        ``uns['cluster_annotation_details']`` (best score, runner-up,
        margin per cluster). Skipped when None.

    Returns
    -------
    dict
        cluster label -> ClusterFlags; every cluster is present, flagged
        or not.
    """
    obs = adata.obs
    clusters = obs[cluster_col]
    has_cluster = clusters.notna().to_numpy()
    labels = sorted(set(clusters[has_cluster].astype(str)),
                    key=lambda x: (0, int(x)) if x.isdigit() else (1, x))
    out = {c: ClusterFlags(c) for c in labels}
    details = details or {}

    # Marker genes, when the cluster step ranked them.
    rank = adata.uns.get('rank_genes_groups')
    names = pd.DataFrame(rank['names']) if rank is not None and 'names' in rank else None

    dbl_col = doublet_score_column(obs)
    study_dbl_median = (float(np.nanmedian(obs.loc[has_cluster, dbl_col]))
                        if dbl_col else None)
    genes_col = 'n_genes_by_counts' if 'n_genes_by_counts' in obs.columns else None
    types = obs[type_col].astype(str).to_numpy() if type_col in obs.columns else None
    authors = author_label_series(obs[author_col]) if author_col and author_col in obs.columns else None

    cl_str = clusters.astype(str).to_numpy()
    for c in labels:
        mask = has_cluster & (cl_str == c)
        f = out[c]

        d = details.get(c, {})
        score = d.get('best_score')
        margin = d.get('margin')
        if score is not None and float(score) < WEAK_SCORE:
            f.reasons.append(f"weak match ({float(score):.2f})")
        elif margin is not None and float(margin) < WEAK_MARGIN:
            f.reasons.append(f"runner-up close ({d.get('runner_up')}, margin {float(margin):.2f})")

        if names is not None and c in names.columns:
            n_stress = _stress_marker_count(names[c].head(top_n_markers))
            if n_stress >= STRESS_MARKERS_MIN:
                f.reasons.append(f"{n_stress} of top {top_n_markers} markers are MT/ribosomal")

        if dbl_col and study_dbl_median is not None and study_dbl_median > 0:
            med = float(np.nanmedian(obs.loc[mask, dbl_col]))
            if med >= DOUBLET_SCORE_RATIO * study_dbl_median:
                f.reasons.append(f"doublet score {med / study_dbl_median:.1f}x study median ({dbl_col})")
        if authors is not None and mask.any():
            b = author_breakdown(authors[mask])
            n = int(mask.sum())
            if b.n_labelled == 0 or b.n_unlabelled > n / 2:
                f.reasons.append("authors left most cells unlabelled")
            elif b.top[1] < AUTHOR_MIXED_FRACTION * b.n_labelled:
                f.reasons.append(f"authors' labels mixed (top {100 * b.top[1] / b.n_labelled:.0f}%)")

        # Gene count is supporting evidence only: two clusters of the same
        # type can differ in depth for honest reasons (Chaffin's two
        # cardiomyocyte clusters differ 2x), so it never flags on its own.
        if f.reasons and genes_col and types is not None and mask.any():
            t = pd.Series(types[mask]).mode()
            t = t.iloc[0] if len(t) else None
            others = has_cluster & ~mask & (types == t)
            if t is not None and others.any():
                ref = float(obs.loc[others, genes_col].median())
                med = float(obs.loc[mask, genes_col].median())
                if ref > 0 and med >= DOUBLET_GENES_RATIO * ref:
                    f.reasons.append(f"{med / ref:.1f}x the genes of other {t}")

    return out
