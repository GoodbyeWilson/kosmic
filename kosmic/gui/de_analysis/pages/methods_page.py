# DE Analysis Methods Page: dataset stats, parameters, filters,
# results, and Vancouver-numbered references for the tools that
# actually ran on this data.
from __future__ import annotations

from PyQt6.QtWidgets import (
    QHBoxLayout, QLabel, QPushButton,
    QTextEdit,
)
from PyQt6.QtGui import QFont, QGuiApplication
import numpy as np
from kosmic import DEFAULT_FDR, DEFAULT_LFC_THRESHOLD
from kosmic.gui.shared.widgets import SecondaryLabel, SimplePage


# Reference catalogue. Each entry maps a "feature" key (something the
# pipeline either did or didn't do) to a Vancouver-style citation. The
# methods page only emits the ones that fired.
_REFS = {
    'scanpy': "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19(1):15.",
    'anndata': "Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Annotated data. bioRxiv. 2021. doi:10.1101/2021.12.16.473007.",
    'scrublet': "Wolock SL, Lopez R, Klein AM. Scrublet: computational identification of cell doublets in single-cell transcriptomic data. Cell Syst. 2019;8(4):281-291.e9.",
    'soupx': "Young MD, Behjati S. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. Gigascience. 2020;9(12):giaa151.",
    'harmony': "Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods. 2019;16(12):1289-1296.",
    'leiden': "Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. Sci Rep. 2019;9(1):5233.",
    'pseudobulk': "Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. Nat Commun. 2021;12(1):5692.",
    'deseq2': "Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550.",
    'apeglm': "Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. 2019;35(12):2084-2092.",
    'bh': "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc B. 1995;57(1):289-300.",
    'decontx': "Yang S, Corbett SE, Koga Y, et al. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. Genome Biol. 2020;21(1):57.",
    'celltypist': "Domínguez Conde C, Xu C, Jarvis LB, et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science. 2022;376(6594):eabl5197.",
    'singscore': "Foroutan M, Bhuva DD, Lyu R, et al. Single sample scoring of molecular phenotypes. BMC Bioinformatics. 2018;19(1):404.",
    'tirosh': "Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science. 2016;352(6282):189-196.",
}


class MethodsPage(SimplePage):
    """Values + references fact sheet."""

    help_id = "de/methods"
    CONTENT_MARGINS = (8, 8, 8, 8)
    CONTENT_SPACING = 6

    def __init__(self, workspace):
        super().__init__()
        self.ws = workspace
        self._setup_ui()

    def _setup_ui(self):
        layout = self.body_layout

        title = QLabel("Methods")
        title_f = title.font()
        title_f.setPointSize(title_f.pointSize() + 4)
        title_f.setBold(True)
        title.setFont(title_f)
        layout.addWidget(title)

        subtitle = SecondaryLabel(
            "Values and parameters that describe this analysis, plus "
            "references for the tools used. Updates live as filters change."
        )
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        row = QHBoxLayout()
        row.setSpacing(6)
        self._copy_btn = QPushButton("Copy to clipboard")
        self._copy_btn.clicked.connect(self._copy_text)
        row.addWidget(self._copy_btn)
        self._refresh_btn = QPushButton("Refresh")
        self._refresh_btn.clicked.connect(self._refresh)
        row.addWidget(self._refresh_btn)
        row.addStretch()
        layout.addLayout(row)

        self._text = QTextEdit()
        self._text.setReadOnly(True)
        font = QFont("Consolas, 'Courier New', monospace")
        font.setStyleHint(QFont.StyleHint.Monospace)
        self._text.setFont(font)
        layout.addWidget(self._text, 1)

    def on_activated(self):
        self._refresh()

    def _copy_text(self):
        QGuiApplication.clipboard().setText(self._text.toPlainText())

    # ------------------------------------------------------------------
    # Compose the sheet
    # ------------------------------------------------------------------

    def _refresh(self):
        adata = self.ws.current_adata
        if adata is None:
            self._text.setPlainText("No dataset loaded.")
            return
        sections = []
        used_refs = set()

        def _safe(header, fn):
            # One failing section must not blank the whole sheet.
            try:
                facts = fn()
            except Exception as e:  # noqa: BLE001
                sections.append((header, [f"(unavailable: {e})"]))
                return
            if facts:
                sections.append((header, facts))

        # References fire on evidence in the data, independent of whether the
        # preprocessing text came from provenance or inference -- so a recorded
        # run still cites Leiden, Scrublet, DecontX, etc.
        self._scan_evidence_refs(adata, used_refs)

        _safe("DATASET", self._dataset_facts)
        _safe("PREPROCESSING", lambda: self._preprocessing_facts(used_refs))
        _safe("GENE DIFFERENTIAL EXPRESSION", lambda: self._de_facts(used_refs))
        _safe("FILTERS", lambda: self._filter_facts(used_refs))
        _safe("RESULTS", self._results_facts)
        _safe("PATHWAY DIFFERENTIAL EXPRESSION",
              lambda: self._pathway_de_facts(used_refs))

        ref_lines = self._reference_lines(used_refs)
        if ref_lines:
            sections.append(("REFERENCES", ref_lines))

        sw_lines = self._software_lines()
        if sw_lines:
            sections.append(("SOFTWARE", sw_lines))

        self._text.setPlainText(self._render(sections))

    @staticmethod
    def _render(sections) -> str:
        """
        Render sections of (header, list-of-(key, value)-or-string) to text.

        'items' may be a list of '(key, value)' tuples (rendered as
        aligned 'key: value' rows) or plain strings (rendered as-is,
        for free-form lines like reference citations).
        """
        out = []
        for header, items in sections:
            if not items:
                continue
            out.append(header)
            # Align key column for tuple entries
            tuple_items = [it for it in items if isinstance(it, tuple)]
            key_width = max((len(k) for k, _ in tuple_items), default=0)
            for it in items:
                if isinstance(it, tuple):
                    k, v = it
                    out.append(f"  {k:<{key_width}}  {v}")
                else:
                    out.append(f"  {it}")
            out.append("")
        return "\n".join(out).rstrip() + "\n"

    # ------------------------------------------------------------------
    # Sections
    # ------------------------------------------------------------------

    def _dataset_facts(self) -> list:
        ws = self.ws
        adata = ws.current_adata
        items = []
        from pathlib import Path
        if ws.h5ad_path:
            items.append(("File:", Path(ws.h5ad_path).name))
        if ws.gse_accession:
            # Pull out just the GSE/GSM accession if the user typed the
            # full filename. Avoids "Accession: GSE183852_..._endothelial_cell"
            # noise in the methods sheet.
            import re as _re
            m = _re.match(r'(GSE\d+|GSM\d+|E-MTAB-\d+|PRJNA\d+)',
                          str(ws.gse_accession), _re.IGNORECASE)
            acc_clean = m.group(1).upper() if m else ws.gse_accession
            items.append(("Accession:", acc_clean))
        items.append(("Cells × genes:", f"{adata.n_obs:,} × {adata.n_vars:,}"))

        # Cell-type subset
        if 'cell_type' in adata.obs.columns:
            ct_unique = adata.obs['cell_type'].dropna().astype(str).unique()
            if len(ct_unique) == 1:
                items.append(("Cell type:", f"{ct_unique[0]} (subset)"))

        # Role split
        if '_role' in adata.obs.columns:
            roles = adata.obs['_role'].astype(str)
            n_d = int((roles == 'disease').sum())
            n_c = int((roles == 'control').sum())
            n_x = int((roles == 'exclude').sum())
            if ws.sample_col and ws.sample_col in adata.obs.columns:
                d_donors = int(adata.obs.loc[roles == 'disease',
                                              ws.sample_col].nunique())
                c_donors = int(adata.obs.loc[roles == 'control',
                                              ws.sample_col].nunique())
                items.append(
                    ("Disease group:", f"{d_donors} donors, {n_d:,} cells"))
                items.append(
                    ("Control group:", f"{c_donors} donors, {n_c:,} cells"))
            else:
                items.append(("Disease group:", f"{n_d:,} cells"))
                items.append(("Control group:", f"{n_c:,} cells"))
            if n_x > 0:
                items.append(("Excluded by role:", f"{n_x:,} cells"))

        # Role source
        role_col = adata.uns.get('role_condition_col')
        role_map = adata.uns.get('role_map')
        if role_col and isinstance(role_map, dict):
            d_keys = sorted(k for k, v in role_map.items() if v == 'disease')
            c_keys = sorted(k for k, v in role_map.items() if v == 'control')
            items.append(
                ("Role source column:", str(role_col)))
            items.append(
                ("Role values disease:", ", ".join(map(str, d_keys)) or "—"))
            items.append(
                ("Role values control:", ", ".join(map(str, c_keys)) or "—"))
        return items

    def _load_provenance(self):
        """The study's provenance record (next to the h5ad), or None."""
        h5ad = getattr(self.ws, 'h5ad_path', None)
        if not h5ad:
            return None
        from pathlib import Path
        from kosmic import provenance
        return provenance.load(Path(h5ad).parent)

    def _count_source_label(self):
        """Which counts fed the DE: prefer the recorded run, fall back to the
        live selection. Returns a display string or None."""
        def _label(cs):
            if cs in ('decontX_counts', 'decontx'):
                return "decontaminated (DecontX)"
            if cs in (None, 'raw'):
                return "raw counts"
            return str(cs)

        prov = self._load_provenance()
        if prov:
            from kosmic import provenance
            st = provenance.last_stage(prov, 'gene_de')
            if st is not None:
                cs = (st.get('params') or {}).get('count_source')
                if cs is not None:
                    return _label(cs)

        gp = getattr(self.ws, 'gene_de_page', None)
        if gp is not None and hasattr(gp, '_selected_counts_layer'):
            return _label(gp._selected_counts_layer())
        return None

    def _detection_pct_value(self):
        """Pre-DESeq2 detection filter as a percentage, provenance-first."""
        prov = self._load_provenance()
        if prov:
            from kosmic import provenance
            st = provenance.last_stage(prov, 'gene_de')
            if st is not None:
                v = (st.get('params') or {}).get('detection_min_pct')
                if v is not None:
                    v = float(v)
                    # Stored as a fraction (0-1); a value > 1 is a percentage.
                    return v * 100.0 if v <= 1.0 else v
        gp = getattr(self.ws, 'gene_de_page', None)
        if gp is not None and hasattr(gp, 'detection_pct_spin'):
            return float(gp.detection_pct_spin.value())  # already a percentage
        return None

    def _deseq2_filter_states(self):
        """(independent_on, cooks_on) for the DESeq2 run, provenance-first."""
        ind = cook = None
        prov = self._load_provenance()
        if prov:
            from kosmic import provenance
            st = provenance.last_stage(prov, 'gene_de')
            if st is not None:
                p = st.get('params') or {}
                ind = p.get('independent_filter')
                cook = p.get('cooks_filter')
        gp = getattr(self.ws, 'gene_de_page', None)
        if ind is None and gp is not None and hasattr(gp, 'indep_filter_check'):
            ind = gp.indep_filter_check.isChecked()
        if cook is None and gp is not None and hasattr(gp, 'cooks_filter_check'):
            cook = gp.cooks_filter_check.isChecked()
        return bool(ind), bool(cook)

    def _filterbyexpr_label(self):
        """Sample-level expression filter settings, provenance-first."""
        mc = ms = None
        prov = self._load_provenance()
        if prov:
            from kosmic import provenance
            st = provenance.last_stage(prov, 'gene_de')
            if st is not None:
                p = st.get('params') or {}
                mc = p.get('filter_min_count')
                ms = p.get('filter_min_samples')
        gp = getattr(self.ws, 'gene_de_page', None)
        if mc is None and gp is not None and hasattr(gp, 'filter_min_count_spin'):
            mc = gp.filter_min_count_spin.value()
        if ms is None and gp is not None and hasattr(gp, 'filter_min_samples_spin'):
            ms = int(gp.filter_min_samples_spin.value())
        if mc is None:
            return None
        samples = ("≥ smaller group size (auto)" if not ms or ms <= 0
                   else f"≥ {int(ms)} donors")
        return (f"CPM above the min-count {mc:g} cutoff in {samples}")

    def _pathway_de_facts(self, used_refs: set) -> list:
        """Pathway-scoring settings + coverage outcome, when a pathway DE run
        exists. Provenance-first, falling back to the live page state."""
        ws = self.ws
        if getattr(ws, 'pathway_de_results', None) is None:
            return []
        from kosmic.de.pathway_scoring import SCORING_METHODS

        method_key = cs = None
        gdp = cov = n_sets = None
        prov = self._load_provenance()
        st = None
        if prov:
            from kosmic import provenance
            st = provenance.last_stage(prov, 'pathway_de')
        if st is not None:
            p = st.get('params') or {}
            method_key = p.get('scoring_method')
            cs = p.get('count_source')
            gdp = p.get('gene_detection_pct')
            cov = p.get('min_coverage')
            n_sets = p.get('n_gene_sets')
        else:
            page = getattr(ws, 'pathway_de_page', None)
            if page is not None:
                method_key = page._get_scoring_method()
                cs = page._selected_counts_layer() or 'raw'
                gdp = page.gene_detection_spin.value() / 100.0
                cov = page.min_coverage_spin.value() / 100.0
                n_sets = len(getattr(ws, 'pathway_gene_sets', None) or {})

        items = []
        if method_key is not None:
            items.append(("Scoring method:",
                          SCORING_METHODS.get(method_key, str(method_key))))
            # Cite the published scorers only when they actually ran; the
            # DESeq2 per-donor methods need no extra citation beyond DESeq2.
            if method_key == 'singscore':
                used_refs.add('singscore')
            elif method_key == 'score_genes':
                used_refs.add('tirosh')
            elif method_key in ('deseq2', 'deseq2_vst', 'deseq2_vst_zscore'):
                used_refs.add('deseq2')
        if cs is not None:
            items.append(("Count source:",
                          "decontaminated (DecontX)"
                          if cs in ('decontX_counts', 'decontx') else "raw counts"))
            if cs in ('decontX_counts', 'decontx'):
                used_refs.add('decontx')
        if gdp is not None:
            items.append(("Gene detection filter:",
                          f"≥ {gdp * 100:g}% of cells"))
        if cov is not None:
            items.append(("Min pathway coverage:",
                          f"≥ {cov * 100:g}% of annotated genes"))
        if n_sets:
            items.append(("Pathways in gene set:", f"{n_sets:,}"))

        report = getattr(ws, 'pathway_coverage_report', None)
        if report is not None and not report.empty and 'status' in report.columns:
            scored = report['status'] == 'scored'
            dropped = report[~scored]
            items.append(("Coverage outcome:",
                          f"{int(scored.sum()):,} scored, {len(dropped):,} dropped"))
            # Name the excluded pathways (+ coverage) -- needed for the paper.
            if not dropped.empty and 'pathway' in dropped.columns:
                names = []
                for _, r in dropped.iterrows():
                    cov = r.get('coverage')
                    try:
                        cov_str = (f" ({float(cov) * 100:.0f}% coverage)"
                                   if np.isfinite(float(cov)) else "")
                    except (TypeError, ValueError):
                        cov_str = ""
                    names.append(f"{r['pathway']}{cov_str}")
                items.append(("Dropped pathways:", "; ".join(names)))
        return items

    @staticmethod
    def _scan_evidence_refs(adata, used_refs: set) -> None:
        """Add references for pipeline steps evidenced in the data itself, so
        citations don't depend on the preprocessing text path."""
        obs, var = adata.obs.columns, adata.var.columns
        uns, layers, obsm = adata.uns, adata.layers, adata.obsm
        used_refs.add('anndata')
        if 'log1p' in uns or 'highly_variable' in var or 'X_pca' in obsm:
            used_refs.add('scanpy')
        if 'leiden' in obs:
            used_refs.add('scanpy')
            used_refs.add('leiden')
        if 'predicted_doublet' in obs:
            used_refs.add('scrublet')
        if 'counts_pre_soupx' in layers:
            used_refs.add('soupx')
        if 'X_pca_harmony' in obsm:
            used_refs.add('harmony')
        if 'decontX_counts' in layers:
            used_refs.add('decontx')
        if 'cell_type_auto' in obs:
            used_refs.add('celltypist')

    @staticmethod
    def _provenance_preprocessing(prov, used_refs: set) -> list:
        """Preprocessing lines from the recorded scRNA stages (accurate)."""
        scrna = ('load', 'qc', 'normalize', 'cluster', 'decontx', 'subset',
                 'setup', 'roles')
        lines = []
        for st in prov.get('stages', []):
            if st.get('stage') not in scrna:
                continue
            if st['stage'] in ('normalize', 'cluster'):
                used_refs.add('scanpy')
            if st['stage'] == 'cluster':
                used_refs.add('leiden')
            if st['stage'] == 'decontx':
                used_refs.add('decontx')
            from kosmic.provenance import _STAGE_TITLES
            title = _STAGE_TITLES.get(st['stage'],
                                      st['stage'].replace('_', ' ').title())
            parts = []
            for k, v in (st.get('params') or {}).items():
                if isinstance(v, float):
                    v = f"{v:.4g}"
                elif isinstance(v, dict):
                    v = ", ".join(f"{a}={b}" for a, b in v.items()) or "(none)"
                elif isinstance(v, (list, tuple)):
                    v = ", ".join(str(x) for x in v) if v else "(none)"
                parts.append(f"{k}={v}")
            lines.append(f"{title}: {'; '.join(parts)}" if parts else title)
        return lines

    def _preprocessing_facts(self, used_refs: set) -> list:
        adata = self.ws.current_adata

        # Prefer the recorded provenance, following the 'source' pointer so a
        # cell-type subset shows its parent dataset's full processing first.
        h5ad = getattr(self.ws, 'h5ad_path', None)
        if h5ad:
            from pathlib import Path
            from kosmic import provenance
            lineage = provenance.load_lineage(Path(h5ad).parent)
            lines = []
            for i, (_dir, rec) in enumerate(lineage):
                if not rec:
                    continue
                if i < len(lineage) - 1:  # an upstream ancestor
                    src = rec.get('study') or '?'
                    lines.append(f"-- upstream dataset: {src} --")
                lines.extend(self._provenance_preprocessing(rec, used_refs))
            if lines:
                return lines

        obs = adata.obs.columns
        var = adata.var.columns
        items = []

        # QC
        if 'filtering_info' in adata.uns:
            info = adata.uns['filtering_info']
            if isinstance(info, dict):
                for key in ('min_genes', 'max_genes', 'min_counts',
                            'max_counts', 'max_mt', 'min_cells'):
                    if key in info and info[key] not in (None, 0, 0.0):
                        items.append(
                            (f"QC {key.replace('_', ' ')}:", str(info[key])))

        # Doublets
        if 'predicted_doublet' in obs:
            n_doub = int(adata.obs['predicted_doublet'].astype(bool).sum())
            pct = n_doub / max(adata.n_obs, 1) * 100
            items.append(
                ("Doublet detection:",
                 f"Scrublet ({n_doub:,} flagged, {pct:.1f}%)"))
            used_refs.add('scrublet')

        # SoupX
        if 'counts_pre_soupx' in adata.layers:
            items.append(
                ("Ambient correction:", "SoupX (counts replaced; "
                 "originals in layers['counts_pre_soupx'])"))
            used_refs.add('soupx')

        # Normalisation + HVG (scanpy). scanpy stores base=None in
        # uns['log1p'] to mean natural log (the sc.pp.log1p default);
        # spell that out rather than printing "None".
        if 'log1p' in adata.uns:
            base = adata.uns['log1p'].get('base')
            if base is None:
                base_str = "natural log"
            elif base == 2:
                base_str = "log2"
            elif base == 10:
                base_str = "log10"
            else:
                base_str = f"log base {base}"
            items.append(
                ("Normalisation:", f"library-size + {base_str}"))
            used_refs.add('scanpy')

        if 'highly_variable' in var:
            n_hvg = int(adata.var['highly_variable'].astype(bool).sum())
            items.append(
                ("HVG selection:", f"{n_hvg:,} genes (Seurat flavour)"))
            used_refs.add('scanpy')

        # PCA
        if 'X_pca' in adata.obsm:
            n_pcs = adata.obsm['X_pca'].shape[1]
            items.append(("PCA:", f"{n_pcs} components"))

        # Harmony
        if 'X_pca_harmony' in adata.obsm:
            items.append(("Batch correction:", "Harmony (harmonypy)"))
            used_refs.add('harmony')

        # Leiden
        if 'leiden' in obs:
            n_clusters = adata.obs['leiden'].nunique()
            res = None
            leiden_uns = adata.uns.get('leiden', {})
            if isinstance(leiden_uns, dict):
                params = leiden_uns.get('params', {})
                if isinstance(params, dict):
                    res = params.get('resolution')
            res_str = f", resolution {res}" if res is not None else ""
            items.append(
                ("Clustering:",
                 f"Leiden, {n_clusters} clusters{res_str} "
                 "(igraph backend, n_iterations=2)"))
            used_refs.add('leiden')

        # Cell-type annotation
        if 'cell_type_auto' in obs:
            items.append(
                ("Cell-type annotation:",
                 "CellTypist automated labelling (obs['cell_type_auto'])"))
            used_refs.add('celltypist')

        return items

    def _de_facts(self, used_refs: set) -> list:
        ws = self.ws
        items = []

        mode = getattr(ws, 'analysis_mode', None)
        if mode:
            items.append(("Analysis mode:",
                          "hypothesis-driven (pathway-restricted FDR)"
                          if mode == 'scoring'
                          else "discovery (genome-wide FDR)"))
        cs = self._count_source_label()
        if cs:
            items.append(("Count source:", cs))
            if 'DecontX' in cs:
                used_refs.add('decontx')

        if getattr(ws, 'unit', 'sample') == 'cell':
            items.append(("Approach:",
                          "cell-level (each cell treated as an observation; "
                          "no pseudobulk aggregation)"))
            items.append(("Method:",
                          "scanpy rank_genes_groups, Wilcoxon rank-sum, "
                          "disease vs control"))
            items.append(("Multiple testing:", "Benjamini-Hochberg FDR"))
            items.append(("Caveat:",
                          "cells within a donor are correlated, so cell-level "
                          "DE is pseudoreplicated and anticonservative relative "
                          "to pseudobulk (Squair et al. 2021)"))
            used_refs.add('scanpy')
            used_refs.add('pseudobulk')
            if ws.de_results is not None:
                items.append(("Genes tested:", f"{len(ws.de_results):,}"))
            return items

        items.append(("Approach:", "per-donor pseudobulk (sum of cell counts)"))
        used_refs.add('pseudobulk')

        if ws.de_method == 'deseq2':
            try:
                import pydeseq2
                ver = getattr(pydeseq2, '__version__', None)
            except Exception:
                ver = None
            method_str = "DESeq2"
            if ver:
                method_str += f" (pydeseq2 v{ver})"
            items.append(("Method:", method_str))
            items.append(("Design formula:", "~condition"))
            items.append(("Size factors:", "median-of-ratios"))
            items.append(("Dispersion:", "Cox-Reid maximum likelihood"))
            items.append(("Hypothesis test:", "Wald test"))
            ind_on, cook_on = self._deseq2_filter_states()
            items.append((
                "Cook's outlier filter:",
                "on (per-gene flagging of extreme-outlier samples)"
                if cook_on else "off"))
            items.append((
                "Independent filter:",
                "on (mean-expression-based independent filtering before "
                "adjusted p-value calculation)" if ind_on else "off"))
            items.append(("LFC shrinkage:",
                          "apeGLM heavy-tailed Cauchy prior"))
            items.append(("Effect estimates retained:",
                          "shrunk LFC for display; unshrunk LFC + SE "
                          "for downstream meta-analysis"))
            used_refs.add('deseq2')
            used_refs.add('apeglm')
        else:
            items.append(("Method:",
                          "Welch's t-test on log2(CPM+1) pseudobulk"))

        items.append(("Min cells per donor:", str(ws.min_cells)))
        if ws.min_counts:
            items.append(("Min transcripts per donor:", f"{ws.min_counts:,}"))
        det = self._detection_pct_value()
        if det is not None:
            items.append((
                "Detection filter:",
                f"expressed in ≥ {det:g}% of cells in a condition" if det > 0
                else "off (all genes tested)"))
        fbe = self._filterbyexpr_label()
        if fbe:
            items.append(("Sample-level filter:", fbe))
        if ws.de_results is not None:
            items.append(("Genes tested:", f"{len(ws.de_results):,}"))
        return items

    def _filter_facts(self, used_refs: set) -> list:
        gene_de_page = getattr(self.ws, 'gene_de_page', None)
        if gene_de_page is None:
            return []
        items = []

        # 'Pathway genes only' is a pure display filter: it hides non-pathway
        # dots but never touches FDR or significance.
        pw_check = getattr(gene_de_page, '_pathway_only_check', None)
        pathway_active = (
            pw_check is not None and pw_check.isChecked()
            and bool(getattr(self.ws, 'pathway_gene_sets', None))
        )
        if pathway_active:
            n_pw = len(self.ws.pathway_gene_sets)
            n_pw_genes = len({g for gs in self.ws.pathway_gene_sets.values()
                              for g in gs})
            items.append(
                ("Restrict to pathways (display):",
                 f"on ({n_pw} pathways, {n_pw_genes:,} unique genes)"))

        # Significance gate. FDR is computed once by the DE engine over the
        # mode-appropriate scope (genome-wide in discovery, pathway-restricted
        # in hypothesis); the gate reads that adjusted p, it does not re-run BH.
        use_raw = gene_de_page._raw_pval_check.isChecked()
        try:
            pval_thresh = float(gene_de_page.pval_filter.value())
            fc_thresh = float(gene_de_page.fc_filter.value())
        except Exception:
            pval_thresh = DEFAULT_FDR
            fc_thresh = DEFAULT_LFC_THRESHOLD
        if use_raw:
            items.append(
                ("Significance gate:", f"raw p < {pval_thresh:g} (no MTC)"))
        else:
            items.append(
                ("Significance gate:", f"BH-adjusted p (FDR) < {pval_thresh:g}"))
            # State the multiple-testing universe explicitly, so a
            # pathway-restricted FDR is not mistaken for a post-hoc genome-wide
            # correction.
            mode = getattr(self.ws, 'analysis_mode', None)
            pw = getattr(self.ws, 'pathway_gene_sets', None)
            if mode == 'scoring' and pw:
                n_pw_genes = len({g for gs in pw.values() for g in gs})
                universe = ("within the predefined pathway gene set "
                            f"({n_pw_genes:,} annotated genes)")
            else:
                universe = "across all tested genes (genome-wide)"
            items.append(("FDR universe:", f"Benjamini-Hochberg {universe}"))
            used_refs.add('bh')
        items.append(("|log2FC| threshold:", f"> {fc_thresh:g}"))

        return items

    def _results_facts(self) -> list:
        """Report both scopes: genome-wide (authoritative stored FDR) and, when
        a pathway set is loaded, the pathway-restricted result (BH within the
        set -- exactly what hypothesis mode computes). Independent of which mode
        the run used, so the sheet reads the same either way."""
        ws = self.ws
        gene_de_page = getattr(ws, 'gene_de_page', None)
        if gene_de_page is None or ws.de_results is None:
            return []

        use_raw = gene_de_page._raw_pval_check.isChecked()
        try:
            pval_thresh = float(gene_de_page.pval_filter.value())
            fc_thresh = float(gene_de_page.fc_filter.value())
        except Exception:
            pval_thresh, fc_thresh = DEFAULT_FDR, DEFAULT_LFC_THRESHOLD

        genome = self._genome_frame(ws.de_results)
        items = [("Genes tested:", f"{len(genome):,}")]

        gw = self._scope_stats(genome, use_raw, pval_thresh, fc_thresh,
                               recompute_bh=False)
        if gw:
            if gw['scope'] != len(genome):
                items.append(("Genome-wide FDR scope:", f"{gw['scope']:,}"))
            items.append(
                ("Genome-wide significant:",
                 f"{gw['sig']:,} ({gw['up']:,} up / {gw['down']:,} down)"))
            if gw['top']:
                items.append(("  top |log2FC|:", ", ".join(gw['top'])))

        pathway_genes = self._pathway_gene_set()
        if pathway_genes and 'names' in genome.columns:
            subset = genome[genome['names'].isin(pathway_genes)]
            pw = self._scope_stats(subset, use_raw, pval_thresh, fc_thresh,
                                   recompute_bh=True)
            if pw:
                items.append(("Pathway set tested:",
                              f"{pw['scope']:,} genes (BH within set)"))
                items.append(
                    ("Pathway-set significant:",
                     f"{pw['sig']:,} ({pw['up']:,} up / {pw['down']:,} down)"))
                if pw['top']:
                    items.append(("  top |log2FC|:", ", ".join(pw['top'])))
        return items

    # ------------------------------------------------------------------
    # References + software
    # ------------------------------------------------------------------

    def _reference_lines(self, used_refs: set) -> list:
        if not used_refs:
            return []
        # Stable order from the catalogue
        out = []
        n = 0
        for key in ('pseudobulk', 'deseq2', 'apeglm', 'bh',
                    'singscore', 'tirosh',
                    'scanpy', 'anndata',
                    'scrublet', 'soupx', 'decontx', 'harmony', 'leiden',
                    'celltypist'):
            if key in used_refs and key in _REFS:
                n += 1
                out.append(f"{n}. {_REFS[key]}")
        return out

    def _software_lines(self) -> list:
        bits = []

        def _ver(modname: str):
            try:
                import importlib
                m = importlib.import_module(modname)
                v = getattr(m, '__version__', None)
                if v is not None:
                    return f"{modname} {v}"
            except Exception:
                pass
            return None

        for mod in ('scanpy', 'anndata', 'pydeseq2', 'numpy', 'scipy',
                    'pandas', 'soupx', 'scrublet', 'harmonypy',
                    'celltypist', 'decontx'):
            line = _ver(mod)
            if line:
                bits.append(line)
        import sys
        bits.append(f"Python {sys.version.split()[0]}")
        return ["; ".join(bits)] if bits else []

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    def _genome_frame(self, de_results):
        """The genome-wide result frame regardless of run mode: hypothesis mode
        stashes it under attrs['genome_wide_ranking'] (its own 'de_results' is
        pathway-restricted); discovery mode's 'de_results' is already genome-wide.
        """
        gw = de_results.attrs.get('genome_wide_ranking')
        return gw if gw is not None else de_results

    def _pathway_gene_set(self) -> set:
        """Union of tested genes across the loaded pathway coverage, or empty."""
        cov = getattr(self.ws, 'pathway_coverage', None)
        genes = set()
        if cov:
            for info in cov.values():
                genes.update(info.get('genes', []))
        return genes

    def _scope_stats(self, df, use_raw, pval_thresh, fc_thresh, *, recompute_bh):
        """Significance tally over 'df'. Adjusted p is BH over df's own raw p
        when 'recompute_bh' (pathway-restricted scope, matching hypothesis
        mode); otherwise the stored 'pvals_adj' is used (authoritative
        genome-wide FDR). Returns a dict or None."""
        if df is None or df.empty or 'logfoldchanges' not in df.columns:
            return None
        lfc = df['logfoldchanges'].astype(float).to_numpy()
        if use_raw:
            if 'pvals' not in df.columns:
                return None
            padj = df['pvals'].astype(float).to_numpy()
        elif recompute_bh:
            if 'pvals' not in df.columns:
                return None
            from kosmic.numerical import bh_fdr
            padj = bh_fdr(df['pvals'].astype(float).to_numpy())
        else:
            if 'pvals_adj' not in df.columns:
                return None
            padj = df['pvals_adj'].astype(float).to_numpy()

        sig = (padj < pval_thresh) & (np.abs(lfc) > fc_thresh)
        top = []
        if sig.any() and 'names' in df.columns:
            names = df['names'].astype(str).to_numpy()
            for idx in np.argsort(-np.abs(lfc)):
                if sig[idx]:
                    top.append(names[idx])
                if len(top) >= 10:
                    break
        return {
            'scope': int(np.isfinite(padj).sum()),
            'sig': int(sig.sum()),
            'up': int((sig & (lfc > 0)).sum()),
            'down': int((sig & (lfc < 0)).sum()),
            'top': top,
        }
