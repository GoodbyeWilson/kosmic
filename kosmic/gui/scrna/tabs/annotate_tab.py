# Annotate Tab
# Cell type annotation: PanglaoDB, CellMarker 2.0, CellTypist.
# Cluster summary table with manual re-annotation.

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QTableWidget, QTableWidgetItem, QHeaderView,
    QCheckBox, QSpinBox,
    QListWidget, QListWidgetItem, QStackedWidget, QSplitter, QMenu, QSizePolicy,
)
from PyQt6.QtCore import QSettings, Qt, pyqtSignal
from pathlib import Path
from typing import Optional
import numpy as np
import pandas as pd

from kosmic.gui.shared.theme import NoScrollComboBox, NoScrollDoubleSpinBox
from kosmic.gui.shared.widgets import BaseWorker, Column, HintLabel, ResultsTable, SettingsGroup, StatusLabel
from kosmic.paths import processed_h5ad_path
from kosmic.reference.cell_type_focus import load_focus_presets
from kosmic import DEFAULT_FDR
from kosmic.gui.shared import dialogs, run_worker
from kosmic.scrna.annotate.author_labels import author_breakdown, author_label_series
from kosmic.scrna.annotate.cluster_qc import flag_clusters


_FOCUS_PRESET_QSETTING_KEY = 'annotate/last_focus_preset'


# ---------------------------------------------------------------------------
# Workers
# ---------------------------------------------------------------------------

class ReferenceAnnotationWorker(BaseWorker):
    """
    Worker for reference-based cell type annotation.

    Emits 'finished_ok' with a tuple '(adata, message)'.
    """

    def __init__(self, adata, reference_name_or_path: str, output_path: str,
                 min_correlation: float = 0.3):
        super().__init__()
        self.adata = adata
        self.reference_name_or_path = reference_name_or_path
        self.output_path = output_path
        self.min_correlation = min_correlation

    def _run(self):
        from pathlib import Path
        from kosmic.scrna.annotate.reference import run_reference_annotation

        # Preserve author annotations before overwriting
        if 'cell_type' in self.adata.obs.columns and 'cell_type_author' not in self.adata.obs.columns:
            self.adata.obs['cell_type_author'] = self.adata.obs['cell_type'].copy()
            self.progress.emit("Preserved author annotations in 'cell_type_author' column")

        self.progress.emit(f"Loading reference: {self.reference_name_or_path}...")
        self.progress_pct.emit(5)

        def _progress(msg):
            self.progress.emit(msg)

        self.adata, cluster_types, annotation_details = run_reference_annotation(
            self.adata,
            self.reference_name_or_path,
            min_correlation=self.min_correlation,
            progress_callback=_progress,
        )

        # Store details in uns for the cluster table
        self.adata.uns['cluster_annotation_details'] = annotation_details
        self.adata.obs['cell_type_auto'] = self.adata.obs['cell_type'].copy()

        self.progress.emit("Saving...")
        self.progress_pct.emit(95)
        Path(self.output_path).parent.mkdir(parents=True, exist_ok=True)
        self.adata.write_h5ad(self.output_path)

        n_assigned = sum(1 for v in cluster_types.values() if v != 'Unknown')
        n_types = len(set(cluster_types.values()) - {'Unknown'})
        self.progress_pct.emit(100)
        message = (
            f"Reference annotation complete: {n_assigned}/{len(cluster_types)} "
            f"clusters assigned ({n_types} cell types)")
        return (self.adata, message)


class MarkerScoringWorker(BaseWorker):
    """
    Worker for marker-based cell type annotation.

    Emits 'finished_ok' with a tuple '(adata, message)'.
    """

    def __init__(self, adata, marker_dict: dict, output_path: str, params: dict = None):
        super().__init__()
        self.adata = adata
        self.marker_dict = marker_dict
        self.output_path = output_path
        self.params = params or {}

    def _run(self):
        from pathlib import Path
        from kosmic.scrna.annotate.score import score_marker_genes, assign_cell_types_per_cluster

        score_threshold = self.params.get('score_threshold', 0.1)
        confidence_margin = self.params.get('confidence_margin', 0.05)
        skip_conversion = self.params.get('skip_gene_conversion', False)
        force_assignment = self.params.get('force_assignment', False)

        # Preserve author annotations before overwriting
        if 'cell_type' in self.adata.obs.columns and 'cell_type_author' not in self.adata.obs.columns:
            self.adata.obs['cell_type_author'] = self.adata.obs['cell_type'].copy()
            self.progress.emit("Preserved author annotations in 'cell_type_author' column")

        self.progress.emit("Scoring marker genes...")
        self.progress.emit(f"  Score threshold: {score_threshold}, Confidence margin: {confidence_margin}")
        self.progress_pct.emit(10)

        def _score_progress(msg, fraction=None):
            self.progress.emit(msg)
            if fraction is not None:
                self.progress_pct.emit(10 + int(fraction * 60))

        min_markers = self.params.get('min_markers', 3)

        self.adata, scored_types = score_marker_genes(
            self.adata, markers=self.marker_dict,
            skip_gene_conversion=skip_conversion,
            progress_callback=_score_progress,
            min_markers=min_markers,
        )
        self.progress_pct.emit(70)

        for ct, n_found in scored_types.items():
            total_markers = len(self.marker_dict.get(ct, []))
            if n_found > 0:
                self.progress.emit(f"Scored {ct}: {n_found}/{total_markers} markers")
            else:
                self.progress.emit(f"No markers found for {ct} (0/{total_markers})")

        skipped_types = [f"{ct} ({n})" for ct, n in scored_types.items() if n == 0]
        low_marker_types = [f"{ct} ({n})" for ct, n in scored_types.items() if 0 < n < 3]
        if skipped_types:
            self.progress.emit(f"Skipped {len(skipped_types)} types (no markers in data): {', '.join(skipped_types[:5])}")
        if low_marker_types:
            self.progress.emit(f"Warning: {len(low_marker_types)} types scored with <3 markers: {', '.join(low_marker_types[:5])}")

        # Run ORA if requested
        run_ora = self.params.get('run_ora', False)
        ora_results = {}
        if run_ora:
            self.progress.emit("Running Over-Representation Analysis (decoupler)...")
            self.progress_pct.emit(72)
            try:
                from kosmic.scrna.annotate.score import run_ora_per_cluster
                ora_results = run_ora_per_cluster(
                    self.adata, self.marker_dict,
                    skip_gene_conversion=skip_conversion,
                )
                self.progress.emit(f"ORA complete for {len(ora_results)} clusters")
            except ImportError:
                self.progress.emit("Warning: decoupler not installed, skipping ORA")
            except Exception as e:
                self.progress.emit(f"Warning: ORA failed: {e}")

        self.progress.emit("Assigning cell types per cluster...")
        self.progress_pct.emit(75)

        cluster_types, annotation_details = assign_cell_types_per_cluster(
            self.adata,
            score_threshold=score_threshold,
            confidence_margin=confidence_margin,
            force_assignment=force_assignment,
        )

        # Merge ORA p-values into annotation details
        if ora_results and annotation_details:
            for cluster_id, details in annotation_details.items():
                cluster_ora = ora_results.get(cluster_id, {})
                assigned_type = cluster_types.get(cluster_id, 'Unknown')
                details['ora_pval'] = cluster_ora.get(assigned_type)
                details['ora_all'] = cluster_ora
            self.adata.uns['cluster_ora_results'] = ora_results

        if cluster_types:
            cluster_col = 'leiden'
            for alt in ['leiden', 'clusters', 'seurat_clusters', 'louvain']:
                if alt in self.adata.obs.columns:
                    cluster_col = alt
                    break

            clusters = self.adata.obs[cluster_col]
            self.adata.obs['cell_type'] = clusters.astype(str).map(cluster_types)
            # 'Unknown' means "clustered, but no confident call". A cell
            # left out of the embedding was never clustered at all and has
            # no opinion attached to it -- calling that Unknown would look
            # like an annotation failure over what is a deliberate
            # exclusion, and puts those cells in a real category.
            unclustered = clusters.isna()
            self.adata.obs['cell_type'] = (
                self.adata.obs['cell_type'].where(unclustered,
                                                  self.adata.obs['cell_type']
                                                  .fillna('Unknown')))
            self.adata.obs.loc[unclustered, 'cell_type'] = pd.NA
            self.adata.obs['cell_type_auto'] = self.adata.obs['cell_type']

            score_cols = [c for c in self.adata.obs.columns
                          if c.endswith('_score') and c != 'cell_type_score']
            if score_cols:
                self.adata.obs['cell_type_score'] = self.adata.obs[score_cols].max(axis=1)
                # Drop per-type score columns — summary is in annotation_details
                self.adata.obs.drop(columns=score_cols, inplace=True)

            # Store marker counts per assigned type in each cluster's details
            for cluster_id, details in annotation_details.items():
                assigned = cluster_types.get(cluster_id, 'Unknown')
                total = len(self.marker_dict.get(assigned, []))
                found = scored_types.get(assigned, 0)
                details['markers_found'] = found
                details['markers_total'] = total

            self.adata.uns['cluster_annotation_details'] = annotation_details
            self.adata.uns['scored_types'] = scored_types
            self.adata.uns['annotation_params'] = {
                'score_threshold': score_threshold,
                'confidence_margin': confidence_margin,
            }

            self.progress.emit("\nCluster \u2192 Cell Type mapping:")
            for cluster, ctype in sorted(cluster_types.items(),
                                         key=lambda x: (0, int(x[0]), '') if x[0].isdigit() else (1, 0, x[0])):
                details = annotation_details.get(str(cluster), {})
                n_cells = details.get('n_cells', 0)
                best_score = details.get('best_score', 0)
                conf = details.get('confidence', '?')
                runner = details.get('runner_up', 'N/A')
                runner_score = details.get('runner_up_score', 0)
                ora_pval = details.get('ora_pval')
                ora_str = f", ORA p={ora_pval:.2e}" if ora_pval is not None else ""
                self.progress.emit(
                    f"  Cluster {cluster}: {ctype} (score={best_score:.3f}, conf={conf}{ora_str}, "
                    f"runner-up={runner} [{runner_score:.3f}], {n_cells:,} cells)"
                )

            ambiguous = [c for c, d in annotation_details.items() if d.get('confidence') == 'Ambiguous']
            if ambiguous:
                self.progress.emit(f"\n{len(ambiguous)} ambiguous cluster(s) (margin < {confidence_margin})")
        else:
            self.progress.emit("No clusters found, using per-cell assignment...")
            score_cols = [c for c in self.adata.obs.columns
                          if c.endswith('_score') and c != 'cell_type_score']
            if score_cols:
                scores = self.adata.obs[score_cols]
                self.adata.obs['cell_type'] = scores.idxmax(axis=1).str.replace('_score', '')
                self.adata.obs['cell_type_auto'] = self.adata.obs['cell_type']
                self.adata.obs['cell_type_score'] = scores.max(axis=1)

        dist = self.adata.obs['cell_type'].value_counts()
        self.progress.emit(f"\nCell type distribution:\n{dist.head(10).to_string()}")

        self.progress.emit("Saving...")
        self.progress_pct.emit(95)
        Path(self.output_path).parent.mkdir(parents=True, exist_ok=True)
        self.adata.write_h5ad(self.output_path)

        self.progress_pct.emit(100)
        cell_types = list(self.marker_dict.keys())
        message = f"Annotation complete: {len(cell_types)} cell types scored"
        return (self.adata, message)


class CellTypistWorker(BaseWorker):
    """
    Worker for CellTypist ML-based cell type annotation.

    Emits 'finished_ok' with a tuple '(adata, message)'.
    """

    def __init__(self, adata, model_name: str, majority_voting: bool, output_path: str,
                 broad_labels: bool = False):
        super().__init__()
        self.adata = adata
        self.model_name = model_name
        self.majority_voting = majority_voting
        self.output_path = output_path
        self.broad_labels = broad_labels

    def _run(self):
        try:
            from kosmic.scrna.annotate.score import run_celltypist
        except ImportError as e:
            raise RuntimeError(
                "CellTypist is not installed. Install it with:\n  pip install celltypist") from e

        # Preserve author annotations before overwriting
        if 'cell_type' in self.adata.obs.columns and 'cell_type_author' not in self.adata.obs.columns:
            self.adata.obs['cell_type_author'] = self.adata.obs['cell_type'].copy()
            self.progress.emit("Preserved author annotations in 'cell_type_author' column")

        self.progress.emit(f"Running CellTypist with model: {self.model_name}...")
        self.progress_pct.emit(10)

        def _progress(msg):
            self.progress.emit(msg)

        self.adata = run_celltypist(
            self.adata,
            self.model_name,
            majority_voting=self.majority_voting,
            progress_callback=_progress,
        )

        # Always store fine labels; optionally collapse cell_type to broad
        if 'cell_type' in self.adata.obs.columns:
            self.adata.obs['cell_type_fine'] = self.adata.obs['cell_type'].copy()
            if self.broad_labels:
                from kosmic.scrna.annotate.score import collapse_celltypist_label
                self.adata.obs['cell_type'] = (
                    self.adata.obs['cell_type_fine'].map(collapse_celltypist_label)
                )
                n_fine = self.adata.obs['cell_type_fine'].nunique()
                n_broad = self.adata.obs['cell_type'].nunique()
                self.progress.emit(f"Collapsed {n_fine} fine types → {n_broad} broad types")
            self.adata.obs['cell_type_auto'] = self.adata.obs['cell_type'].copy()

        # Build cluster annotation details for the table
        cluster_col = None
        for alt in ['leiden', 'clusters', 'seurat_clusters', 'louvain']:
            if alt in self.adata.obs.columns:
                cluster_col = alt
                break

        if cluster_col and 'cell_type' in self.adata.obs.columns:
            annotation_details = {}
            for cluster in self.adata.obs[cluster_col].unique():
                mask = self.adata.obs[cluster_col] == cluster
                cluster_cells = self.adata.obs.loc[mask]
                n_cells = int(mask.sum())

                # Get dominant type and its confidence
                type_counts = cluster_cells['cell_type'].value_counts()
                best_pct = type_counts.iloc[0] / n_cells if n_cells > 0 else 0

                # Runner-up type
                runner_up = type_counts.index[1] if len(type_counts) > 1 else ''
                runner_up_pct = type_counts.iloc[1] / n_cells if len(type_counts) > 1 else 0

                # Mean confidence score
                conf = cluster_cells['cell_type_score'].mean() if 'cell_type_score' in cluster_cells.columns else 0

                # Determine confidence level
                if best_pct > 0.8:
                    confidence = 'High'
                elif best_pct > 0.5:
                    confidence = 'Ambiguous'
                else:
                    confidence = 'Below threshold'

                annotation_details[str(cluster)] = {
                    'best_score': conf,
                    'confidence': confidence,
                    'n_cells': n_cells,
                    'runner_up': runner_up,
                    'runner_up_score': runner_up_pct,
                    'margin': best_pct - runner_up_pct,
                }

            self.adata.uns['cluster_annotation_details'] = annotation_details

        self.progress.emit("Saving...")
        self.progress_pct.emit(95)
        Path(self.output_path).parent.mkdir(parents=True, exist_ok=True)
        self.adata.write_h5ad(self.output_path)

        n_types = self.adata.obs['cell_type'].nunique()
        self.progress_pct.emit(100)
        message = f"CellTypist annotation complete: {n_types} cell types found"
        return (self.adata, message)


class RankGenesGroupsWorker(BaseWorker):
    """
    Worker for cluster-level differential expression (Wilcoxon).

    Runs once per dataset, caches result on 'adata.uns['rank_genes_groups']'
    so subsequent cluster clicks read straight from the cache. Result
    persists through h5ad writes -- no recompute on next session.

    Emits 'finished_ok' with '(adata, message)'.
    """

    def __init__(self, adata, cluster_col: str = 'leiden'):
        super().__init__()
        self.adata = adata
        self.cluster_col = cluster_col

    def _run(self):
        from kosmic.scrna.annotate.marker_genes import compute_cluster_marker_genes

        n_cells = self.adata.n_obs
        n_clusters = self.adata.obs[self.cluster_col].nunique()
        self.progress.emit(
            f"Computing top genes per cluster ({n_clusters} clusters, "
            f"{n_cells:,} cells, Wilcoxon)...")
        self.progress_pct.emit(20)

        compute_cluster_marker_genes(
            self.adata, cluster_col=self.cluster_col, method='wilcoxon')

        self.progress_pct.emit(100)
        return (self.adata, f"Top genes computed for {n_clusters} clusters")


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CELLTYPIST_MODELS = {
    "Pan-Immune": [
        ("Immune_All_Low.pkl", "98 types \u2014 immune sub-populations, 20 tissues"),
        ("Immune_All_High.pkl", "32 types \u2014 immune populations, 20 tissues"),
    ],
    "Heart": [
        ("Healthy_Adult_Heart.pkl", "75 types \u2014 healthy adult human heart"),
    ],
    "Brain": [
        ("Adult_Human_MTG.pkl", "127 types \u2014 adult human middle temporal gyrus"),
        ("Adult_Human_PrefrontalCortex.pkl", "113 types \u2014 adult human prefrontal cortex"),
        ("Developing_Human_Brain.pkl", "129 types \u2014 first-trimester developing human brain"),
        ("Developing_Mouse_Brain.pkl", "174 types \u2014 embryonic mouse brain"),
        ("Mouse_Whole_Brain.pkl", "334 types \u2014 whole adult mouse brain"),
        ("Mouse_Isocortex_Hippocampus.pkl", "42 types \u2014 adult mouse isocortex & hippocampus"),
        ("Human_AdultAged_Hippocampus.pkl", "15 types \u2014 adult/aged human hippocampus"),
        ("Human_Longitudinal_Hippocampus.pkl", "24 types \u2014 anterior/posterior hippocampus"),
        ("Developing_Human_Hippocampus.pkl", "15 types \u2014 developing human hippocampus"),
        ("Adult_CynomolgusMacaque_Hippocampus.pkl", "34 types \u2014 cynomolgus macaque hippocampus"),
        ("Adult_RhesusMacaque_Hippocampus.pkl", "14 types \u2014 rhesus macaque hippocampus"),
        ("Adult_Pig_Hippocampus.pkl", "13 types \u2014 adult pig hippocampus"),
        ("Mouse_Dentate_Gyrus.pkl", "24 types \u2014 perinatal/juvenile/adult mouse"),
        ("Mouse_Postnatal_DentateGyrus.pkl", "22 types \u2014 postnatal mouse dentate gyrus"),
        ("Developing_Mouse_Hippocampus.pkl", "11 types \u2014 postnatal day 7 mouse"),
        ("Adult_Mouse_OlfactoryBulb.pkl", "40 types \u2014 adult mouse olfactory bulb"),
    ],
    "Lung & Airways": [
        ("Human_Lung_Atlas.pkl", "61 types \u2014 integrated Human Lung Cell Atlas"),
        ("Cells_Lung_Airway.pkl", "78 types \u2014 human lungs/airways (scRNA)"),
        ("Nuclei_Lung_Airway.pkl", "78 types \u2014 human lungs/airways (snRNA)"),
        ("Cells_Fetal_Lung.pkl", "144 types \u2014 human embryonic/fetal lungs"),
        ("Human_IPF_Lung.pkl", "38 types \u2014 IPF, COPD, and healthy lungs"),
        ("Human_PF_Lung.pkl", "31 types \u2014 pulmonary fibrosis lungs"),
        ("Autopsy_COVID19_Lung.pkl", "34 types \u2014 COVID-19 autopsy lungs"),
        ("Lethal_COVID19_Lung.pkl", "41 types \u2014 lethal COVID-19 lungs"),
        ("COVID19_Immune_Landscape.pkl", "64 types \u2014 lung & blood COVID-19 immune"),
        ("PaediatricAdult_COVID19_Airway.pkl", "59 types \u2014 COVID-19 airways"),
    ],
    "Blood & PBMC": [
        ("COVID19_HumanChallenge_Blood.pkl", "83 types \u2014 SARS-CoV-2 challenge blood"),
        ("Healthy_COVID19_PBMC.pkl", "51 types \u2014 healthy & COVID-19 PBMC"),
        ("Adult_COVID19_PBMC.pkl", "20 types \u2014 COVID-19 PBMC"),
        ("PaediatricAdult_COVID19_PBMC.pkl", "42 types \u2014 paediatric/adult COVID-19 PBMC"),
        ("Adult_cHSPCs_Illumina.pkl", "20 types \u2014 circulating HSPCs (Illumina)"),
        ("Adult_cHSPCs_Ultima.pkl", "20 types \u2014 circulating HSPCs (Ultima)"),
    ],
    "Gut & Intestinal": [
        ("Cells_Intestinal_Tract.pkl", "134 types \u2014 fetal/pediatric/adult human gut"),
        ("Adult_Mouse_Gut.pkl", "126 types \u2014 adult mouse gut"),
        ("Human_Colorectal_Cancer.pkl", "36 types \u2014 colorectal cancer colon tissue"),
    ],
    "Liver": [
        ("Healthy_Human_Liver.pkl", "17 types \u2014 adult human liver"),
        ("Healthy_Mouse_Liver.pkl", "17 types \u2014 healthy murine liver"),
    ],
    "Skin": [
        ("Adult_Human_Skin.pkl", "34 types \u2014 healthy adult skin"),
        ("Fetal_Human_Skin.pkl", "14 types \u2014 developing fetal skin"),
    ],
    "Pancreas": [
        ("Adult_Human_PancreaticIslet.pkl", "12 types \u2014 healthy adult pancreatic islets"),
        ("Fetal_Human_Pancreas.pkl", "19 types \u2014 embryonic pancreas 9-19 weeks"),
    ],
    "Vascular": [
        ("Adult_Human_Vascular.pkl", "42 types \u2014 vascular populations multi-organ"),
    ],
    "Breast": [
        ("Cells_Adult_Breast.pkl", "58 types \u2014 adult human breast"),
    ],
    "Reproductive & Placental": [
        ("Human_Endometrium_Atlas.pkl", "36 types \u2014 endometrial, menstrual cycle"),
        ("Human_Placenta_Decidua.pkl", "32 types \u2014 first-trimester placenta/decidua"),
        ("Developing_Human_Gonads.pkl", "53 types \u2014 gonadal/extragonadal tissue"),
    ],
    "Thymus & Tonsil": [
        ("Developing_Human_Thymus.pkl", "55 types \u2014 embryonic to adult thymus"),
        ("Cells_Human_Tonsil.pkl", "117 types \u2014 human tonsil 3-65 years"),
    ],
    "Eye": [
        ("Fetal_Human_Retina.pkl", "25 types \u2014 fetal neural retina"),
        ("Human_Developmental_Retina.pkl", "16 types \u2014 fetal retina"),
    ],
    "Other Developmental": [
        ("Pan_Fetal_Human.pkl", "138 types \u2014 stromal/immune from human fetus"),
        ("Developing_Human_Organs.pkl", "27 types \u2014 5 endoderm organs 7-21 weeks"),
        ("Fetal_Human_AdrenalGlands.pkl", "9 types \u2014 fetal adrenal glands"),
        ("Fetal_Human_Pituitary.pkl", "14 types \u2014 fetal pituitary"),
        ("Human_Embryonic_YolkSac.pkl", "49 types \u2014 yolk sac 4-8 weeks"),
        ("Nuclei_Human_InnerEar.pkl", "18 types \u2014 fetal/adult inner ear"),
    ],
    "Mouse Immune": [
        ("Mouse_Dendritic_Subtypes.pkl", "10 types \u2014 dendritic cells multi-tissue"),
    ],
}

DEFAULT_MARKERS = {
    'Endothelial': ['PECAM1', 'CDH5', 'VWF', 'KDR', 'FLT1', 'CLDN5', 'ERG'],
    'Cardiomyocyte': ['TNNT2', 'MYH7', 'MYH6', 'ACTC1', 'TTN', 'MYBPC3'],
    'Fibroblast': ['DCN', 'COL1A1', 'COL1A2', 'LUM', 'PDGFRA', 'THY1'],
    'Macrophage': ['CD68', 'CD14', 'FCGR3A', 'CSF1R', 'AIF1', 'MARCO'],
    'T_cell': ['CD3D', 'CD3E', 'CD4', 'CD8A', 'IL7R', 'TRAC'],
    'B_cell': ['CD19', 'MS4A1', 'CD79A', 'CD79B', 'PAX5'],
    'Smooth_Muscle': ['ACTA2', 'MYH11', 'TAGLN', 'CNN1', 'MYOCD'],
    'Pericyte': ['PDGFRB', 'RGS5', 'KCNJ8', 'ABCC9', 'NOTCH3'],
    'Neuron': ['MAP2', 'RBFOX3', 'SYP', 'SNAP25', 'TUBB3'],
    'Astrocyte': ['GFAP', 'AQP4', 'S100B', 'SLC1A2', 'ALDH1L1'],
    'Microglia': ['CX3CR1', 'P2RY12', 'TMEM119', 'AIF1', 'ITGAM'],
    'Oligodendrocyte': ['MBP', 'MOG', 'PLP1', 'OLIG2', 'MAG'],
}

# ---------------------------------------------------------------------------
# AnnotateTab
# ---------------------------------------------------------------------------

class AnnotateTab(QWidget):
    """Tab for cell type annotation and cluster re-labelling."""

    log_message = pyqtSignal(str)
    annotation_complete = pyqtSignal(str)  # emitted with status message

    def __init__(self, main_window, embedded=True):
        # AnnotateTab is structurally embedded: ClusterTab mounts
        # 'controls_widget' in its sidebar and treats 'self' as the
        # right-side panel content. The 'embedded' arg is unused in
        # the body; it exists so callers can stay explicit about the
        # mount mode.
        super().__init__()
        self.main_window = main_window
        self.project_dir = None
        self.adata = None
        self.h5ad_path = None
        self.marker_worker = None
        self.celltypist_worker = None
        self.ref_worker = None

        # Cell-type focus presets (Cardiac, Stroke, plus any user-added
        # JSON in '~/.kosmic/cell_type_focus/'). Loaded once on tab init.
        self._focus_presets = load_focus_presets()

        self._setup_ui()
        self._restore_last_preset()

    def reset_state(self):
        """Clear all cached state (called on raw data reload)."""
        # Kill running workers
        for attr in ('marker_worker', 'celltypist_worker', 'ref_worker'):
            w = getattr(self, attr, None)
            if w is not None and w.isRunning():
                w.quit()
                w.wait(2000)
        self.marker_worker = None
        self.celltypist_worker = None
        self.ref_worker = None
        self.adata = None
        self.h5ad_path = None
        self._last_adata_version = -1
        # Reset UI
        self.annot_status.setText("")
        self.cluster_table.setRowCount(0)

    def _setup_ui(self):
        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(0, 0, 0, 0)

        # --- Annotation method controls (mounted in ClusterTab's sidebar) ---
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(4, 4, 4, 4)
        left_layout.setSpacing(8)

        # Status
        self.annot_status = StatusLabel("")
        self.annot_status.setWordWrap(True)
        left_layout.addWidget(self.annot_status)

        self.use_existing_check = QCheckBox("Use existing cell_type")
        self.use_existing_check.setChecked(True)
        left_layout.addWidget(self.use_existing_check)

        left_layout.addWidget(QLabel("Method"))
        self.method_combo = NoScrollComboBox()
        self.method_combo.addItem("CellTypist")
        self.method_combo.addItem("Reference Atlas")
        self.method_combo.addItem("CellMarker 2.0")
        self.method_combo.addItem("PanglaoDB")
        self.method_combo.currentIndexChanged.connect(self._on_method_changed)
        left_layout.addWidget(self.method_combo)

        # Stacked widget for method-specific panels
        self.annotation_stack = QStackedWidget()
        self.annotation_stack.currentChanged.connect(self._resize_stack)

        # --- Page 0: CellTypist panel ---
        celltypist_page = QWidget()
        ct_layout = QVBoxLayout(celltypist_page)
        ct_layout.setContentsMargins(0, 4, 0, 0)

        ct_layout.addWidget(QLabel("Model"))
        self.model_combo = NoScrollComboBox()
        self.model_combo.setMaxVisibleItems(20)
        self.model_combo.setMinimumContentsLength(8)
        self.model_combo.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed)
        self._populate_model_combo()
        ct_layout.addWidget(self.model_combo)

        self.majority_voting_check = QCheckBox("Majority voting")
        self.majority_voting_check.setChecked(True)
        self.majority_voting_check.setToolTip("Consensus cell type label per cluster.")
        ct_layout.addWidget(self.majority_voting_check)

        self.broad_labels_check = QCheckBox("Broad labels")
        self.broad_labels_check.setChecked(True)
        self.broad_labels_check.setToolTip("Collapse fine-grained labels to broad categories.")
        ct_layout.addWidget(self.broad_labels_check)

        self.run_celltypist_btn = QPushButton("Run CellTypist")
        self.run_celltypist_btn.clicked.connect(self._run_celltypist)
        self.run_celltypist_btn.setEnabled(False)
        ct_layout.addWidget(self.run_celltypist_btn)

        ct_layout.addStretch()
        self.annotation_stack.addWidget(celltypist_page)

        # --- Page 1: Reference Atlas panel ---
        ref_page = QWidget()
        ref_layout = QVBoxLayout(ref_page)
        ref_layout.setContentsMargins(0, 4, 0, 0)

        ref_layout.addWidget(QLabel("Reference"))
        self.ref_combo = NoScrollComboBox()
        self.ref_combo.setMinimumContentsLength(8)
        self.ref_combo.setMaxVisibleItems(15)
        self._populate_ref_combo()
        ref_layout.addWidget(self.ref_combo)

        self.ref_info_label = HintLabel("")
        self.ref_info_label.setWordWrap(True)
        ref_layout.addWidget(self.ref_info_label)
        self.ref_combo.currentIndexChanged.connect(self._on_ref_changed)
        self._on_ref_changed()

        ref_quality_group = SettingsGroup("Settings", collapsible=True, expanded=True)
        self.ref_min_corr_spin = NoScrollDoubleSpinBox()
        self.ref_min_corr_spin.setRange(0.0, 1.0)
        self.ref_min_corr_spin.setValue(0.3)
        self.ref_min_corr_spin.setSingleStep(0.05)
        self.ref_min_corr_spin.setToolTip(
            "Minimum Pearson correlation to assign a cell type.\n"
            "Clusters below this threshold are labelled 'Unknown'.\n"
            "Typical range: 0.2-0.5"
        )
        ref_quality_group.add_row("Min correlation:", self.ref_min_corr_spin)
        ref_layout.addWidget(ref_quality_group)

        self.run_ref_btn = QPushButton("Run Reference")
        self.run_ref_btn.clicked.connect(self._run_reference_annotation)
        self.run_ref_btn.setEnabled(False)
        ref_layout.addWidget(self.run_ref_btn)

        self.build_ref_btn = QPushButton("Build from h5ad...")
        self.build_ref_btn.setToolTip(
            "Build a lightweight reference from an annotated h5ad file.\n"
            "The atlas must have a cell type column in obs."
        )
        self.build_ref_btn.clicked.connect(self._build_custom_reference)
        ref_layout.addWidget(self.build_ref_btn)

        ref_layout.addStretch()
        self.annotation_stack.addWidget(ref_page)

        # --- Page 2: CellMarker 2.0 panel ---
        cm2_page = QWidget()
        cm2_layout = QVBoxLayout(cm2_page)
        cm2_layout.setContentsMargins(0, 4, 0, 0)

        cm2_species_layout = QHBoxLayout()
        cm2_species_layout.addWidget(QLabel("Species:"))
        self.cm2_species_combo = NoScrollComboBox()
        self.cm2_species_combo.addItems(["Human", "Mouse"])
        self.cm2_species_combo.setToolTip("CellMarker 2.0 markers are species-specific.\nSelect the species matching your dataset.")
        self.cm2_species_combo.currentTextChanged.connect(self._on_cm2_filter_changed)
        cm2_species_layout.addWidget(self.cm2_species_combo, 1)
        cm2_layout.addLayout(cm2_species_layout)

        cm2_tissue_layout = QHBoxLayout()
        cm2_tissue_layout.addWidget(QLabel("Tissue:"))
        self.cm2_tissue_combo = NoScrollComboBox()
        self.cm2_tissue_combo.setMinimumContentsLength(10)
        self.cm2_tissue_combo.addItem("All")
        self._populate_cm2_tissue_combo()
        self.cm2_tissue_combo.currentTextChanged.connect(self._on_cm2_filter_changed)
        cm2_tissue_layout.addWidget(self.cm2_tissue_combo, 1)
        cm2_layout.addLayout(cm2_tissue_layout)

        cm2_cancer_layout = QHBoxLayout()
        cm2_cancer_layout.addWidget(QLabel("Cancer type:"))
        self.cm2_cancer_combo = NoScrollComboBox()
        self.cm2_cancer_combo.addItems(["Normal", "All"])
        self.cm2_cancer_combo.setToolTip("'Normal' excludes cancer-specific markers.\n'All' includes everything.")
        cm2_cancer_layout.addWidget(self.cm2_cancer_combo, 1)
        cm2_layout.addLayout(cm2_cancer_layout)

        cm2_layout.addWidget(QLabel("Cell types"))
        self.cm2_cell_list = QListWidget()
        self.cm2_cell_list.setSelectionMode(QListWidget.SelectionMode.MultiSelection)
        self.cm2_cell_list.setMaximumHeight(250)
        self._populate_cm2_cell_types()
        cm2_layout.addWidget(self.cm2_cell_list)

        cm2_select_grid = QGridLayout()
        cm2_select_grid.setSpacing(4)
        self.cm2_select_all_btn = QPushButton("Select All")
        self.cm2_select_all_btn.clicked.connect(self._cm2_select_all)
        cm2_select_grid.addWidget(self.cm2_select_all_btn, 0, 0)
        self.cm2_select_none_btn = QPushButton("Select None")
        self.cm2_select_none_btn.clicked.connect(self._cm2_select_none)
        cm2_select_grid.addWidget(self.cm2_select_none_btn, 0, 1)
        self.cm2_preset_combo = NoScrollComboBox()
        self.cm2_preset_combo.addItem("— Apply preset —")
        for name, preset in self._focus_presets.items():
            if preset.get('cellmarker2_human') or preset.get('cellmarker2_mouse'):
                self.cm2_preset_combo.addItem(name)
        self.cm2_preset_combo.setToolTip(
            "Select cell types from a tissue-focus preset (uses the "
            "current Species). Bundled: Cardiac. Drop your own JSON "
            "in ~/.kosmic/cell_type_focus/ to add more.")
        self.cm2_preset_combo.currentIndexChanged.connect(
            self._on_cm2_preset_changed)
        cm2_select_grid.addWidget(self.cm2_preset_combo, 1, 0, 1, 2)
        cm2_layout.addLayout(cm2_select_grid)

        cm2_quality_group = SettingsGroup("Quality Settings", collapsible=True, expanded=True)
        self.cm2_score_threshold_spin = NoScrollDoubleSpinBox()
        self.cm2_score_threshold_spin.setRange(-1.0, 1.0)
        self.cm2_score_threshold_spin.setValue(0.0)
        self.cm2_score_threshold_spin.setSingleStep(0.05)
        self.cm2_score_threshold_spin.setToolTip("Minimum UCell score (0–1) to assign a cell type.\nClusters below this are labelled 'Unknown'.\n0.0 recommended (use confidence margin for quality control).")
        cm2_quality_group.add_row("Score threshold:", self.cm2_score_threshold_spin)

        self.cm2_confidence_margin_spin = NoScrollDoubleSpinBox()
        self.cm2_confidence_margin_spin.setRange(0.0, 0.5)
        self.cm2_confidence_margin_spin.setValue(0.05)
        self.cm2_confidence_margin_spin.setSingleStep(0.01)
        self.cm2_confidence_margin_spin.setToolTip("Minimum gap between 1st and 2nd best score.\nIf gap < margin, cluster is flagged as 'Ambiguous'.")
        cm2_quality_group.add_row("Confidence margin:", self.cm2_confidence_margin_spin)

        self.cm2_min_markers_spin = QSpinBox()
        self.cm2_min_markers_spin.setRange(1, 20)
        self.cm2_min_markers_spin.setValue(3)
        self.cm2_min_markers_spin.setToolTip("Minimum markers found in data to score a cell type.\nTypes with fewer markers are skipped (unreliable).")
        cm2_quality_group.add_row("Min markers:", self.cm2_min_markers_spin)

        self.cm2_ora_check = QCheckBox("Run ORA")
        self.cm2_ora_check.setToolTip("Run decoupler ORA alongside marker scoring.\nAdds statistical p-values to the cluster table.")
        cm2_quality_group.add_widget(self.cm2_ora_check)

        cm2_layout.addWidget(cm2_quality_group)

        self.run_cm2_btn = QPushButton("Run CellMarker 2.0")
        self.run_cm2_btn.clicked.connect(self._run_cellmarker2_annotation)
        self.run_cm2_btn.setEnabled(False)
        cm2_layout.addWidget(self.run_cm2_btn)

        cm2_layout.addStretch()
        self.annotation_stack.addWidget(cm2_page)

        # --- Page 2: PanglaoDB panel ---
        panglaodb_page = QWidget()
        pdb_layout = QVBoxLayout(panglaodb_page)
        pdb_layout.setContentsMargins(0, 4, 0, 0)

        tissue_layout = QHBoxLayout()
        tissue_layout.addWidget(QLabel("Tissue/Organ:"))
        self.tissue_combo = NoScrollComboBox()
        self.tissue_combo.addItem("All")
        self._populate_tissue_combo()
        self.tissue_combo.currentTextChanged.connect(self._on_tissue_changed)
        tissue_layout.addWidget(self.tissue_combo, 1)
        pdb_layout.addLayout(tissue_layout)

        pdb_layout.addWidget(QLabel("Cell types"))
        self.marker_list = QListWidget()
        self.marker_list.setSelectionMode(QListWidget.SelectionMode.MultiSelection)
        self.marker_list.setMaximumHeight(250)
        self._populate_cell_types()
        pdb_layout.addWidget(self.marker_list)

        select_grid = QGridLayout()
        select_grid.setSpacing(4)
        self.select_all_btn = QPushButton("Select All")
        self.select_all_btn.clicked.connect(self._select_all_cell_types)
        select_grid.addWidget(self.select_all_btn, 0, 0)
        self.select_none_btn = QPushButton("Select None")
        self.select_none_btn.clicked.connect(self._select_no_cell_types)
        select_grid.addWidget(self.select_none_btn, 0, 1)
        self.pdb_preset_combo = NoScrollComboBox()
        self.pdb_preset_combo.addItem("— Apply preset —")
        for name, preset in self._focus_presets.items():
            if preset.get('panglaodb'):
                self.pdb_preset_combo.addItem(name)
        self.pdb_preset_combo.setToolTip(
            "Select cell types from a tissue-focus preset. Bundled: "
            "Cardiac, Stroke. Drop your own JSON in "
            "~/.kosmic/cell_type_focus/ to add more.")
        self.pdb_preset_combo.currentIndexChanged.connect(
            self._on_pdb_preset_changed)
        select_grid.addWidget(self.pdb_preset_combo, 1, 0, 1, 2)
        pdb_layout.addLayout(select_grid)

        quality_group = SettingsGroup("Quality Settings", collapsible=True, expanded=True)
        self.species_combo = NoScrollComboBox()
        self.species_combo.addItems(["Both", "Human", "Mouse"])
        self.species_combo.setToolTip("Filter PanglaoDB markers by species.\n'Both' uses average specificity across human and mouse.")
        quality_group.add_row("Species:", self.species_combo)

        self.specificity_spin = NoScrollDoubleSpinBox()
        self.specificity_spin.setRange(0.0, 1.0)
        self.specificity_spin.setValue(0.0)
        self.specificity_spin.setSingleStep(0.05)
        self.specificity_spin.setToolTip("Minimum marker specificity (0-1).\nHigher values = fewer but more specific markers.\nStart at 0, increase to 0.1-0.3 if seeing misannotations.")
        quality_group.add_row("Min specificity:", self.specificity_spin)

        self.score_threshold_spin = NoScrollDoubleSpinBox()
        self.score_threshold_spin.setRange(0.0, 1.0)
        self.score_threshold_spin.setValue(0.1)
        self.score_threshold_spin.setSingleStep(0.05)
        self.score_threshold_spin.setToolTip("Minimum UCell score (0–1) to assign a cell type.\nClusters below this are labelled 'Unknown'.")
        quality_group.add_row("Score threshold:", self.score_threshold_spin)

        self.confidence_margin_spin = NoScrollDoubleSpinBox()
        self.confidence_margin_spin.setRange(0.0, 0.5)
        self.confidence_margin_spin.setValue(0.05)
        self.confidence_margin_spin.setSingleStep(0.01)
        self.confidence_margin_spin.setToolTip("Minimum gap between 1st and 2nd best score.\nIf gap < margin, cluster is flagged as 'Ambiguous'.")
        quality_group.add_row("Confidence margin:", self.confidence_margin_spin)

        self.pdb_min_markers_spin = QSpinBox()
        self.pdb_min_markers_spin.setRange(1, 20)
        self.pdb_min_markers_spin.setValue(3)
        self.pdb_min_markers_spin.setToolTip("Minimum markers found in data to score a cell type.\nTypes with fewer markers are skipped (unreliable).")
        quality_group.add_row("Min markers:", self.pdb_min_markers_spin)

        self.canonical_check = QCheckBox("Canonical markers only")
        self.canonical_check.setToolTip("Use only canonical PanglaoDB markers.\nFewer markers but higher confidence.")
        quality_group.add_widget(self.canonical_check)

        self.pdb_ora_check = QCheckBox("Run ORA")
        self.pdb_ora_check.setToolTip("Run decoupler ORA alongside marker scoring.\nAdds statistical p-values to the cluster table.")
        quality_group.add_widget(self.pdb_ora_check)

        pdb_layout.addWidget(quality_group)

        self.run_annot_btn = QPushButton("Run PanglaoDB")
        self.run_annot_btn.clicked.connect(self._run_annotation)
        self.run_annot_btn.setEnabled(False)
        pdb_layout.addWidget(self.run_annot_btn)

        pdb_layout.addStretch()
        self.annotation_stack.addWidget(panglaodb_page)

        left_layout.addWidget(self.annotation_stack)
        self._resize_stack(0)  # size to first page

        # ClusterTab mounts this in its own sidebar.
        self.controls_widget = left_widget

        # --- Right: Cluster summary table ---
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)

        summary_group = SettingsGroup("Cluster Summary")
        summary_layout = QVBoxLayout()
        summary_layout.setContentsMargins(0, 0, 0, 0)
        summary_layout.setSpacing(6)
        summary_group.add_layout(summary_layout)

        self.cluster_table = ResultsTable()
        self.cluster_table.set_schema([
            Column("Cluster",      "cluster",     "s"),
            Column("Cells",        "count",       ",d"),
            Column("Author Type",  "author",      "s"),
            Column("Assigned Type", "assigned",   "s"),
            Column("Score",        "score_text",  "s"),
            Column("Fine Type",    "detail_text", "s"),
            Column("ORA p-val",    "ora_text",    "s", na_text="—"),
            Column("Confidence",   "conf",        "s"),
            Column("Runner-up",    "runner_text", "s"),
            Column("Flags",        "flags_text",  "s"),
        ])
        header = self.cluster_table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        header.setStretchLastSection(True)
        self.cluster_table.setHorizontalScrollMode(QTableWidget.ScrollMode.ScrollPerPixel)
        self.cluster_table.setAlternatingRowColors(True)
        self.cluster_table.setToolTip("Right-click for options, double-click 'Assigned Type' to rename")
        self.cluster_table.cellDoubleClicked.connect(self._on_cluster_table_double_click)
        self.cluster_table.cellChanged.connect(self._on_cluster_table_cell_changed)
        # Single click on any cell -> populate the top-genes panel below
        self.cluster_table.cellClicked.connect(self._on_cluster_row_clicked)
        self.cluster_table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.cluster_table.customContextMenuRequested.connect(self._on_cluster_table_context_menu)
        # How well the clusters line up with the labels the authors
        # deposited, when the file carries any (cell_type_author).
        self.author_agreement_label = HintLabel("")
        self.author_agreement_label.setWordWrap(True)
        self.author_agreement_label.setVisible(False)
        summary_layout.addWidget(self.author_agreement_label)
        summary_layout.addWidget(self.cluster_table)

        hint_label = HintLabel(
            "Right-click a row to rename, merge, or mark clusters as Unknown")
        summary_layout.addWidget(hint_label)

        btn_layout = QHBoxLayout()
        self.save_annot_btn = QPushButton("Save Changes")
        self.save_annot_btn.setToolTip("Save current annotations to h5ad file")
        self.save_annot_btn.clicked.connect(self._save_annotations)
        self.save_annot_btn.setEnabled(False)
        btn_layout.addWidget(self.save_annot_btn)

        self.undo_annot_btn = QPushButton("Undo Changes")
        self.undo_annot_btn.setToolTip("Revert to the last automated annotation")
        self.undo_annot_btn.clicked.connect(self._undo_annotation_changes)
        self.undo_annot_btn.setEnabled(False)
        btn_layout.addWidget(self.undo_annot_btn)
        summary_layout.addLayout(btn_layout)

        top_genes_group = SettingsGroup("Top Differential Genes (selected cluster)")
        top_genes_layout = QVBoxLayout()
        top_genes_layout.setContentsMargins(0, 0, 0, 0)
        self._setup_top_genes_panel(top_genes_layout)
        top_genes_group.add_layout(top_genes_layout)

        right_v_splitter = QSplitter(Qt.Orientation.Vertical)
        right_v_splitter.addWidget(summary_group)
        right_v_splitter.addWidget(top_genes_group)
        right_v_splitter.setStretchFactor(0, 2)
        right_v_splitter.setStretchFactor(1, 1)
        right_v_splitter.setChildrenCollapsible(False)
        right_layout.addWidget(right_v_splitter)

        outer_layout.addWidget(right_widget)

        self.progress_bar = None  # wired to sidebar by main.py
        self.progress_label = QLabel("Ready")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_project_directory(self, directory: str):
        self.project_dir = Path(directory)

    def set_adata(self, adata):
        """Set the AnnData object from main window."""
        self.adata = adata
        self._update_status()

    def on_tab_activated(self):
        """Called when tab becomes visible — pull latest adata from workspace."""
        ws = self.main_window
        if ws.current_adata is None:
            return
        version = getattr(ws, '_adata_version', 0)
        if version != getattr(self, '_last_adata_version', -1):
            self._last_adata_version = version
            self.set_adata(ws.current_adata)

    # ------------------------------------------------------------------
    # Status
    # ------------------------------------------------------------------

    def _on_status(self, status: str):
        self.progress_label.setText(status)
        self.log_message.emit(status)

    def _on_progress(self, value: int):
        if self.progress_bar:
            self.progress_bar.setValue(value)

    def _update_status(self):
        """Update button enable/disable based on current data."""
        if self.adata is None:
            self.annot_status.setText("No data loaded. Run clustering first.")
            self.annot_status.set_state('warning')
            self.run_annot_btn.setEnabled(False)
            self.run_celltypist_btn.setEnabled(False)
            self.run_cm2_btn.setEnabled(False)
            self.run_ref_btn.setEnabled(False)
            return

        has_leiden = 'leiden' in self.adata.obs.columns
        has_clusters = 'clusters' in self.adata.obs.columns or 'seurat_clusters' in self.adata.obs.columns
        has_cell_type = 'cell_type' in self.adata.obs.columns

        can_annotate = has_leiden or has_clusters

        if has_cell_type:
            n_types = self.adata.obs['cell_type'].nunique()
            self.annot_status.setText(f"cell_type column exists ({n_types} types)")
            self.annot_status.set_state('success')
        elif can_annotate:
            self.annot_status.setText("Clusters found \u2014 ready for annotation")
            self.annot_status.set_state('info')
        else:
            self.annot_status.setText("No clusters found \u2014 run clustering first")
            self.annot_status.set_state('warning')

        self.run_annot_btn.setEnabled(can_annotate)
        self.run_cm2_btn.setEnabled(can_annotate)
        self.run_ref_btn.setEnabled(can_annotate)
        self.run_celltypist_btn.setEnabled(True)  # CellTypist works without clusters

        if can_annotate or has_cell_type:
            self._populate_cluster_table()

    # ------------------------------------------------------------------
    # Cluster table
    # ------------------------------------------------------------------

    def _populate_cluster_table(self):
        """Populate the cluster summary table with annotation diagnostics."""
        self._table_updating = True

        if self.adata is None:
            self.cluster_table.set_data(pd.DataFrame())
            self._table_updating = False
            return

        if 'leiden' in self.adata.obs.columns:
            cluster_col = 'leiden'
        elif 'clusters' in self.adata.obs.columns:
            cluster_col = 'clusters'
        elif 'seurat_clusters' in self.adata.obs.columns:
            cluster_col = 'seurat_clusters'
        else:
            self.cluster_table.set_data(pd.DataFrame())
            self._table_updating = False
            return

        if 'cell_type' in self.adata.obs.columns:
            ct_col = 'cell_type'
        elif 'cell_type_auto' in self.adata.obs.columns:
            ct_col = 'cell_type_auto'
        elif 'Names' in self.adata.obs.columns:
            ct_col = 'Names'
        else:
            ct_col = None

        author_col = ('cell_type_author'
                      if 'cell_type_author' in self.adata.obs.columns
                      else None)
        has_fine = 'cell_type_fine' in self.adata.obs.columns

        details = self.adata.uns.get('cluster_annotation_details', {})
        cluster_counts = (
            self.adata.obs[cluster_col].value_counts().sort_index())

        # Cluster-level quality signals (weak match, MT markers, doublet
        # score, authors' unlabelled cells): the per-cell QC cannot see
        # these, and every one of them was reported 'High' before.
        try:
            flags = flag_clusters(
                self.adata, cluster_col=cluster_col,
                type_col=ct_col or 'cell_type', author_col=author_col,
                details=details)
        except Exception as exc:  # flags are advisory; never block the table
            self._on_status(f"Cluster flags unavailable: {exc}")
            flags = {}

        author_labels = None
        if author_col:
            author_labels = author_label_series(self.adata.obs[author_col])
        n_labelled = n_agree = n_unlabelled = 0

        rows = []
        for cluster, count in cluster_counts.items():
            cluster_mask = self.adata.obs[cluster_col] == cluster

            author_text = ""
            author_tip = ""
            if author_labels is not None and cluster_mask.sum() > 0:
                breakdown = author_breakdown(author_labels[cluster_mask.values])
                if breakdown.n_labelled:
                    top_label, top_n = breakdown.top
                    author_text = (f"{top_label} "
                                   f"({100 * top_n / breakdown.n_labelled:.0f}%)")
                    author_tip = breakdown.describe()
                    n_labelled += breakdown.n_labelled
                    n_agree += top_n
                n_unlabelled += breakdown.n_unlabelled

            ct_text = ""
            if ct_col and cluster_mask.sum() > 0:
                dominant_ct = self.adata.obs.loc[cluster_mask, ct_col].mode()
                if len(dominant_ct) > 0:
                    ct_text = str(dominant_ct.iloc[0])

            cluster_details = details.get(str(cluster), {})

            score_val = cluster_details.get('best_score', '')
            score_text = (f"{float(score_val):.3f}"
                          if score_val != '' and np.isreal(score_val)
                          else "")

            m_found = cluster_details.get('markers_found')
            m_total = cluster_details.get('markers_total')
            if has_fine:
                fine_types = self.adata.obs.loc[
                    cluster_mask, 'cell_type_fine'].value_counts()
                detail_text = (fine_types.index[0]
                               if len(fine_types) > 0 else "")
            elif m_found is not None and m_total is not None:
                detail_text = f"{m_found}/{m_total} markers"
            else:
                detail_text = ""

            ora_pval = cluster_details.get('ora_pval')
            ora_text = f"{ora_pval:.2e}" if ora_pval is not None else "—"

            conf_text = cluster_details.get('confidence', '')

            runner_up = cluster_details.get('runner_up', '')
            runner_score = cluster_details.get('runner_up_score', '')
            runner_text = (f"{runner_up} ({float(runner_score):.3f})"
                           if runner_up and np.isreal(runner_score)
                           else "")

            rows.append({
                'cluster':     str(cluster),
                'count':       int(count),
                'author':      author_text,
                '_author_tip': author_tip,
                'assigned':    ct_text,
                'score_text':  score_text,
                'detail_text': detail_text,
                'ora_text':    ora_text,
                'conf':        conf_text,
                'runner_text': runner_text,
                'flags_text':  (flags[str(cluster)].text()
                                if str(cluster) in flags else ""),
                # Row-level flags for decorate (hidden from rendering)
                '_m_found':    m_found if m_found is not None else -1,
                '_has_fine':   has_fine,
                '_ora_pval':   (float(ora_pval)
                                if ora_pval is not None
                                else float('nan')),
            })

        df = pd.DataFrame(rows)

        yellow = Qt.GlobalColor.yellow
        green = Qt.GlobalColor.green
        red = Qt.GlobalColor.red

        def _decorate(item, row, col_i):
            if col_i == 2:
                if row['_author_tip']:
                    item.setToolTip(row['_author_tip'])
            elif col_i == 3:
                # Assigned Type is user-editable (context-menu renames and
                # double-click cell edits both feed _on_cluster_table_cell_changed).
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsEditable)
            elif col_i == 5:
                if (not row['_has_fine']
                        and 0 <= row['_m_found'] < 3):
                    item.setForeground(yellow)
                    item.setToolTip(
                        "Few markers — assignment may be unreliable")
            elif col_i == 6:
                p = row['_ora_pval']
                if p == p:  # not NaN
                    item.setForeground(green if p < DEFAULT_FDR else red)
            elif col_i == 7:
                conf = row['conf']
                if conf in ("Ambiguous", "Weak"):
                    item.setForeground(yellow)
                elif conf == "Below threshold":
                    item.setForeground(red)
                elif conf == "High":
                    item.setForeground(green)
            elif col_i == 9:
                if row['flags_text']:
                    item.setForeground(yellow)
                    item.setToolTip(
                        row['flags_text'].replace('; ', chr(10))
                        + chr(10) * 2
                        + "Look at this cluster before trusting its label: "
                        "these are the signs of doublets or damaged cells, "
                        "which the per-cell QC cannot see. Right-click to "
                        "mark it Unknown if the markers confirm it.")

        self.cluster_table.set_data(df, decorate=_decorate)

        if author_labels is not None and n_labelled:
            text = (f"Authors' labels ({author_col}): {100 * n_agree / n_labelled:.1f}% "
                    f"of {n_labelled:,} labelled cells carry the label that is "
                    f"dominant in their cluster")
            if n_unlabelled:
                text += f"; {n_unlabelled:,} cells have no author label"
            self.author_agreement_label.setText(text + ".")
            self.author_agreement_label.setVisible(True)
        else:
            self.author_agreement_label.setVisible(False)

        has_annotations = (ct_col is not None
                           and 'cell_type_auto' in self.adata.obs.columns)
        self.save_annot_btn.setEnabled(has_annotations)
        self.undo_annot_btn.setEnabled(has_annotations)

        self._table_updating = False

    def _on_cluster_table_double_click(self, row, col):
        """Handle double-click on cluster table — open dropdown for cell type column."""
        if col != 3:  # "Assigned Type" column
            return

        available_types = ["Unknown"]
        for i in range(self.marker_list.count()):
            available_types.append(self.marker_list.item(i).text())

        if self.adata is not None and 'cell_type' in self.adata.obs.columns:
            for ct in self.adata.obs['cell_type'].unique():
                if ct not in available_types:
                    available_types.append(str(ct))

        combo = NoScrollComboBox()
        combo.addItems(sorted(available_types))

        current_text = self.cluster_table.item(row, col).text() if self.cluster_table.item(row, col) else ""
        idx = combo.findText(current_text)
        if idx >= 0:
            combo.setCurrentIndex(idx)

        self.cluster_table.setCellWidget(row, col, combo)
        combo.showPopup()

        def on_selection_changed():
            new_type = combo.currentText()
            self.cluster_table.removeCellWidget(row, col)
            self.cluster_table.setItem(row, col, QTableWidgetItem(new_type))
            self._apply_cluster_reannotation(row, new_type)

        combo.activated.connect(on_selection_changed)

    def _on_cluster_table_cell_changed(self, row, col):
        """Handle direct cell edits in the cluster table."""
        if getattr(self, '_table_updating', False):
            return
        if col != 3:  # "Assigned Type" column
            return
        item = self.cluster_table.item(row, col)
        if item:
            self._apply_cluster_reannotation(row, item.text())

    def _on_cluster_table_context_menu(self, pos):
        """Show right-click context menu on cluster table."""
        row = self.cluster_table.rowAt(pos.y())
        menu = QMenu(self)

        has_data = self.adata is not None and 'cell_type' in self.adata.obs.columns

        if row >= 0 and has_data:
            cluster_item = self.cluster_table.item(row, 0)
            cluster_id = cluster_item.text() if cluster_item else ""
            current_type = ""
            type_item = self.cluster_table.item(row, 3)
            if type_item:
                current_type = type_item.text()

            # Rename — opens the same dropdown as double-click
            rename_action = menu.addAction(f"Rename cluster {cluster_id}...")
            rename_action.triggered.connect(
                lambda: self._on_cluster_table_double_click(row, 3)
            )

            # Merge with — submenu of other types in the table
            merge_menu = menu.addMenu("Merge with...")
            other_types = set()
            for r in range(self.cluster_table.rowCount()):
                item = self.cluster_table.item(r, 3)
                if item and item.text() and item.text() != current_type:
                    other_types.add(item.text())
            for t in sorted(other_types):
                action = merge_menu.addAction(t)
                action.triggered.connect(
                    lambda checked, new_type=t: self._apply_cluster_reannotation(row, new_type)
                )
            if not other_types:
                merge_menu.setEnabled(False)

            # Mark as Unknown
            unknown_action = menu.addAction("Mark as Unknown")
            unknown_action.triggered.connect(
                lambda: self._apply_cluster_reannotation(row, "Unknown")
            )
            if current_type == "Unknown":
                unknown_action.setEnabled(False)

        menu.exec(self.cluster_table.viewport().mapToGlobal(pos))

    def _apply_cluster_reannotation(self, row, new_type):
        """Apply a manual re-annotation for a specific cluster."""
        if self.adata is None:
            return

        cluster_item = self.cluster_table.item(row, 0)
        if not cluster_item:
            return
        cluster_id = cluster_item.text()

        if 'leiden' in self.adata.obs.columns:
            cluster_col = 'leiden'
        elif 'clusters' in self.adata.obs.columns:
            cluster_col = 'clusters'
        elif 'seurat_clusters' in self.adata.obs.columns:
            cluster_col = 'seurat_clusters'
        else:
            return

        mask = self.adata.obs[cluster_col] == cluster_id
        self.adata.obs.loc[mask, 'cell_type'] = new_type

        # Provenance: tag this cluster as a manual override so the
        # cluster table can show which calls came from auto vs. user.
        provenance = dict(self.adata.uns.get('cluster_type_source', {}))
        provenance[str(cluster_id)] = 'manual'
        self.adata.uns['cluster_type_source'] = provenance

        self.main_window.current_adata = self.adata
        self.main_window._adata_version += 1
        self._last_adata_version = self.main_window._adata_version
        self._on_status(f"Cluster {cluster_id} re-annotated as '{new_type}' ({mask.sum():,} cells)")

    # ------------------------------------------------------------------
    # Top differential genes panel
    # ------------------------------------------------------------------

    def _setup_top_genes_panel(self, parent_layout):
        """
        Build the top-N differential genes table that sits below the
        cluster table. Populated when the user clicks a cluster row.
        """
        from PyQt6.QtWidgets import QSpinBox

        controls_row = QHBoxLayout()
        controls_row.setSpacing(4)
        self._top_genes_status = HintLabel(
            "Click a cluster row to see its top differential genes.")
        self._top_genes_status.setWordWrap(True)
        controls_row.addWidget(self._top_genes_status, 1)

        controls_row.addWidget(QLabel("Top N"))
        self._top_genes_n_spin = QSpinBox()
        self._top_genes_n_spin.setRange(5, 100)
        self._top_genes_n_spin.setValue(20)
        self._top_genes_n_spin.setFixedWidth(60)
        self._top_genes_n_spin.valueChanged.connect(
            lambda _v: self._refresh_top_genes_for_current_cluster())
        controls_row.addWidget(self._top_genes_n_spin)
        parent_layout.addLayout(controls_row)

        self._top_genes_table = ResultsTable()
        self._top_genes_table.set_schema([
            Column('Gene',    'gene',     's'),
            Column('log2FC',  'log2fc',   '.2f'),
            Column('adj p',   'pval_adj', '.2e', na_text='—'),
            Column('pct in',  'pct_in',   '.0%', na_text='—'),
            Column('pct out', 'pct_out',  '.0%', na_text='—'),
            Column('score',   'score',    '.2f'),
        ])
        self._top_genes_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents)
        self._top_genes_table.horizontalHeader().setStretchLastSection(True)
        self._top_genes_table.setAlternatingRowColors(True)
        self._top_genes_table.setSelectionBehavior(
            QTableWidget.SelectionBehavior.SelectRows)
        parent_layout.addWidget(self._top_genes_table, 1)

        self._rank_genes_worker = None
        self._current_top_genes_cluster = None

    def _on_cluster_row_clicked(self, row, _col):
        """Cluster table click -> show top genes for that cluster."""
        if self.adata is None:
            return
        cluster_item = self.cluster_table.item(row, 0)
        if cluster_item is None:
            return
        cluster_id = cluster_item.text()
        self._current_top_genes_cluster = cluster_id

        # If rank_genes_groups already cached, render straight away.
        if 'rank_genes_groups' in self.adata.uns:
            self._populate_top_genes_for_cluster(cluster_id)
            return

        # Need to compute first. Skip if a worker is already running.
        if self._rank_genes_worker is not None and self._rank_genes_worker.isRunning():
            self._top_genes_status.setText(
                "Computing top genes... will populate when done.")
            return

        cluster_col = self._cluster_col()
        if cluster_col is None:
            self._top_genes_status.setText(
                "No cluster column found. Run the Cluster tab first.")
            return

        self._top_genes_status.setText(
            "Computing top genes per cluster (one-time, ~1 min on big data)...")
        self._top_genes_table.setRowCount(0)

        self._rank_genes_worker = RankGenesGroupsWorker(self.adata, cluster_col)
        run_worker(
            self._rank_genes_worker,
            on_finished=self._on_rank_genes_finished,
            on_failed=self._on_rank_genes_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_progress,
        )

    def _on_rank_genes_finished(self, payload):
        adata, message = payload
        self.adata = adata
        self._on_status(message)
        # Render whichever cluster the user clicked (or had selected).
        if self._current_top_genes_cluster is not None:
            self._populate_top_genes_for_cluster(self._current_top_genes_cluster)
        else:
            self._top_genes_status.setText(
                "Click a cluster row to see its top differential genes.")

    def _on_rank_genes_failed(self, message: str):
        self._top_genes_status.setText(f"Top-genes computation failed: {message}")
        dialogs.warning(self, "Top genes failed", message)

    def _refresh_top_genes_for_current_cluster(self):
        if self._current_top_genes_cluster is None:
            return
        if self.adata is None or 'rank_genes_groups' not in self.adata.uns:
            return
        self._populate_top_genes_for_cluster(self._current_top_genes_cluster)

    def _populate_top_genes_for_cluster(self, cluster_id):
        """
        Pull the top-N DE genes for 'cluster_id' from cache and
        render them in the bottom table.
        """
        from kosmic.scrna.annotate.marker_genes import get_top_marker_genes

        n = self._top_genes_n_spin.value()
        try:
            df = get_top_marker_genes(self.adata, cluster_id, n=n)
        except (KeyError, ValueError) as exc:
            self._top_genes_status.setText(
                f"No top genes for cluster {cluster_id}: {exc}")
            self._top_genes_table.setRowCount(0)
            return

        self._top_genes_status.setText(
            f"Top {len(df)} differential genes for cluster {cluster_id} "
            "(vs all other clusters, Wilcoxon)")

        self._top_genes_table.set_data(df)

    def _cluster_col(self) -> Optional[str]:
        """Return the cluster column to use, preferring 'leiden'."""
        if self.adata is None:
            return None
        for cand in ('leiden', 'clusters', 'seurat_clusters', 'louvain'):
            if cand in self.adata.obs.columns:
                return cand
        return None

    def _save_annotations(self):
        """Save current annotations to h5ad file."""
        if self.adata is None:
            return

        if self.h5ad_path is None:
            if hasattr(self.main_window, 'inspect_tab') and self.main_window.inspect_tab.h5ad_path:
                self.h5ad_path = self.main_window.inspect_tab.h5ad_path
            elif self.project_dir:
                self.h5ad_path = processed_h5ad_path(self.project_dir, "annotated.h5ad")
            else:
                dialogs.warning(self, "No Path", "No save path available.")
                return

        try:
            Path(str(self.h5ad_path)).parent.mkdir(parents=True, exist_ok=True)
            self.adata.write_h5ad(str(self.h5ad_path))
            self._on_status(f"Annotations saved to {self.h5ad_path}")
            dialogs.info(self, "Saved", f"Annotations saved to:\n{self.h5ad_path}")
        except Exception as e:
            dialogs.warning(self, "Save Failed", f"Could not save: {str(e)}")

    def _undo_annotation_changes(self):
        """Revert cell_type to the last automated annotation (cell_type_auto)."""
        if self.adata is None:
            return
        if 'cell_type_auto' not in self.adata.obs.columns:
            dialogs.warning(self, "Cannot Undo", "No automated annotation to revert to.")
            return

        self.adata.obs['cell_type'] = self.adata.obs['cell_type_auto'].copy()
        self.main_window.current_adata = self.adata
        self.main_window._adata_version += 1
        self._last_adata_version = self.main_window._adata_version
        self._populate_cluster_table()
        self._on_status("Reverted to automated annotation")


    # ------------------------------------------------------------------
    # Annotation methods
    # ------------------------------------------------------------------

    def _on_method_changed(self, index):
        """Switch between annotation method panels."""
        self.annotation_stack.setCurrentIndex(index)

    def _resize_stack(self, index):
        """Resize stacked widget to fit only the current page."""
        for i in range(self.annotation_stack.count()):
            w = self.annotation_stack.widget(i)
            if i == index:
                w.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)
            else:
                w.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Ignored)
        self.annotation_stack.adjustSize()

    def _resolve_h5ad_path(self):
        """Get save path — uses central path from main_window, overwriting same file."""
        if hasattr(self.main_window, 'current_h5ad_path') and self.main_window.current_h5ad_path:
            self.h5ad_path = self.main_window.current_h5ad_path
            return str(self.h5ad_path)
        if self.h5ad_path is not None:
            return str(self.h5ad_path)
        if hasattr(self.main_window, 'inspect_tab') and self.main_window.inspect_tab.h5ad_path:
            self.h5ad_path = self.main_window.inspect_tab.h5ad_path
            return str(self.h5ad_path)
        dialogs.warning(self, "No Path", "No save path available. Run clustering first or load data.")
        return None

    def _run_annotation(self):
        """Run marker-based cell type annotation using PanglaoDB markers."""
        if self.adata is None:
            dialogs.warning(self, "No Data", "No data loaded. Run clustering first or load data in Inspect tab.")
            return

        selected_items = self.marker_list.selectedItems()
        if not selected_items:
            dialogs.warning(self, "No Cell Types", "Select cell types to annotate")
            return

        selected_cell_types = [item.text() for item in selected_items]
        organ = self.tissue_combo.currentText()

        try:
            from kosmic.reference.markers.panglaodb import get_markers
            species = self.species_combo.currentText().lower()
            marker_dict = get_markers(
                cell_types=selected_cell_types,
                organ=organ if organ != "All" else None,
                species=species,
                min_specificity=self.specificity_spin.value(),
                canonical_only=self.canonical_check.isChecked()
            )
        except Exception as e:
            self.log_message.emit(f"PanglaoDB error: {e}, falling back to defaults")
            marker_dict = {}
            for ct in selected_cell_types:
                if ct in DEFAULT_MARKERS:
                    marker_dict[ct] = DEFAULT_MARKERS[ct]

        if not marker_dict:
            dialogs.warning(self, "No Markers", "No markers found for selected cell types.")
            return

        output_path = self._resolve_h5ad_path()
        if output_path is None:
            return

        self.run_annot_btn.setEnabled(False)
        if self.progress_bar:
            self.progress_bar.setValue(0)

        annotation_params = {
            'score_threshold': self.score_threshold_spin.value(),
            'confidence_margin': self.confidence_margin_spin.value(),
            'run_ora': self.pdb_ora_check.isChecked(),
            'min_markers': self.pdb_min_markers_spin.value(),
        }

        self.marker_worker = MarkerScoringWorker(self.adata.copy(), marker_dict, output_path, annotation_params)
        run_worker(
            self.marker_worker,
            on_finished=self._on_annotation_finished,
            on_failed=self._on_annotation_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_progress,
        )

    def _on_annotation_finished(self, payload):
        adata, message = payload
        self.run_annot_btn.setEnabled(True)
        self.run_cm2_btn.setEnabled(True)
        self.run_ref_btn.setEnabled(True)

        self.adata = adata
        # Update central reference without full broadcast to all tabs
        self.main_window.current_adata = adata
        self.main_window._adata_version += 1
        self._last_adata_version = self.main_window._adata_version
        self._update_status()
        self.annotation_complete.emit(message)
        dialogs.info(self, "Annotation Complete", message)

    def _on_annotation_failed(self, message: str):
        self.run_annot_btn.setEnabled(True)
        self.run_cm2_btn.setEnabled(True)
        self.run_ref_btn.setEnabled(True)
        dialogs.warning(self, "Annotation Failed", message)

    def _run_celltypist(self):
        """Run CellTypist cell type annotation."""
        if self.adata is None:
            dialogs.warning(self, "No Data", "No data loaded. Run clustering first or load data in Inspect tab.")
            return

        model_name = self.model_combo.currentData()
        if not model_name:
            dialogs.warning(self, "No Model", "Select a CellTypist model (not a category header).")
            return

        output_path = self._resolve_h5ad_path()
        if output_path is None:
            return

        majority_voting = self.majority_voting_check.isChecked()

        self.run_celltypist_btn.setEnabled(False)
        if self.progress_bar:
            self.progress_bar.setValue(0)

        broad_labels = self.broad_labels_check.isChecked()

        self.celltypist_worker = CellTypistWorker(
            self.adata.copy(), model_name, majority_voting, output_path,
            broad_labels=broad_labels,
        )
        run_worker(
            self.celltypist_worker,
            on_finished=self._on_celltypist_finished,
            on_failed=self._on_celltypist_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_progress,
        )

    def _on_celltypist_finished(self, payload):
        """Handle CellTypist annotation completion."""
        adata, message = payload
        self.run_celltypist_btn.setEnabled(True)

        self.adata = adata
        self.main_window.current_adata = adata
        self.main_window._adata_version += 1
        self._last_adata_version = self.main_window._adata_version
        self._update_status()
        self.annotation_complete.emit(message)
        dialogs.info(self, "Annotation Complete", message)

    def _on_celltypist_failed(self, message: str):
        self.run_celltypist_btn.setEnabled(True)
        dialogs.warning(self, "Annotation Failed", message)

    def _run_cellmarker2_annotation(self):
        """Run cell type annotation using CellMarker 2.0 markers."""
        if self.adata is None:
            dialogs.warning(self, "No Data", "No data loaded. Run clustering first or load data in Inspect tab.")
            return

        selected_items = self.cm2_cell_list.selectedItems()
        if not selected_items:
            dialogs.warning(self, "No Cell Types", "Select cell types to annotate")
            return

        selected_cell_types = [item.text() for item in selected_items]
        species = self.cm2_species_combo.currentText()
        tissue = self.cm2_tissue_combo.currentText()
        cancer_type = self.cm2_cancer_combo.currentText()

        try:
            from kosmic.reference.markers.cellmarker2 import get_markers
            marker_dict = get_markers(
                cell_types=selected_cell_types,
                tissue=tissue if tissue != "All" else None,
                species=species,
                cancer_type=cancer_type
            )
        except Exception as e:
            dialogs.warning(self, "CellMarker 2.0 Error", f"Could not load markers: {str(e)}")
            return

        if not marker_dict:
            dialogs.warning(self, "No Markers", "No markers found for selected cell types.")
            return

        output_path = self._resolve_h5ad_path()
        if output_path is None:
            return

        self.run_cm2_btn.setEnabled(False)
        if self.progress_bar:
            self.progress_bar.setValue(0)

        annotation_params = {
            'score_threshold': self.cm2_score_threshold_spin.value(),
            'confidence_margin': self.cm2_confidence_margin_spin.value(),
            'skip_gene_conversion': True,
            'force_assignment': True,
            'run_ora': self.cm2_ora_check.isChecked(),
            'min_markers': self.cm2_min_markers_spin.value(),
        }

        self.marker_worker = MarkerScoringWorker(self.adata.copy(), marker_dict, output_path, annotation_params)
        run_worker(
            self.marker_worker,
            on_finished=self._on_annotation_finished,
            on_failed=self._on_annotation_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_progress,
        )

    # ------------------------------------------------------------------
    # Reference Atlas methods
    # ------------------------------------------------------------------

    def _populate_ref_combo(self):
        """Populate the reference atlas combo box."""
        from kosmic.scrna.annotate.reference import list_available_references
        refs = list_available_references()

        self.ref_combo.clear()
        # Group by species
        for species_label in ['human', 'mouse']:
            species_refs = {k: v for k, v in refs.items() if v['species'] == species_label}
            if not species_refs:
                continue
            header = f"\u2500\u2500 {species_label.title()} \u2500\u2500"
            self.ref_combo.addItem(header, None)
            idx = self.ref_combo.count() - 1
            self.ref_combo.model().item(idx).setEnabled(False)

            for name, info in sorted(species_refs.items()):
                suffix = ""
                if not info['installed']:
                    suffix = " (not installed)"
                self.ref_combo.addItem(f"  {name}{suffix}", name)

        # Custom / load file option
        self.ref_combo.addItem("\u2500\u2500 Custom \u2500\u2500", None)
        idx = self.ref_combo.count() - 1
        self.ref_combo.model().item(idx).setEnabled(False)
        self.ref_combo.addItem("  Load custom reference file...", "__custom__")

        # Select first valid item
        for i in range(self.ref_combo.count()):
            if self.ref_combo.itemData(i) is not None:
                self.ref_combo.setCurrentIndex(i)
                break

    def _on_ref_changed(self):
        """Update info label when reference selection changes."""
        from kosmic.scrna.annotate.reference import list_available_references
        ref_name = self.ref_combo.currentData()
        refs = list_available_references()
        if ref_name and ref_name != '__custom__' and ref_name in refs:
            info = refs[ref_name]
            desc = info.get('description', '')
            # A build record stored as JSON is shown as its headline counts only.
            if desc.startswith('{'):
                try:
                    import json as _json
                    rec = _json.loads(desc)
                    if 'n_kept' in rec and 'n_total' in rec:
                        desc = (f"{rec['n_kept']:,} nuclei kept of "
                                f"{rec['n_total']:,}; built {rec.get('built', '')[:10]}")
                except (ValueError, TypeError):
                    pass
            self.ref_info_label.setText(f"{desc}\nSource: {info.get('source', '')}")
        else:
            self.ref_info_label.setText("")

    def _run_reference_annotation(self):
        """Run reference-based cell type annotation."""
        if self.adata is None:
            dialogs.warning(self, "No Data", "No data loaded. Run clustering first.")
            return

        ref_name = self.ref_combo.currentData()
        if not ref_name:
            dialogs.warning(self, "No Reference", "Select a reference atlas (not a category header).")
            return

        if ref_name == '__custom__':
            from PyQt6.QtWidgets import QFileDialog
            path, _ = QFileDialog.getOpenFileName(
                self, "Load Reference File", "",
                "Reference Files (*.json *.json.gz);;All Files (*)"
            )
            if not path:
                return
            ref_name = path

        # Check if reference is installed
        from kosmic.scrna.annotate.reference import BUILTIN_REFERENCES, REFERENCES_DIR
        if ref_name in BUILTIN_REFERENCES:
            ref_path = REFERENCES_DIR / BUILTIN_REFERENCES[ref_name]['file']
            if not ref_path.exists():
                dialogs.warning(
                    self, "Reference Not Installed",
                    f"'{ref_name}' is not yet installed.\n\n"
                    f"To build a reference, use 'Build from h5ad...' with an annotated atlas,\n"
                    f"or place a reference file at:\n{ref_path}"
                )
                return

        output_path = self._resolve_h5ad_path()
        if output_path is None:
            return

        self.run_ref_btn.setEnabled(False)
        if self.progress_bar:
            self.progress_bar.setValue(0)

        self.ref_worker = ReferenceAnnotationWorker(
            self.adata.copy(),
            ref_name,
            output_path,
            min_correlation=self.ref_min_corr_spin.value(),
        )
        run_worker(
            self.ref_worker,
            on_finished=self._on_ref_annotation_finished,
            on_failed=self._on_ref_annotation_failed,
            on_progress=self._on_status,
            on_progress_pct=self._on_progress,
        )

    def _on_ref_annotation_finished(self, payload):
        """Handle reference annotation completion."""
        adata, message = payload
        self.run_ref_btn.setEnabled(True)

        self.adata = adata
        self.main_window.current_adata = adata
        self.main_window._adata_version += 1
        self._last_adata_version = self.main_window._adata_version
        self._update_status()
        self.annotation_complete.emit(message)
        dialogs.info(self, "Annotation Complete", message)

    def _on_ref_annotation_failed(self, message: str):
        self.run_ref_btn.setEnabled(True)
        dialogs.warning(self, "Annotation Failed", message)

    def _build_custom_reference(self):
        """Build a lightweight reference from an annotated h5ad file."""
        from PyQt6.QtWidgets import QFileDialog, QInputDialog

        h5ad_path, _ = QFileDialog.getOpenFileName(
            self, "Select Annotated Atlas h5ad", "",
            "H5AD Files (*.h5ad);;All Files (*)"
        )
        if not h5ad_path:
            return

        # Ask for cell type column name
        col_name, ok = QInputDialog.getText(
            self, "Cell Type Column",
            "Column name in obs containing cell type labels:",
            text="cell_type"
        )
        if not ok or not col_name:
            return

        # Ask for reference name
        from pathlib import Path as P
        default_name = P(h5ad_path).stem
        ref_name, ok = QInputDialog.getText(
            self, "Reference Name",
            "Name for this reference:",
            text=default_name
        )
        if not ok or not ref_name:
            return

        # Ask for species
        species, ok = QInputDialog.getItem(
            self, "Species", "Species:", ["human", "mouse"], 0, False
        )
        if not ok:
            return

        from kosmic.scrna.annotate.reference import REFERENCES_DIR
        output_path = REFERENCES_DIR / f"{ref_name.lower().replace(' ', '_')}.json.gz"

        self._on_status(f"Building reference from {P(h5ad_path).name}...")
        self.build_ref_btn.setEnabled(False)

        try:
            from kosmic.scrna.annotate.reference import build_reference_from_h5ad
            meta = build_reference_from_h5ad(
                h5ad_path, str(output_path),
                cell_type_col=col_name,
                name=ref_name,
                species=species,
            )
            self._on_status(
                f"Reference built: {meta['n_cell_types']} types, "
                f"{meta['n_genes']} genes \u2192 {output_path.name}"
            )
            self._populate_ref_combo()
            dialogs.info(
                self, "Reference Built",
                f"Reference '{ref_name}' created with {meta['n_cell_types']} cell types "
                f"and {meta['n_genes']} genes.\n\nSaved to: {output_path}"
            )
        except Exception as e:
            dialogs.warning(self, "Build Failed", f"Could not build reference: {str(e)}")
        finally:
            self.build_ref_btn.setEnabled(True)

    # ------------------------------------------------------------------
    # Combo box population helpers
    # ------------------------------------------------------------------

    def _populate_model_combo(self):
        """Populate the CellTypist model combo with categorized entries."""
        for category, models in CELLTYPIST_MODELS.items():
            self.model_combo.addItem(f"— {category} —", None)
            idx = self.model_combo.count() - 1
            self.model_combo.model().item(idx).setEnabled(False)

            for model_file, description in models:
                # Short display: extract type count from description
                short = model_file.replace('.pkl', '').replace('_', ' ')
                # Truncate long names
                if len(short) > 25:
                    short = short[:23] + '…'
                self.model_combo.addItem(short, model_file)
                full = f"{model_file.replace('.pkl', '')} — {description}"
                self.model_combo.setItemData(
                    self.model_combo.count() - 1, full, Qt.ItemDataRole.ToolTipRole
                )

        # Default to Healthy_Adult_Heart -- KOSMIC's primary use case is
        # cardiac scRNA. Falls through to the first non-header model
        # if that one is ever removed from the catalog.
        default_model = "Healthy_Adult_Heart.pkl"
        for i in range(self.model_combo.count()):
            if self.model_combo.itemData(i) == default_model:
                self.model_combo.setCurrentIndex(i)
                break
        else:
            for i in range(self.model_combo.count()):
                if self.model_combo.itemData(i) is not None:
                    self.model_combo.setCurrentIndex(i)
                    break

    def _populate_tissue_combo(self):
        """Populate tissue/organ combo box from PanglaoDB."""
        try:
            from kosmic.reference.markers.panglaodb import get_organs
            organs = get_organs()
            for organ in organs:
                self.tissue_combo.addItem(organ)
            idx = self.tissue_combo.findText("Heart")
            if idx >= 0:
                self.tissue_combo.setCurrentIndex(idx)
        except Exception as e:
            self.log_message.emit(f"Warning: Could not load PanglaoDB organs: {e}")

    def _populate_cell_types(self, organ: str = None):
        """Populate cell type list from PanglaoDB."""
        self.marker_list.clear()
        try:
            from kosmic.reference.markers.panglaodb import get_cell_types
            if organ and organ != "All":
                cell_types = get_cell_types(organ)
            else:
                cell_types = get_cell_types()

            for ct in cell_types:
                item = QListWidgetItem(ct)
                self.marker_list.addItem(item)
                item.setSelected(True)

            self.annot_status.setText(f"{len(cell_types)} cell types available")
        except Exception as e:
            self.log_message.emit(f"Warning: Could not load PanglaoDB cell types: {e}")
            for cell_type in DEFAULT_MARKERS.keys():
                item = QListWidgetItem(cell_type)
                self.marker_list.addItem(item)
                item.setSelected(True)

    def _on_tissue_changed(self, tissue: str):
        self._populate_cell_types(tissue)

    def _select_all_cell_types(self):
        for i in range(self.marker_list.count()):
            self.marker_list.item(i).setSelected(True)

    def _select_no_cell_types(self):
        for i in range(self.marker_list.count()):
            self.marker_list.item(i).setSelected(False)

    def _on_pdb_preset_changed(self, index: int):
        if index <= 0:
            return  # placeholder row
        name = self.pdb_preset_combo.itemText(index)
        preset = self._focus_presets.get(name, {})
        cell_types = set(preset.get('panglaodb', []))
        if not cell_types:
            return
        self.tissue_combo.setCurrentText("All")
        self._select_no_cell_types()
        selected_count = 0
        for i in range(self.marker_list.count()):
            item = self.marker_list.item(i)
            if item.text() in cell_types:
                item.setSelected(True)
                selected_count += 1
        self.annot_status.setText(
            f"{name} preset: {selected_count} cell types selected")
        QSettings("KOSMIC", "KOSMIC").setValue(_FOCUS_PRESET_QSETTING_KEY, name)

    # --- CellMarker 2.0 helpers ---

    def _populate_cm2_tissue_combo(self):
        try:
            from kosmic.reference.markers.cellmarker2 import get_tissues
            for tissue in get_tissues():
                self.cm2_tissue_combo.addItem(tissue)
            idx = self.cm2_tissue_combo.findText("Heart")
            if idx >= 0:
                self.cm2_tissue_combo.setCurrentIndex(idx)
        except Exception as e:
            self.log_message.emit(f"Warning: Could not load CellMarker 2.0 tissues: {e}")

    def _populate_cm2_cell_types(self):
        self.cm2_cell_list.clear()
        try:
            from kosmic.reference.markers.cellmarker2 import get_cell_types
            species = self.cm2_species_combo.currentText()
            tissue = self.cm2_tissue_combo.currentText()
            cell_types = get_cell_types(
                tissue=tissue if tissue != "All" else None,
                species=species
            )
            for ct in cell_types:
                item = QListWidgetItem(ct)
                self.cm2_cell_list.addItem(item)
                item.setSelected(True)
            self.annot_status.setText(f"CellMarker 2.0: {len(cell_types)} cell types available")
        except Exception as e:
            self.log_message.emit(f"Warning: Could not load CellMarker 2.0 cell types: {e}")

    def _on_cm2_filter_changed(self, _text=None):
        self._populate_cm2_cell_types()

    def _cm2_select_all(self):
        for i in range(self.cm2_cell_list.count()):
            self.cm2_cell_list.item(i).setSelected(True)

    def _cm2_select_none(self):
        for i in range(self.cm2_cell_list.count()):
            self.cm2_cell_list.item(i).setSelected(False)

    def _on_cm2_preset_changed(self, index: int):
        if index <= 0:
            return  # placeholder row
        name = self.cm2_preset_combo.itemText(index)
        preset = self._focus_presets.get(name, {})
        species = self.cm2_species_combo.currentText()
        key = 'cellmarker2_human' if species == 'Human' else 'cellmarker2_mouse'
        cell_types = set(preset.get(key, []))
        if not cell_types:
            self.annot_status.setText(
                f"{name} preset has no entries for {species}.")
            return
        self.cm2_tissue_combo.setCurrentText("Heart")
        self._cm2_select_none()
        selected_count = 0
        for i in range(self.cm2_cell_list.count()):
            item = self.cm2_cell_list.item(i)
            if item.text() in cell_types:
                item.setSelected(True)
                selected_count += 1
        self.annot_status.setText(
            f"CellMarker 2.0 {name} ({species}): "
            f"{selected_count} cell types selected")
        QSettings("KOSMIC", "KOSMIC").setValue(_FOCUS_PRESET_QSETTING_KEY, name)

    def _restore_last_preset(self):
        # Restore the saved preset name into both comboboxes (if present
        # in the current loaded preset list). Does not auto-apply -- the
        # marker lists are typically empty at this point; the user
        # re-picks once they've populated them.
        last = QSettings("KOSMIC", "KOSMIC").value(_FOCUS_PRESET_QSETTING_KEY, '')
        if not last:
            return
        for combo in (self.pdb_preset_combo, self.cm2_preset_combo):
            idx = combo.findText(last)
            if idx > 0:
                combo.blockSignals(True)
                combo.setCurrentIndex(idx)
                combo.blockSignals(False)

    def refresh_theme(self):
        """Re-apply theme to status labels."""
        self._update_status()
