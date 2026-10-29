# StudySource -- where Figure Export's study figures get their data (ADR-014).
#
# A figure draws from the workspace that holds its study when one does,
# and otherwise from the study's folder: the processed h5ad (read in parts
# by kosmic.scrna.load.h5ad_meta) and the DE outputs (kosmic.de.saved).
# The figure pages ask this object and do not know which served them.
#
# - scRNA: 'scrna_adata()' is the scRNA workspace's dataset when it holds
#   this study (it can carry unsaved state, such as an open subclustering),
#   else obs, var, uns and the embeddings read from the file. That read
#   takes 1-2 s on a large atlas, so it runs on a background thread; until
#   it finishes 'scrna_adata()' is None and 'scrna_status()' says so, and
#   'changed' is emitted when it lands. Expression is read on demand
#   ('scrna_expression', 'scrna_cells'), from a figure's render thread.
# - DE: 'de_state()' answers to the DE workspace's attribute names. Each
#   attribute is the DE workspace's value when it is on this study and has
#   one, else what the study's folder holds; a saved run chosen with
#   'set_de_run' serves the run-specific attributes from its table and
#   record instead. 'available(name)' tells whether an attribute has a
#   value without reading a table, so filling the figure list stays cheap.
#
# Everything read from a file is cached against the file's modification
# time and size, so a file saved again is read again. A read that meets a
# save of the file ('StudyFileChanged') is retried once.
from __future__ import annotations

import threading
from pathlib import Path
from typing import Any, Callable, Optional

from PyQt6.QtCore import QObject, pyqtSignal


def _same_folder(a, b) -> bool:
    if not a or not b:
        return False
    try:
        return Path(a).resolve() == Path(b).resolve()
    except OSError:
        return False


def _signature(path: Path) -> Optional[tuple]:
    try:
        st = path.stat()
    except OSError:
        return None
    return (str(path), st.st_mtime_ns, st.st_size)


def _empty(value) -> bool:
    if value is None:
        return True
    if hasattr(value, 'empty'):
        return bool(value.empty)
    try:
        return len(value) == 0
    except TypeError:
        return False


class StudySource(QObject):
    """The data Figure Export's study figures draw, for one study."""

    #: The study changed, or a background read finished: availability may
    #: have changed.
    changed = pyqtSignal()

    def __init__(self, scrna_ws: Callable[[], Any], de_ws: Callable[[], Any], parent=None):
        super().__init__(parent)
        self._scrna_ws, self._de_ws = scrna_ws, de_ws
        self._study: Optional[Path] = None
        self._lock = threading.Lock()
        self._meta = None                  # (signature, AnnData) of the file read
        self._meta_error: Optional[str] = None
        self._meta_thread: Optional[threading.Thread] = None
        self._reads: dict = {}             # cached expression reads
        self._saved: dict = {}             # cached saved DE outputs
        self._de_run = None                # a kosmic.de.saved.GeneDERun, or None

    # -- The study ---------------------------------------------------------

    def study(self) -> Optional[Path]:
        return self._study

    def set_study(self, study_dir) -> bool:
        """Show 'study_dir' (None for none). Returns True when it changed."""
        new = Path(study_dir) if study_dir else None
        if _same_folder(new, self._study) or (new is None and self._study is None):
            return False
        with self._lock:
            self._study = new
            self._meta, self._meta_error = None, None
            self._reads.clear()
            self._saved.clear()
            self._de_run = None
        self.changed.emit()
        return True

    def h5ad_path(self) -> Optional[Path]:
        """The study's registered h5ad, or None."""
        if self._study is None or not self._study.is_dir():
            return None
        from kosmic.scrna.load.h5ad_meta import registered_h5ad
        try:
            return registered_h5ad(self._study)
        except OSError:
            return None

    # -- scRNA ---------------------------------------------------------------

    def _held_scrna(self):
        ws = self._scrna_ws()
        if ws is None or getattr(ws, 'current_adata', None) is None:
            return None
        if not _same_folder(getattr(ws, 'current_project_dir', None), self._study):
            return None
        return ws.current_adata

    def holds_scrna(self) -> bool:
        """True when the scRNA workspace holds this study."""
        return self._held_scrna() is not None

    def scrna_adata(self):
        """The study's dataset for figures that need obs, var, uns and the
        embeddings: the scRNA workspace's when it holds the study, else the
        file's (None until the background read has finished)."""
        held = self._held_scrna()
        if held is not None:
            return held
        path = self.h5ad_path()
        if path is None:
            return None
        sig = _signature(path)
        with self._lock:
            if self._meta is not None and self._meta[0] == sig:
                return self._meta[1]
        self._start_meta_read(path, sig)
        return None

    def scrna_status(self) -> Optional[str]:
        """Why 'scrna_adata()' is None, for a figure's tooltip; None when it
        is not (or there is no study)."""
        if self._study is None or self._held_scrna() is not None:
            return None
        if not self._study.is_dir():
            return f"The study folder {self._study} cannot be found."
        path = self.h5ad_path()
        if path is None:
            return f"{self._study.name} has no processed dataset yet."
        with self._lock:
            if self._meta is not None and self._meta[0] == _signature(path):
                return None
            if self._meta_error:
                return self._meta_error
        return f"Reading {path.name}..."

    def _start_meta_read(self, path: Path, sig) -> None:
        with self._lock:
            if self._meta_thread is not None and self._meta_thread.is_alive():
                return
            if self._meta_error and self._meta_error_sig == sig:
                return
            study = self._study
            t = threading.Thread(target=self._read_meta, args=(study, path, sig),
                                 name='figure-study-read', daemon=True)
            self._meta_thread = t
        t.start()

    def _read_meta(self, study, path: Path, sig) -> None:
        from kosmic.scrna.load.h5ad_meta import StudyFileChanged, read_figure_meta
        error, adata = None, None
        for attempt in (1, 2):
            try:
                adata = read_figure_meta(path)
                break
            except StudyFileChanged:
                sig = _signature(path)
                if attempt == 2:
                    error = f"{path.name} is being saved; open the figure again."
            except Exception as exc:  # noqa: BLE001 -- shown as the figure's tooltip
                error = f"Could not read {path.name}: {exc}"
                break
        with self._lock:
            if not _same_folder(study, self._study):
                return                      # the study changed meanwhile
            if adata is not None:
                self._meta, self._meta_error = (sig, adata), None
            else:
                self._meta_error, self._meta_error_sig = error, sig
        try:
            self.changed.emit()
        except RuntimeError:
            pass                            # the workspace has been closed

    _meta_error_sig = None

    def wait(self, timeout: float = 30) -> None:
        """Wait for a background read (tests, and before quitting)."""
        t = self._meta_thread
        if t is not None:
            t.join(timeout)

    def _cached_read(self, key, path: Path, read: Callable):
        from kosmic.scrna.load.h5ad_meta import StudyFileChanged
        full = (key, _signature(path))
        with self._lock:
            if full in self._reads:
                return self._reads[full]
        try:
            value = read()
        except StudyFileChanged:
            value = read()                  # once more, on the file as now saved
            full = (key, _signature(path))
        with self._lock:
            if len(self._reads) >= 4:
                self._reads.pop(next(iter(self._reads)))
            self._reads[full] = value
        return value

    def scrna_expression(self, genes):
        """Log-normalised expression of 'genes' for every cell: the held
        dataset itself (the figure reads the genes from it, as before), or
        those genes read from the file. Call from a render thread."""
        held = self._held_scrna()
        if held is not None:
            return held
        path = self.h5ad_path()
        if path is None:
            return None
        from kosmic.scrna.load.h5ad_meta import read_genes
        genes = [str(g) for g in genes]
        return self._cached_read(('genes', tuple(genes)), path,
                                 lambda: read_genes(path, genes, source='log'))

    def scrna_cells(self, rows):
        """Every gene for the cells at positions 'rows' (sorted, unique),
        log-normalised. Call from a render thread."""
        import numpy as np
        rows = np.unique(np.asarray(rows, dtype=np.int64))
        held = self._held_scrna()
        if held is not None:
            return held[rows]
        path = self.h5ad_path()
        if path is None:
            return None
        from kosmic.scrna.load.h5ad_meta import read_cells
        key = ('cells', len(rows), hash(rows.tobytes()))
        return self._cached_read(key, path, lambda: read_cells(path, rows, source='log'))

    # -- DE ----------------------------------------------------------------------

    def _live_de(self):
        ws = self._de_ws()
        if ws is not None and _same_folder(getattr(ws, 'project_dir', None), self._study):
            return ws
        return None

    def de_runs(self) -> list:
        """The study's saved gene DE runs (kosmic.de.saved.GeneDERun)."""
        if self._study is None:
            return []
        from kosmic.de.saved import list_gene_de_runs
        return list_gene_de_runs(self._study)

    def set_de_run(self, run) -> None:
        """Serve the run-specific DE attributes from this saved run (None:
        the DE workspace's run when it has one, else the first saved)."""
        self._de_run = run

    def de_run(self):
        return self._de_run

    def de_state(self) -> 'DEState':
        return DEState(self)

    def de_counts(self, genes, counts_source: str = 'raw'):
        """The counts of 'genes' for every cell, for the Gene Heatmap: the DE
        workspace's dataset when it holds this study, else read from the
        file from the count source the run used ('raw' or a layer name).
        Call from a render thread."""
        live = self._live_de()
        if live is not None and getattr(live, 'current_adata', None) is not None:
            return live.current_adata
        path = self.h5ad_path()
        if path is None:
            return None
        from kosmic.scrna.load.h5ad_meta import read_genes
        genes = [str(g) for g in genes]
        layer = None if counts_source in (None, '', 'raw') else counts_source
        return self._cached_read(
            ('counts', tuple(genes), layer), path,
            lambda: read_genes(path, genes, source='counts', counts_layer=layer))

    def _saved_value(self, key, path: Optional[Path], read: Callable):
        """A saved output, cached against its file."""
        sig = _signature(path) if path is not None else None
        with self._lock:
            if key in self._saved and self._saved[key][0] == sig:
                return self._saved[key][1]
        value = read() if sig is not None else None
        with self._lock:
            self._saved[key] = (sig, value)
        return value


class DEState:
    """The DE workspace's attribute names, answered for one study (see the
    module header). Created per use by 'StudySource.de_state()'."""

    #: Attributes that belong to one gene DE run.
    RUN_ATTRS = ('de_results', 'pathway_gene_sets', 'geneset_label', 'sample_col',
                 'condition_col', 'disease_label', 'control_label', 'min_cells',
                 'pathway_coverage', 'counts_source')

    def __init__(self, source: StudySource):
        self._source = source
        self.project_dir = source.study()
        live = source._live_de()
        self._run = source.de_run()
        # The DE workspace's run is the one shown unless a saved run was
        # chosen, or it has none and the study has a saved one.
        self._live = live
        self._live_run = (self._run is None and live is not None
                          and getattr(live, 'de_results', None) is not None)
        if self._run is None and not self._live_run:
            runs = source.de_runs()
            self._run = runs[0] if runs else None

    # -- What the DE workspace would answer --------------------------------------

    def _live_value(self, name):
        if self._live is None:
            return None
        if name in self.RUN_ATTRS and not self._live_run:
            return None
        if name == 'counts_source':
            from kosmic.scrna.counts import loaded_count_source
            adata = getattr(self._live, 'current_adata', None)
            return loaded_count_source(adata) if adata is not None else None
        value = getattr(self._live, name, None)
        return None if _empty(value) else value

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        value = self._live_value(name)
        if value is not None:
            return value
        return self._saved(name)

    @property
    def h5ad_path(self):
        live = self._live_value('h5ad_path')
        return live if live is not None else self._source.h5ad_path()

    @property
    def current_adata(self):
        """Only the DE workspace's dataset; a figure that needs counts from
        the file asks 'StudySource.de_counts'."""
        return self._live_value('current_adata')

    @property
    def run_label(self) -> Optional[str]:
        if self._live_run:
            return "Run in the DE workspace"
        return self._run.label if self._run is not None else None

    # -- What the study's folder holds -------------------------------------------

    def _record(self) -> dict:
        if self._run is None:
            return {}
        from kosmic.de.saved import read_run_record
        rec = self._source._saved_value(('record', str(self._run.path)),
                                        self._run.record_path,
                                        lambda: read_run_record(self._run.path))
        return rec or {}

    def _semantics(self) -> dict:
        from kosmic.manifest import read_manifest
        study = self.project_dir
        if study is None:
            return {}
        return ((read_manifest(study) or {}).get('semantics') or {})

    def _label(self, role: str) -> str:
        rec = self._record().get(f'{role}_label')
        if rec:
            return rec
        roles = self._semantics().get('role_map') or {}
        return next((cond for cond, r in roles.items() if r == role), '')

    def _saved(self, name):
        from kosmic.de import saved
        from kosmic.paths import de_pathway_scoring_dir
        src, study = self._source, self.project_dir
        if study is None:
            return None
        rec = self._record()
        if name == 'de_results':
            if self._run is None:
                return None
            run = self._run
            return src._saved_value(('de', str(run.path)), run.path,
                                    lambda: saved.read_gene_de(run))
        if name == 'pathway_gene_sets':
            return rec.get('pathway_gene_sets') or {}
        if name == 'geneset_label':
            return rec.get('geneset_label', '')
        if name == 'sample_col':
            return rec.get('sample_column') or self._semantics().get('sample_column', '')
        if name == 'condition_col':
            return rec.get('condition_column') or self._semantics().get('condition_column', '')
        if name == 'disease_label':
            return self._label('disease')
        if name == 'control_label':
            return self._label('control')
        if name == 'min_cells':
            from kosmic import DE_MIN_CELLS
            return rec.get('min_cells') or DE_MIN_CELLS
        if name == 'counts_source':
            return rec.get('counts_source', 'raw')
        if name == 'pathway_coverage':
            return self._coverage(rec.get('pathway_gene_sets') or {})
        scoring = de_pathway_scoring_dir(study)
        if name == 'pathway_de_results':
            return src._saved_value(('pathway_de', ''), scoring / 'pathway_de_results.csv',
                                    lambda: saved.read_pathway_de(study))
        if name == 'pathway_score_data':
            from kosmic.de.pathway_scoring import PER_DONOR_SCORES_FILE
            return src._saved_value(('scores', ''), scoring / PER_DONOR_SCORES_FILE,
                                    lambda: saved.read_per_donor_scores(study))
        if name == 'fgsea_results':
            from kosmic.de.fgsea import FGSEA_RESULTS_FILE
            return src._saved_value(('fgsea',), scoring / FGSEA_RESULTS_FILE,
                                    lambda: saved.read_fgsea(study))
        if name == 'gsea_library_results':
            return src._saved_value(('gsea',), saved.gsea_library_path(study),
                                    lambda: saved.read_gsea_library(study))
        if name in ('gsea_library_details', 'gsea_library_name'):
            gsea = saved.gsea_library_path(study)
            details, library = src._saved_value(
                ('gsea_details',), saved.gsea_details_path(gsea),
                lambda: saved.read_gsea_details(gsea)) or (None, '')
            return details if name == 'gsea_library_details' else (library or None)
        if name == 'enrichment_results':
            return src._saved_value(('ora',), saved.enrichment_path(scoring),
                                    lambda: saved.read_enrichment(scoring))
        return None

    def _coverage(self, gene_sets: dict):
        """Coverage of the run's gene sets in the study's genes, computed as
        the DE page computes it when it loads saved results."""
        if not gene_sets:
            return None
        path = self._source.h5ad_path()
        if path is None:
            return None

        def _compute():
            from kosmic.de.de_analysis import prepare_gene_coverage
            from kosmic.scrna.inspect.detection import detect_species
            from kosmic.scrna.load.h5ad_meta import read_var
            names = [str(g) for g in read_var(path).index]
            return prepare_gene_coverage(gene_sets, set(names), detect_species(names))[1]

        key = ('coverage', tuple(sorted(gene_sets)), str(self._run.path if self._run else ''))
        return self._source._saved_value(key, path, _compute)

    # -- Cheap availability ------------------------------------------------------

    def available(self, name: str) -> bool:
        """Whether 'name' has a value, without reading a table: the DE
        workspace's value, else whether the saved file exists."""
        if self._live_value(name) is not None:
            return True
        study = self.project_dir
        if study is None:
            return False
        from kosmic.de import saved
        from kosmic.paths import de_pathway_scoring_dir
        scoring = de_pathway_scoring_dir(study)
        if name == 'de_results':
            return self._run is not None and self._run.path.is_file()
        if name in ('pathway_gene_sets', 'geneset_label'):
            return bool(self._record().get('pathway_gene_sets'))
        if name == 'pathway_coverage':
            return (bool(self._record().get('pathway_gene_sets'))
                    and self._source.h5ad_path() is not None)
        if name in ('disease_label', 'control_label'):
            return bool(self._label(name.split('_')[0]))
        if name in ('sample_col', 'condition_col'):
            return bool(getattr(self, name))
        if name == 'pathway_de_results':
            return (scoring / 'pathway_de_results.csv').is_file()
        if name == 'pathway_score_data':
            from kosmic.de.pathway_scoring import PER_DONOR_SCORES_FILE
            return (scoring / PER_DONOR_SCORES_FILE).is_file()
        if name == 'fgsea_results':
            from kosmic.de.fgsea import FGSEA_RESULTS_FILE
            return (scoring / FGSEA_RESULTS_FILE).is_file()
        if name == 'gsea_library_results':
            return saved.gsea_library_path(study).is_file()
        if name == 'gsea_library_details':
            return saved.gsea_details_path(saved.gsea_library_path(study)).is_file()
        if name == 'enrichment_results':
            return saved.enrichment_path(scoring).is_file()
        return False
