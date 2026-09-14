# Intake workers: download from GEO and convert raw formats to h5ad.
# Moved verbatim from the scRNA download tab (ADR-001 Phase E).
from pathlib import Path
from typing import TYPE_CHECKING
import urllib.request

from PyQt6.QtCore import pyqtSignal

from kosmic.gui.shared.widgets import BaseWorker

if TYPE_CHECKING:
    import anndata


class GEOQueryWorker(BaseWorker):
    """
    Worker to query GEO for available files.

    Emits 'finished_ok' with '(files, message)'.
    """

    def __init__(self, gse: str):
        super().__init__()
        self.gse = gse.upper()

    def _run(self):
        from kosmic.scrna.load.download import query_geo_files
        try:
            data_files = query_geo_files(self.gse)
        except urllib.error.HTTPError as e:
            raise RuntimeError(
                f"GEO query failed: {e.code} - Check GSE number") from e

        if not data_files:
            raise RuntimeError("No supplementary files found")
        return (data_files, f"Found {len(data_files)} files")


class DownloadWorker(BaseWorker):
    """
    Worker thread for downloading files.

    Emits 'finished_ok' with a status message. Additionally declares
    'bytes_progress(int, int)' -- (downloaded_bytes, total_bytes) -- for
    byte-accurate progress bars alongside the inherited text 'progress'.
    """

    bytes_progress = pyqtSignal(int, int)

    def __init__(self, url: str, output_path: str, expected_size: int = 0):
        super().__init__()
        self.url = url
        self.output_path = output_path
        self.expected_size = expected_size
        self._cancelled = False

    def _run(self):
        from kosmic.scrna.load.download import download_file

        self.progress.emit("Connecting...")

        def on_progress(downloaded, total):
            self.bytes_progress.emit(downloaded, total)
            if total > 0:
                pct = downloaded * 100 // total
                self.progress.emit(f"{downloaded/1e6:.1f}/{total/1e6:.1f} MB ({pct}%)")
            else:
                self.progress.emit(f"Downloaded {downloaded/1e6:.1f} MB...")

        download_file(
            self.url, self.output_path,
            progress_callback=on_progress,
            cancel_check=lambda: self._cancelled,
            expected_size=self.expected_size,
        )
        return f"Downloaded: {Path(self.output_path).name}"

    def cancel(self):
        self._cancelled = True


def _study_dataset_path(output_dir, suffix: str = ".h5ad") -> Path:
    """Name an imported dataset after its study, not its source file.

    ``output_dir`` is a study's ``processed_data/``, so its parent is the
    accession. Source names carry the depositor's conventions
    ("human_dcm_hcm_scportal_03.17.2022.h5ad"), and those propagate: the
    Meta workspace recovers each study's accession from its DE CSV
    filename, and per-cell-type subsets inherit the stem. Naming from the
    study folder keeps that chain readable.
    """
    out = Path(output_dir)
    return out / f"{out.parent.name}{suffix}"


class H5adImportWorker(BaseWorker):
    """
    Worker to import an h5ad from raw_data/ into processed_data/.

    No format conversion is needed, but the counts matrix still has to be
    chosen: depositors disagree about which slot holds raw integer counts.
    CellxGene normalises 'X' and stores counts in '.raw'; Broad SCP keeps
    CellBender counts in 'X' and CellRanger counts in a layer. Handing
    log-normalised values to the pseudobulk DE step fails silently, so the
    slot is selected explicitly and recorded in the log.

    'slot' of None selects the best candidate automatically.

    Emits 'finished_ok' with '(output_path, message)'.
    """

    def __init__(self, h5ad_path: str, output_dir: str, slot: str = None):
        super().__init__()
        self.h5ad_path = h5ad_path
        self.output_dir = output_dir
        self.slot = slot

    def _run(self):
        from kosmic.scrna.load.matrices import (
            RAW_COUNTS, describe_count_matrices_path, read_counts,
        )

        src = Path(self.h5ad_path)
        out = _study_dataset_path(self.output_dir)

        # Inspect the file without loading it: a large object can hold the
        # same expression twice (normalised X plus raw counts), and reading
        # both to choose one exhausts memory on the bigger datasets.
        self.progress.emit(f"Inspecting {src.name}...")
        self.progress_pct.emit(5)
        described = describe_count_matrices_path(src)
        for d in described:
            self.progress.emit(
                f"  {d['label']}: {d['verdict']} (max {d['max']:,.1f})")

        by_slot = {d['slot']: d for d in described}
        slot = self.slot
        if slot is None:
            if by_slot.get('X', {}).get('verdict') == RAW_COUNTS:
                slot = 'X'
            else:
                slot = next((d['slot'] for d in described
                             if d['verdict'] == RAW_COUNTS), None)
        if slot is None:
            raise RuntimeError(
                f"{src.name} has no matrix that looks like raw integer "
                "counts. Pseudobulk DE needs counts; check the source file.")

        if slot == "X":
            self.progress.emit("X holds raw counts — using it as-is")
        else:
            self.progress.emit(
                f"X is not raw counts — using "
                f"{by_slot.get(slot, {}).get('label', slot)} instead")

        size_mb = src.stat().st_size / (1024 * 1024)
        self.progress.emit(f"Reading {src.name} ({size_mb:,.0f} MB)...")
        self.progress_pct.emit(25)
        adata = read_counts(src, slot)

        self.progress.emit(f"Writing {out.name}...")
        self.progress_pct.emit(75)
        out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out, compression="gzip")

        self.progress_pct.emit(100)
        message = (
            f"Imported: {out.name} ({adata.n_obs:,} cells x "
            f"{adata.n_vars:,} genes, counts from {slot})")
        return (str(out), message)


class RDSConvertWorker(BaseWorker):
    """
    Worker to convert Seurat RDS/Robj to h5ad using R + Python.

    Emits 'finished_ok' with '(output_path, message)'.
    """

    def __init__(self, rds_path: str, output_dir: str):
        super().__init__()
        self.rds_path = rds_path
        self.output_dir = output_dir

    def _run(self):
        from kosmic.scrna.load.rds_converter import (
            find_r_executable, run_rds_to_intermediate
        )
        from kosmic.scrna.load.converters import assemble_h5ad_from_r_export

        rds_name = Path(self.rds_path).stem
        if rds_name.lower().endswith('.robj'):
            rds_name = rds_name[:-5]
        elif rds_name.lower().endswith('.rds'):
            rds_name = rds_name[:-4]
        h5ad_path = _study_dataset_path(self.output_dir)

        rds_size_mb = Path(self.rds_path).stat().st_size / (1024 * 1024)
        self.progress.emit(f"Step 1/2: Loading R object ({rds_size_mb:.0f} MB)...")
        self.progress_pct.emit(5)

        # Step 1: R export via the analysis package (no path interpolation)
        def on_progress(msg):
            self.progress.emit(f"Step 1/2: {msg}")

        intermediate_dir = run_rds_to_intermediate(
            self.rds_path, self.output_dir, progress_callback=on_progress
        )

        if intermediate_dir is None:
            r_exe = find_r_executable()
            if not r_exe:
                raise RuntimeError(
                    "R not found. To convert RDS/Robj files, please either:\n\n"
                    "1. Install R and the Seurat package, OR\n"
                    "2. Download the CSV count matrix file instead "
                    "(works with pure Python)")
            raise RuntimeError("R export failed")

        # Step 2: assemble h5ad from intermediate files
        self.progress.emit("Step 2/2: Creating h5ad from intermediate files...")
        self.progress_pct.emit(50)

        adata = assemble_h5ad_from_r_export(intermediate_dir, rds_name, h5ad_path)

        self.progress_pct.emit(100)
        message = (
            f"Converted: {h5ad_path.name} "
            f"({adata.n_obs:,} cells x {adata.n_vars:,} genes)")
        return (str(h5ad_path), message)


class ExtractWorker(BaseWorker):
    """
    Worker thread for extracting archives.

    Emits 'finished_ok' with '(extracted_files, message)'.
    """

    def __init__(self, archive_path: str, output_dir: str):
        super().__init__()
        self.archive_path = archive_path
        self.output_dir = output_dir

    def _run(self):
        from kosmic.scrna.load.download import extract_archive

        self.progress.emit("Extracting...")
        extracted = extract_archive(
            self.archive_path, self.output_dir,
            progress_callback=lambda msg: self.progress.emit(msg),
        )
        return (extracted, f"Extracted {len(extracted)} files")


class MTXConvertWorker(BaseWorker):
    """
    Worker to convert MTX + metadata files to h5ad (for SCP and other sources).

    Emits 'finished_ok' with '(output_path, message)'.
    """

    def __init__(self, mtx_path: str, metadata_path: str, output_path: str,
                 genes_path: str = None, barcodes_path: str = None,
                 clusters_path: str = None):
        super().__init__()
        self.mtx_path = Path(mtx_path)
        self.metadata_path = Path(metadata_path) if metadata_path else None
        self.output_path = Path(output_path)
        self.genes_path = Path(genes_path) if genes_path else None
        self.barcodes_path = Path(barcodes_path) if barcodes_path else None
        self.clusters_path = Path(clusters_path) if clusters_path else None

    def _run(self):
        from kosmic.scrna.load.converters import mtx_to_h5ad

        self.progress.emit("Converting MTX to h5ad...")
        self.progress_pct.emit(10)

        adata = mtx_to_h5ad(
            str(self.mtx_path),
            output_path=str(self.output_path),
            genes_path=str(self.genes_path) if self.genes_path else None,
            barcodes_path=str(self.barcodes_path) if self.barcodes_path else None,
            metadata_path=str(self.metadata_path) if self.metadata_path else None,
            clusters_path=str(self.clusters_path) if self.clusters_path else None,
        )

        self.progress_pct.emit(100)
        message = (
            f"Created {self.output_path.name}: "
            f"{adata.n_obs:,} cells x {adata.n_vars:,} genes")
        return (str(self.output_path), message)


class CSVConvertWorker(BaseWorker):
    """
    Worker to convert CSV count matrix to h5ad.

    Emits 'finished_ok' with '(output_path, message)'.
    """

    def __init__(self, csv_path: str, output_dir: str):
        super().__init__()
        self.csv_path = csv_path
        self.output_dir = output_dir

    def _run(self):
        from kosmic.scrna.load.converters import csv_to_h5ad

        csv_path = Path(self.csv_path)
        h5ad_path = _study_dataset_path(self.output_dir)

        self.progress.emit("Converting CSV to h5ad...")
        self.progress_pct.emit(10)

        adata = csv_to_h5ad(str(csv_path), str(h5ad_path))

        self.progress_pct.emit(100)
        message = (
            f"Created {h5ad_path.name}: "
            f"{adata.n_obs:,} cells x {adata.n_vars:,} genes")
        return (str(h5ad_path), message)


class ExtractAndConvert10xWorker(BaseWorker):
    """
    Worker thread for extracting sample tar.gz files and converting to h5ad.

    Emits 'finished_ok' with '(output_path, message)'.
    """

    def __init__(self, raw_data_dir: str, output_path: str):
        super().__init__()
        self.raw_data_dir = Path(raw_data_dir)
        self.output_path = Path(output_path)

    def _run(self):
        from kosmic.scrna.load.download import load_multi_sample_10x

        self.progress.emit("Scanning for sample tar.gz files...")
        self.progress_pct.emit(5)

        adata = load_multi_sample_10x(
            str(self.raw_data_dir), str(self.output_path),
            progress_callback=lambda msg: self.progress.emit(msg),
        )

        self.progress_pct.emit(100)
        n_samples = adata.obs['sample'].nunique() if 'sample' in adata.obs.columns else 1
        msg = (
            f"Extracted & loaded {n_samples} sample(s): "
            f"{adata.n_obs:,} cells x {adata.n_vars:,} genes")
        return (str(self.output_path), msg)


class BatchConvertWorker(BaseWorker):
    """
    Generic batch-convert worker. Calls
    'converter_fn(path, progress_emit) -> AnnData' once per item and
    emits 'finished_ok' with '(adatas_list, message)' where adatas_list
    is '[(sample_name, adata), ...]'. 'progress_emit' is forwarded to
    the UI log. Pass 'sort_items=False' to preserve insertion order
    (used for matrix.txt where filename ordering matters).
    """

    def __init__(self, items, converter_fn, label: str,
                 *, sort_items: bool = True):
        super().__init__()
        self.items = [str(p) for p in items]
        self._converter_fn = converter_fn
        self._label = label
        self._sort = sort_items

    def _run(self):
        items = sorted(self.items) if self._sort else list(self.items)
        total = len(items)
        all_adatas = []

        for i, path in enumerate(items):
            self.progress.emit(
                f"[{i + 1}/{total}] Processing: {Path(path).name}")
            self.progress_pct.emit(int((i / total) * 90))

            adata = self._converter_fn(
                path, progress_emit=self.progress.emit)
            if 'sample' in adata.obs.columns:
                sample_name = adata.obs['sample'].iloc[0]
            else:
                sample_name = Path(path).stem
            all_adatas.append((sample_name, adata))

        self.progress_pct.emit(95)
        return (all_adatas,
                f"Converted {len(all_adatas)} {self._label} file(s)")


def _convert_tsv(path: str, *, progress_emit) -> "anndata.AnnData":
    from kosmic.scrna.load.converters import tsv_to_h5ad
    return tsv_to_h5ad(path, progress_callback=progress_emit)


def _convert_h5(path: str, *, progress_emit) -> "anndata.AnnData":
    from kosmic.scrna.load.converters import h5_to_h5ad
    adata = h5_to_h5ad(path)
    sample_name = (adata.obs['sample'].iloc[0]
                   if 'sample' in adata.obs.columns else Path(path).stem)
    progress_emit(
        f"  Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes "
        f"(sample: {sample_name})")
    return adata


def _convert_matrix_txt(path: str, *, progress_emit) -> "anndata.AnnData":
    from kosmic.scrna.load.converters import convert_matrix_to_h5ad
    return convert_matrix_to_h5ad(
        path, output_file=None, progress_callback=progress_emit)


class TenXFolderConvertWorker(BaseWorker):
    """
    Worker to convert 10X folder data (MTX + features + barcodes) to h5ad.

    Emits 'finished_ok' with '(adatas_list, message)'.
    """

    def __init__(self, folders):
        super().__init__()
        self.folders = [str(f) for f in folders]

    def _run(self):
        from kosmic.scrna.load.converters import load_10x_folder, find_10x_files, load_prefixed_10x
        import scanpy as sc

        total = len(self.folders)
        all_adatas = []

        for i, folder_str in enumerate(self.folders):
            folder = Path(folder_str)
            self.progress.emit(f"[{i+1}/{total}] Processing: {folder.name}")
            self.progress_pct.emit(int((i / total) * 90))

            # Standard 10X files
            mtx_files = list(folder.glob("matrix.mtx.gz")) + list(folder.glob("matrix.mtx"))
            if mtx_files:
                adata = load_10x_folder(folder)
                adata.obs['sample'] = folder.name
                all_adatas.append((folder.name, adata))
                self.progress.emit(
                    f"  Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
            else:
                # Prefixed files
                samples = find_10x_files(folder)
                prefixed = [s for s in samples if s.get('prefix')]
                if prefixed:
                    sample_adatas = []
                    for s_info in prefixed:
                        self.progress.emit(f"  Loading sample: {s_info['sample']}")
                        adata = load_prefixed_10x(s_info)
                        adata.obs['sample'] = s_info['sample']
                        sample_adatas.append(adata)
                        self.progress.emit(
                            f"    {adata.n_obs:,} cells x {adata.n_vars:,} genes")

                    if len(sample_adatas) > 1:
                        combined = sc.concat(sample_adatas, label='sample',
                                             keys=[a.obs['sample'].iloc[0] for a in sample_adatas])
                        combined.obs_names_make_unique()
                        all_adatas.append((folder.name, combined))
                    elif sample_adatas:
                        all_adatas.append((sample_adatas[0].obs['sample'].iloc[0], sample_adatas[0]))

        self.progress_pct.emit(95)
        return (all_adatas, f"Converted {len(all_adatas)} 10X folder(s)")



# ---------------------------------------------------------------------------




class _H5adLoadWorker(BaseWorker):
    """
    Background worker to load an h5ad file without blocking the UI.

    Emits 'finished_ok' with '(adata, path, message)'.
    """

    def __init__(self, h5ad_path: str):
        super().__init__()
        self.h5ad_path = h5ad_path

    def _run(self):
        from kosmic.scrna.load.converters import load_h5ad_with_summary
        adata, summary = load_h5ad_with_summary(self.h5ad_path)
        msg = f"{adata.n_obs:,} cells x {adata.n_vars:,} genes"
        return (adata, self.h5ad_path, msg)
