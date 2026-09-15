# RDS/Seurat File Converter
# Convert Seurat RDS/Robj files to h5ad using three strategies:
# 1. rds2py (pure Python, handles Seurat objects directly)
# 2. pyreadr (pure Python, for simpler RDS files)
# 3. R subprocess (full Seurat support, requires R + Seurat installed)
#
# **Issue #8 fix:** R scripts receive paths via commandArgs(trailingOnly=TRUE)
# instead of f-string interpolation, eliminating script injection risk.

import os
import shutil
import subprocess
from pathlib import Path
from typing import Callable, Optional, Tuple

import numpy as np
import pandas as pd


# R executable finder

_UNIX_R_PATHS = (
    "/usr/local/bin/Rscript",
    "/opt/homebrew/bin/Rscript",
)
_WINDOWS_R_ROOT = Path(r"C:\Program Files\R")


def find_r_executable() -> Optional[str]:
    """
    Find the Rscript executable on the system.

    Checks PATH first, then common installation locations.

    Returns
    -------
    str or None
        Path to Rscript, or None if not found.
    """
    r_exe = shutil.which("Rscript")
    if r_exe:
        return r_exe

    for rpath in _UNIX_R_PATHS:
        if Path(rpath).exists():
            return rpath

    # Windows: pick the highest-versioned R install under Program Files.
    if _WINDOWS_R_ROOT.exists():
        candidates = sorted(_WINDOWS_R_ROOT.glob("R-*/bin/Rscript.exe"),
                            reverse=True)
        if candidates:
            return str(candidates[0])

    return None


# R script template (Issue #8 fix: uses commandArgs, no path interpolation)

# This script receives paths via command-line arguments, NOT embedded in the
# script text. This prevents injection attacks from malicious filenames.
_RDS_EXPORT_SCRIPT = r"""
# Seurat RDS/Robj to intermediate files exporter
# Arguments: <rds_path> <output_dir> <prefix> <is_robj>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript export.R <rds_path> <output_dir> <prefix> <is_robj>")
}

file_path <- args[1]
output_dir <- args[2]
output_prefix <- file.path(output_dir, args[3])
is_robj <- tolower(args[4]) == "true"

library(Seurat)
library(Matrix)

# Helper for garbage collection
gc_full <- function() {
  for(i in 1:3) gc(verbose = FALSE, full = TRUE)
}

message("PROGRESS:5:Loading R object file...")

if (is_robj) {
  # Robj files use load()
  message("Detected Robj format, using load()...")
  env <- new.env()
  is_gzipped <- endsWith(tolower(file_path), ".gz")
  if (is_gzipped) {
    message("Decompressing gzipped file...")
    temp_file <- tempfile(fileext = ".Robj")
    con_in <- gzfile(file_path, "rb")
    con_out <- file(temp_file, "wb")
    while (length(chunk <- readBin(con_in, "raw", 1024*1024)) > 0) {
      writeBin(chunk, con_out)
    }
    close(con_in)
    close(con_out)
    message("Loading decompressed file...")
    load(temp_file, envir = env)
    unlink(temp_file)
  } else {
    load(file_path, envir = env)
  }
  obj_names <- ls(env)
  message("Found objects: ", paste(obj_names, collapse=", "))
  obj <- NULL
  for (name in obj_names) {
    candidate <- get(name, envir = env)
    if (inherits(candidate, "Seurat")) {
      obj <- candidate
      message("Using Seurat object: ", name)
      break
    }
  }
  if (is.null(obj)) {
    stop("No Seurat object found in .Robj file. Found: ", paste(obj_names, collapse=", "))
  }
} else {
  # RDS files
  is_gzipped <- endsWith(tolower(file_path), ".gz")
  if (is_gzipped) {
    message("Decompressing gzipped RDS file...")
    con <- gzfile(file_path, "rb")
    raw_data <- readBin(con, "raw", file.info(file_path)$size * 10)
    close(con)
    # Check for double-gzip
    if (length(raw_data) >= 2 && raw_data[1] == as.raw(0x1f) && raw_data[2] == as.raw(0x8b)) {
      message("Detected double-gzipped file, decompressing again...")
      temp_file <- tempfile()
      writeBin(raw_data, temp_file)
      con <- gzfile(temp_file, "rb")
      raw_data <- readBin(con, "raw", length(raw_data) * 10)
      close(con)
      unlink(temp_file)
    }
    temp_rds <- tempfile(fileext = ".rds")
    writeBin(raw_data, temp_rds)
    message("Loading RDS data...")
    obj <- readRDS(temp_rds)
    unlink(temp_rds)
  } else {
    obj <- readRDS(file_path)
  }
}

if (!inherits(obj, "Seurat")) {
  stop("Object is not a Seurat object")
}

n_cells <- ncol(obj)
n_genes <- nrow(obj)
message("PROGRESS:25:Loaded ", format(n_cells, big.mark=","), " cells x ", format(n_genes, big.mark=","), " genes")

# Save metadata
message("PROGRESS:30:Saving metadata...")
write.csv(obj@meta.data, paste0(output_prefix, "_metadata.csv"), row.names = TRUE)
message("PROGRESS:40:Metadata saved")

# Extract expression matrix from the RNA assay. The object's default
# assay is often SCT, whose "counts" slot holds SCT-corrected counts on a
# reduced gene set, not the raw counts a pseudobulk model needs.
message("PROGRESS:42:Extracting expression matrix...")
assay <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
message("Assays present: ", paste(names(obj@assays), collapse = ", "),
        " | exporting counts from: ", assay)
counts <- NULL
tryCatch({
  counts <- GetAssayData(obj, assay = assay, layer = "counts")
}, error = function(e) {
  message("Trying slot instead of layer...")
})
if (is.null(counts)) {
  counts <- GetAssayData(obj, assay = assay, slot = "counts")
}
if (any(counts@x != floor(counts@x))) {
  stop("The ", assay, " counts slot is not integer counts; refusing to export it.")
}

message("PROGRESS:50:Writing expression matrix...")
writeMM(counts, paste0(output_prefix, "_expression.mtx"))
message("PROGRESS:70:Expression matrix saved")

# Gene/cell names
message("PROGRESS:72:Saving gene and cell names...")
writeLines(rownames(counts), paste0(output_prefix, "_genes.txt"))
writeLines(colnames(counts), paste0(output_prefix, "_cells.txt"))
message("PROGRESS:75:Gene/cell names saved")

# Variable features
message("PROGRESS:76:Checking for variable features...")
tryCatch({
  var_features <- VariableFeatures(obj)
  if (length(var_features) > 0) {
    writeLines(var_features, paste0(output_prefix, "_variable_features.txt"))
    message("PROGRESS:80:Saved ", length(var_features), " variable features")
  }
}, error = function(e) {
  message("PROGRESS:80:No variable features found")
})

# Embeddings
message("PROGRESS:82:Checking for embeddings...")
reductions <- names(obj@reductions)
if (length(reductions) > 0) {
  message("Found: ", paste(reductions, collapse=", "))
  n_red <- length(reductions)
  for (i in seq_along(reductions)) {
    red_name <- reductions[i]
    pct <- 82 + floor((i / n_red) * 13)
    tryCatch({
      embeddings <- Embeddings(obj, reduction = red_name)
      out_file <- paste0(output_prefix, "_reduction_", red_name, ".csv")
      write.csv(embeddings, out_file, row.names = TRUE)
      message("PROGRESS:", pct, ":Saved ", red_name, " (", ncol(embeddings), " dims)")
    }, error = function(e) {
      message("PROGRESS:", pct, ":Could not export ", red_name)
    })
  }
} else {
  message("PROGRESS:95:No embeddings found in object")
}

message("PROGRESS:100:R export complete!")
message("Metadata columns: ", paste(colnames(obj@meta.data), collapse=", "))
"""

# Memory-efficient R script for convert_tab's LoadRDSWorker._convert_with_r()
_RDS_MEMORY_EFFICIENT_SCRIPT = r"""
# Memory-efficient Seurat to intermediate files conversion
# Arguments: <rds_path> <temp_dir>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript convert.R <rds_path> <temp_dir>")
}

file_path <- args[1]
temp_dir <- args[2]

options(warn = -1)
suppressPackageStartupMessages({ library(Matrix) })

gc_full <- function() {
  for(i in 1:3) gc(verbose = FALSE, full = TRUE)
}

cat("Memory-efficient RDS loading...\n")

is_gzipped <- grepl("\\.gz$", file_path, ignore.case = TRUE)

tryCatch({
  if (is_gzipped) {
    cat("Decompressing gzipped file...\n")
    con <- gzfile(file_path, "rb")
    raw_data <- readBin(con, "raw", file.info(file_path)$size * 10)
    close(con)
    if (length(raw_data) >= 2 && raw_data[1] == as.raw(0x1f) && raw_data[2] == as.raw(0x8b)) {
      cat("Detected double-gzipped file, decompressing again...\n")
      temp_file <- tempfile()
      writeBin(raw_data, temp_file)
      rm(raw_data); gc_full()
      con <- gzfile(temp_file, "rb")
      raw_data <- readBin(con, "raw", file.info(temp_file)$size * 10)
      close(con)
      unlink(temp_file)
    }
    cat("Loading object from decompressed data...\n")
    obj <- unserialize(raw_data)
    rm(raw_data); gc_full()
  } else {
    obj <- readRDS(file_path)
  }

  cat("Object class:", class(obj)[1], "\n")
  gc_full()

  if (inherits(obj, "Seurat")) {
    suppressPackageStartupMessages(library(Seurat))
    cat("Seurat object: ", ncol(obj), " cells x ", nrow(obj), " genes\n", sep="")

    # Try SeuratDisk if available
    if (requireNamespace("SeuratDisk", quietly = TRUE)) {
      cat("Using SeuratDisk for direct h5ad export...\n")
      library(SeuratDisk)
      h5seurat_path <- file.path(temp_dir, "temp.h5seurat")
      tryCatch({
        SaveH5Seurat(obj, filename = h5seurat_path, overwrite = TRUE)
        gc_full()
        rm(obj); gc_full()
        Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)
        unlink(h5seurat_path)
        cat("H5AD_SUCCESS\n")
      }, error = function(e) {
        cat("SeuratDisk failed:", conditionMessage(e), "\n")
        cat("Falling back to sparse matrix export...\n")
      })
    }

    if (!file.exists(file.path(temp_dir, "temp.h5ad"))) {
      cat("Extracting sparse matrix (memory-efficient)...\n")
      assay <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
      counts <- tryCatch({
        GetAssayData(obj, assay = assay, slot = "counts")
      }, error = function(e) {
        tryCatch({
          obj[["RNA"]]$counts
        }, error = function(e2) {
          obj@assays$RNA@counts
        })
      })

      if (!inherits(counts, "dgCMatrix")) {
        counts <- as(counts, "dgCMatrix")
      }

      cat("Matrix dimensions:", nrow(counts), "x", ncol(counts), "\n")

      gene_names <- rownames(counts)
      cell_names <- colnames(counts)

      cat("Saving metadata...\n")
      write.csv(obj@meta.data, file.path(temp_dir, "metadata.csv"), row.names = TRUE)
      rm(obj); gc_full()

      cat("Saving sparse matrix (MTX format)...\n")
      writeMM(counts, file.path(temp_dir, "matrix.mtx"))
      writeLines(gene_names, file.path(temp_dir, "genes.txt"))
      writeLines(cell_names, file.path(temp_dir, "barcodes.txt"))
      rm(counts); gc_full()

      cat("MTX_SUCCESS\n")
    }

  } else if (is.matrix(obj) || inherits(obj, "dgCMatrix") || inherits(obj, "Matrix")) {
    cat("Detected matrix object\n")
    if (!inherits(obj, "dgCMatrix")) {
      obj <- as(obj, "dgCMatrix")
    }
    writeMM(obj, file.path(temp_dir, "matrix.mtx"))
    if (!is.null(rownames(obj))) writeLines(rownames(obj), file.path(temp_dir, "genes.txt"))
    if (!is.null(colnames(obj))) writeLines(colnames(obj), file.path(temp_dir, "barcodes.txt"))
    cat("MTX_SUCCESS\n")
  } else {
    cat("Unknown object type:", class(obj)[1], "\n")
    cat("FAILED\n")
  }

}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  cat("FAILED\n")
})
"""


def generate_rds_export_script() -> str:
    """Return the R script template for Seurat RDS export.

    The script uses commandArgs(trailingOnly=TRUE) to receive paths,
    not f-string interpolation. This prevents injection attacks (Issue #8).
    """
    return _RDS_EXPORT_SCRIPT


# Strategy 1: rds2py

def load_rds_rds2py(rds_path: str, progress_callback: Optional[Callable] = None):
    """
    Try to load an RDS file using rds2py (handles Seurat objects directly in Python).

    Parameters
    ----------
    rds_path : str
        Path to the RDS file.
    progress_callback : callable, optional
        Called with status strings.

    Returns
    -------
    anndata.AnnData or None
        Loaded AnnData, or None if rds2py fails.
    """
    try:
        import rds2py
        from anndata import AnnData
        from scipy.sparse import csr_matrix

        if progress_callback:
            progress_callback("Trying rds2py (Seurat-aware Python reader)...")

        rds_obj = rds2py.read_rds(str(rds_path))

        # Already AnnData-like
        if hasattr(rds_obj, 'X') and hasattr(rds_obj, 'obs'):
            return rds_obj

        # Try to extract Seurat components from dict
        if isinstance(rds_obj, dict):
            if 'assays' in rds_obj:
                assays = rds_obj['assays']
                for assay_name in ['RNA', 'Spatial', 'SCT']:
                    if assay_name in assays:
                        assay = assays[assay_name]
                        counts = None
                        for slot in ['counts', 'data', 'layers']:
                            if slot in assay:
                                counts = assay[slot]
                                break

                        if counts is not None:
                            if hasattr(counts, 'toarray'):
                                X = counts
                            elif isinstance(counts, np.ndarray):
                                X = csr_matrix(counts)
                            else:
                                X = csr_matrix(np.array(counts))

                            cell_names = None
                            if 'meta.data' in rds_obj:
                                meta = rds_obj['meta.data']
                                if isinstance(meta, pd.DataFrame):
                                    cell_names = meta.index.tolist()
                                elif hasattr(meta, 'index'):
                                    cell_names = list(meta.index)

                            gene_names = None
                            if 'features' in assay:
                                gene_names = list(assay['features'])
                            elif hasattr(assay, 'var_names'):
                                gene_names = list(assay.var_names)

                            # Transpose if needed (Seurat: genes x cells)
                            if gene_names and X.shape[0] == len(gene_names):
                                X = X.T

                            adata = AnnData(X=X)
                            if cell_names:
                                adata.obs_names = pd.Index(cell_names[:adata.n_obs])
                            if gene_names:
                                adata.var_names = pd.Index(gene_names[:adata.n_vars])

                            if 'meta.data' in rds_obj and isinstance(rds_obj['meta.data'], pd.DataFrame):
                                adata.obs = rds_obj['meta.data'].loc[adata.obs_names]

                            return adata

        return None

    except ImportError:
        return None
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return None


# Strategy 2: pyreadr

def load_rds_pyreadr(rds_path: str, progress_callback: Optional[Callable] = None):
    """
    Try to load an RDS file using pyreadr (pure Python, simple RDS files).

    Parameters
    ----------
    rds_path : str
        Path to the RDS file.
    progress_callback : callable, optional
        Called with status strings.

    Returns
    -------
    anndata.AnnData or None
        Loaded AnnData, or None if pyreadr fails.
    """
    try:
        import pyreadr
        from anndata import AnnData
        from scipy.sparse import csr_matrix

        if progress_callback:
            progress_callback("Trying pyreadr (pure Python)...")

        result = pyreadr.read_r(str(rds_path))

        if result is None or len(result) == 0:
            return None

        for key, df in result.items():
            if isinstance(df, pd.DataFrame):
                if df.shape[0] > 100 and df.shape[1] > 100:
                    X = csr_matrix(df.values.T)
                    adata = AnnData(X)
                    adata.var_names = df.index.astype(str)
                    adata.obs_names = df.columns.astype(str)
                    return adata

        return None

    except ImportError:
        return None
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return None


# Strategy 3: R subprocess (Issue #8 fix)

def run_rds_to_intermediate(
    rds_path: str,
    output_dir: str,
    progress_callback: Optional[Callable] = None,
) -> Optional[Path]:
    """
    Use R subprocess to export a Seurat RDS to intermediate files.

    The R script receives paths via commandArgs, NOT via string interpolation
    in the script body (fixes Issue #8 — R script injection).

    Parameters
    ----------
    rds_path : str
        Path to the RDS file.
    output_dir : str
        Directory for intermediate files.
    progress_callback : callable, optional
        Called with (progress_pct, message) tuples.

    Returns
    -------
    Path or None
        Path to the intermediate directory, or None on failure.
    """

    r_exe = find_r_executable()
    if not r_exe:
        return None

    rds_path = Path(rds_path)
    intermediate_dir = Path(output_dir) / "intermediate"
    intermediate_dir.mkdir(exist_ok=True)

    # Derive prefix from filename
    rds_name = rds_path.stem
    if rds_name.lower().endswith('.robj'):
        rds_name = rds_name[:-5]
    elif rds_name.lower().endswith('.rds'):
        rds_name = rds_name[:-4]

    is_robj = '.robj' in str(rds_path).lower()

    # Write script to a temp file (no path interpolation!)
    script_path = Path(output_dir) / "_convert_rds.R"
    script_path.write_text(_RDS_EXPORT_SCRIPT)

    # Pass paths as command-line arguments (Issue #8 fix)
    process = subprocess.Popen(
        [r_exe, str(script_path), str(rds_path), str(intermediate_dir),
         rds_name, str(is_robj).lower()],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    error_lines = []
    for line in process.stderr:
        line = line.strip()
        if not line:
            continue
        if line.startswith("PROGRESS:"):
            parts = line.split(":", 2)
            if len(parts) >= 3:
                try:
                    pct = int(parts[1])
                    msg = parts[2]
                    if progress_callback:
                        progress_callback(f"R export ({pct}%): {msg}")
                except ValueError:
                    pass
        elif "Error" in line or "error" in line:
            error_lines.append(line)

    process.wait()

    if process.returncode != 0:
        return None

    return intermediate_dir


def load_rds_via_r(
    rds_path: str,
    output_dir: str,
    progress_callback: Optional[Callable] = None,
):
    """
    Convert an RDS file to AnnData using R subprocess with memory-efficient extraction.

    Uses commandArgs for path passing (Issue #8 fix).

    Parameters
    ----------
    rds_path : str
        Path to the RDS file.
    output_dir : str
        Working directory for temp files.
    progress_callback : callable, optional
        Called with status strings.

    Returns
    -------
    anndata.AnnData or None
        Loaded AnnData, or None on failure.
    """
    r_exe = find_r_executable()
    if not r_exe:
        if progress_callback:
            progress_callback("R not found. Please install R.")
        return None

    temp_dir = Path(output_dir) / f".temp_convert_{os.getpid()}"
    temp_dir.mkdir(exist_ok=True)

    # Write script to file, pass paths as args (Issue #8 fix)
    script_path = temp_dir / "_convert.R"
    script_path.write_text(_RDS_MEMORY_EFFICIENT_SCRIPT)

    try:
        if progress_callback:
            progress_callback("Converting with R (memory-efficient mode)...")

        result = subprocess.run(
            [r_exe, str(script_path), str(rds_path), str(temp_dir)],
            capture_output=True,
            text=True,
            timeout=1800,
        )

        output = result.stdout + result.stderr

        if "H5AD_SUCCESS" in output:
            import scanpy as sc
            h5ad_path = temp_dir / "temp.h5ad"
            if h5ad_path.exists():
                adata = sc.read_h5ad(h5ad_path)
                shutil.rmtree(temp_dir, ignore_errors=True)
                return adata

        elif "MTX_SUCCESS" in output:
            from scipy.io import mmread
            from scipy.sparse import csr_matrix
            from anndata import AnnData

            mtx_path = temp_dir / "matrix.mtx"
            genes_path = temp_dir / "genes.txt"
            barcodes_path = temp_dir / "barcodes.txt"
            meta_path = temp_dir / "metadata.csv"

            if mtx_path.exists():
                matrix = mmread(str(mtx_path))
                X = csr_matrix(matrix.T)

                gene_names = None
                if genes_path.exists():
                    with open(genes_path, 'r') as f:
                        gene_names = [line.strip() for line in f]

                cell_names = None
                if barcodes_path.exists():
                    with open(barcodes_path, 'r') as f:
                        cell_names = [line.strip() for line in f]

                adata = AnnData(X=X)
                if gene_names:
                    adata.var_names = pd.Index(gene_names)
                if cell_names:
                    adata.obs_names = pd.Index(cell_names)

                if meta_path.exists():
                    meta = pd.read_csv(meta_path, index_col=0)
                    if cell_names:
                        common_cells = adata.obs_names.intersection(meta.index)
                        if len(common_cells) > 0:
                            adata = adata[common_cells, :].copy()
                            adata.obs = meta.loc[adata.obs_names]

                shutil.rmtree(temp_dir, ignore_errors=True)
                return adata

        shutil.rmtree(temp_dir, ignore_errors=True)
        return None

    except subprocess.TimeoutExpired:
        shutil.rmtree(temp_dir, ignore_errors=True)
        return None
    except (subprocess.CalledProcessError, OSError, ValueError, KeyError):
        shutil.rmtree(temp_dir, ignore_errors=True)
        return None


# Orchestrator: tries all strategies

def load_rds(
    rds_path: str,
    output_path: str,
    progress_callback: Optional[Callable] = None,
) -> Tuple:
    """
    Load an RDS file using the best available strategy.

    Tries: rds2py -> pyreadr -> R subprocess.

    Parameters
    ----------
    rds_path : str
        Path to the RDS file.
    output_path : str
        Path to save the resulting h5ad file.
    progress_callback : callable, optional
        Called with status strings.

    Returns
    -------
    tuple of (AnnData, dict)
        The loaded AnnData and a summary dict.

    Raises
    ------
    RuntimeError
        If all strategies fail.
    """
    from kosmic.scrna.load.converters import build_adata_summary

    rds_path = Path(rds_path)
    output_path = Path(output_path)
    output_dir = str(output_path.parent)

    # Strategy 1: rds2py
    adata = load_rds_rds2py(str(rds_path), progress_callback)

    # Strategy 2: pyreadr
    if adata is None:
        adata = load_rds_pyreadr(str(rds_path), progress_callback)

    # Strategy 3: R subprocess
    if adata is None:
        if progress_callback:
            progress_callback("Complex RDS format detected, using R for conversion...")
        adata = load_rds_via_r(str(rds_path), output_dir, progress_callback)

    if adata is None:
        raise RuntimeError(
            "Failed to load RDS file. Check if R and Seurat are installed."
        )

    # Ensure unique names
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)

    summary = build_adata_summary(adata, output_path)
    return adata, summary
