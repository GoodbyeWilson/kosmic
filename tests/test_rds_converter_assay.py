"""The Seurat import must export the RNA assay's raw counts, not the
default assay's.

A deposited Seurat object often has SCT as its default assay. SCT's
"counts" slot holds corrected counts on a reduced gene set, and
exporting it silently changes the count source of every downstream
step. Koenig (GSE183852) is the case that found this: default assay
SCT with 29,484 genes against 33,694 in RNA.

Needs R with Seurat installed; skipped otherwise.
"""
from __future__ import annotations

import subprocess

import pytest

from kosmic.scrna.load.rds_converter import (
    find_r_executable, run_rds_to_intermediate,
)


def _seurat_available(r_exe):
    try:
        out = subprocess.run(
            [r_exe, "-e", "cat(requireNamespace('Seurat', quietly=TRUE))"],
            capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return "TRUE" in out.stdout


_MAKE_OBJECT = r"""
suppressPackageStartupMessages(library(Seurat))
set.seed(1)
counts <- matrix(rpois(40 * 12, 3), nrow = 40,
                 dimnames = list(paste0("G", 1:40), paste0("c", 1:12)))
obj <- CreateSeuratObject(counts, assay = "RNA")
# A second assay with fewer genes, made the default -- the SCT situation.
sub <- counts[1:25, ]
obj[["SCT"]] <- CreateAssayObject(counts = sub)
DefaultAssay(obj) <- "SCT"
saveRDS(obj, commandArgs(trailingOnly = TRUE)[[1]])
"""


def test_export_uses_rna_assay_not_default(tmp_path):
    r_exe = find_r_executable()
    if not r_exe or not _seurat_available(r_exe):
        pytest.skip("R with Seurat not available")

    script = tmp_path / "make.R"
    script.write_text(_MAKE_OBJECT)
    rds = tmp_path / "toy.rds"
    subprocess.run([r_exe, str(script), str(rds)], check=True,
                   capture_output=True, text=True, timeout=600)

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    messages = []
    inter = run_rds_to_intermediate(str(rds), str(out_dir),
                                    progress_callback=messages.append)
    assert inter is not None, messages

    genes = (inter / "toy_genes.txt").read_text().split()
    assert len(genes) == 40, "exported the default (SCT) assay, not RNA"
