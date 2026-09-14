"""Tests for the preranked GSEA (fgsea) leg.

Pins that `run_fgsea` recovers a planted enrichment (a pathway whose genes
sit at the top of the DESeq2 ranking gets a high positive NES and small p),
a null pathway does not, and inputs are validated.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from kosmic.de.fgsea import run_fgsea


def _gene_de(seed=0):
    rng = np.random.default_rng(seed)
    n = 1500
    genes = [f"G{i:04d}" for i in range(n)]
    lfc = rng.normal(0, 1, n)
    up = genes[:40]          # pathway planted at the top of the ranking
    lfc[:40] += 4.0
    return pd.DataFrame({"names": genes, "logfoldchanges": lfc,
                         "pvals": rng.uniform(0, 1, n)}), up


def test_fgsea_recovers_planted_pathway():
    de, up = _gene_de()
    gene_sets = {"UP": up, "RANDOM": [f"G{i:04d}" for i in range(700, 740)]}
    res = run_fgsea(de, gene_sets, permutation_num=200, seed=0)

    assert set(res["names"]) == {"UP", "RANDOM"}
    up_row = res[res["names"] == "UP"].iloc[0]
    rand_row = res[res["names"] == "RANDOM"].iloc[0]
    # Planted set: strongly positive NES, significant.
    assert up_row["nes"] > 1.5
    assert up_row["pvals"] < 0.05
    # Random set: not significant.
    assert rand_row["pvals"] > 0.1


def test_fgsea_case_insensitive_matching():
    de, up = _gene_de()
    # lower-case the gene-set genes; matching should still work.
    res = run_fgsea(de, {"UP": [g.lower() for g in up]},
                    permutation_num=100, seed=0)
    assert res[res["names"] == "UP"].iloc[0]["nes"] > 1.5


def test_fgsea_return_details():
    de, up = _gene_de()
    res, details = run_fgsea(de, {"UP": up}, permutation_num=100, seed=0,
                             return_details=True)
    assert "UP" in details
    d = details["UP"]
    assert len(d["RES"]) == len(de)          # running ES across all ranked genes
    assert len(d["hits"]) > 0                 # the set's genes have rank positions
    assert "nes" in d and "es" in d


def test_fgsea_empty_input_raises():
    with pytest.raises(ValueError):
        run_fgsea(pd.DataFrame(), {"PW": ["A", "B", "C"]})
