# Synthetic h5ad generator for the in-app tutorial. Simulates a cardiac
# scRNA-seq experiment: 6 patients (3 control, 3 disease) x 5 cell types.
from importlib.resources import files
from pathlib import Path
import numpy as np
import scipy.sparse as sp

# Gene lists

# Metabolic pathway genes (snapshot used by the tutorial walkthrough).
METABOLIC_GENES = {
    'Glycolysis_Gluconeogenesis': [
        'HK1', 'HK2', 'HK3', 'GCK', 'GPI', 'PFKL', 'PFKM', 'PFKP', 'ALDOA', 'ALDOB', 'ALDOC',
        'TPI1', 'GAPDH', 'PGK1', 'PGK2', 'PGAM1', 'PGAM2', 'ENO1', 'ENO2', 'ENO3',
        'PKM', 'PKLR', 'LDHA', 'LDHB', 'LDHC', 'G6PC', 'PCK1', 'PCK2', 'FBP1', 'FBP2',
    ],
    'TCA_Cycle': [
        'CS', 'ACO1', 'ACO2', 'IDH1', 'IDH2', 'IDH3A', 'IDH3B', 'IDH3G', 'OGDH', 'DLST',
        'SUCLA2', 'SUCLG1', 'SUCLG2', 'SDHA', 'SDHB', 'SDHC', 'SDHD', 'FH', 'MDH1', 'MDH2',
    ],
    'Oxidative_Phosphorylation': [
        'MT-ND1', 'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND4L', 'MT-ND5', 'MT-ND6',
        'MT-CYB', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-ATP6', 'MT-ATP8',
        'NDUFA1', 'NDUFA2', 'NDUFA3', 'NDUFA4', 'NDUFA5', 'NDUFA6', 'NDUFA7', 'NDUFA8', 'NDUFA9', 'NDUFA10',
        'NDUFB1', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFB8', 'NDUFB9', 'NDUFB10',
        'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS5', 'NDUFS6', 'NDUFS7', 'NDUFS8',
        'COX4I1', 'COX4I2', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1', 'COX6C', 'COX7A1', 'COX7A2', 'COX7B', 'COX7C', 'COX8A',
        'ATP5F1A', 'ATP5F1B', 'ATP5F1C', 'ATP5F1D', 'ATP5F1E', 'ATP5PB', 'ATP5MC1', 'ATP5MC2', 'ATP5MC3',
        'ATP5PD', 'ATP5PF', 'ATP5PO', 'ATP5IF1', 'ATP5MG', 'ATP5ME',
    ],
    'Fatty_Acid_Oxidation': [
        'CPT1A', 'CPT1B', 'CPT1C', 'CPT2', 'ACSL1', 'ACSL3', 'ACSL4', 'ACSL5', 'ACSL6',
        'ACADM', 'ACADL', 'ACADVL', 'ACADS', 'HADHA', 'HADHB', 'ECHS1', 'HADH',
        'ACAA1', 'ACAA2', 'ECI1', 'ECI2', 'DECR1',
    ],
    'Fatty_Acid_Synthesis': [
        'ACACA', 'ACACB', 'FASN', 'SCD', 'ELOVL1', 'ELOVL2', 'ELOVL3', 'ELOVL4', 'ELOVL5', 'ELOVL6', 'ELOVL7',
        'FADS1', 'FADS2', 'FADS3', 'MCAT', 'ACLY',
    ],
    'Pentose_Phosphate_Pathway': [
        'G6PD', 'PGLS', 'PGD', 'RPIA', 'RPE', 'TKT', 'TALDO1', 'PRPS1', 'PRPS2',
    ],
    'Glutamine_Metabolism': [
        'GLS', 'GLS2', 'GLUL', 'GLUD1', 'GLUD2', 'GOT1', 'GOT2', 'GPT', 'GPT2',
    ],
    'Amino_Acid_Metabolism': [
        'BCAT1', 'BCAT2', 'BCKDHA', 'BCKDHB', 'DBT', 'DLD',
        'PAH', 'TAT', 'HPD', 'HGD', 'FAH',
        'CBS', 'CTH', 'AHCY', 'MAT1A', 'MAT2A', 'MAT2B',
    ],
    'Ketone_Body_Metabolism': [
        'HMGCS1', 'HMGCS2', 'HMGCL', 'BDH1', 'BDH2', 'OXCT1', 'OXCT2', 'ACAT1', 'ACAT2',
    ],
    'Nucleotide_Metabolism': [
        'IMPDH1', 'IMPDH2', 'GMPS', 'ADSL', 'ADSS1', 'ADSS2', 'ATIC',
        'UMPS', 'CAD', 'DHODH', 'CTPS1', 'CTPS2',
    ],
}

# Cell type marker genes (from ClusterTab.DEFAULT_MARKERS)
MARKER_GENES = {
    'Endothelial': ['PECAM1', 'CDH5', 'VWF', 'KDR', 'FLT1', 'CLDN5', 'ERG'],
    'Cardiomyocyte': ['TNNT2', 'MYH7', 'MYH6', 'ACTC1', 'TTN', 'MYBPC3'],
    'Fibroblast': ['DCN', 'COL1A1', 'COL1A2', 'LUM', 'PDGFRA', 'THY1'],
    'Immune': ['CD68', 'CD14', 'FCGR3A', 'CSF1R', 'AIF1', 'MARCO'],
    'Smooth_Muscle': ['ACTA2', 'MYH11', 'TAGLN', 'CNN1', 'MYOCD'],
}

# Broader cell-type programs (moderate enrichment, ~2x) for UMAP separation
CELL_TYPE_PROGRAMS = {
    'Endothelial': [
        'ESAM', 'EMCN', 'CLEC14A', 'PLVAP', 'ENG', 'TIE1', 'TEK',
        'ROBO4', 'ADGRL4', 'ECSCR', 'FLT4', 'NOS3', 'THBD', 'PROCR',
        'ICAM2', 'VCAM1', 'SELE', 'SELP', 'ACKR1', 'DARC',
        'SOX17', 'SOX18', 'ETS1', 'ETS2', 'FLI1', 'EPAS1',
        'GJA4', 'GJA5', 'CALCRL', 'RAMP2', 'RAMP3', 'BMPR2',
        'NRP1', 'NRP2', 'DLL4', 'JAG1', 'HEY2', 'EFNB2',
        'APLN', 'APLNR', 'LYVE1', 'PROX1', 'PDPN', 'MMRN1',
        'VWA1', 'HSPG2', 'LAMA4', 'COL4A1', 'COL4A2',
    ],
    'Cardiomyocyte': [
        'RYR2', 'SCN5A', 'CASQ2', 'PLN', 'KCNJ2', 'KCNQ1', 'KCNH2',
        'CACNA1C', 'SLC8A1', 'ATP2A2', 'SERCA2', 'TNNI3', 'TNNC1',
        'MYL2', 'MYL3', 'MYL4', 'MYL7', 'DES', 'DMD', 'DTNA',
        'CSRP3', 'LDB3', 'TCAP', 'MYOM1', 'MYOM2', 'OBSCN',
        'NEBL', 'ANKRD1', 'ANKRD2', 'XIRP1', 'XIRP2',
        'HAND1', 'HAND2', 'NKX2-5', 'GATA4', 'TBX5', 'MEF2C',
        'NPPA', 'NPPB', 'BMP10', 'MYH14', 'ACTN2', 'SYNPO2L',
        'JPH2', 'TRDN', 'CALM1', 'CAMK2D', 'SLC25A4', 'CKM',
        'CKMT2', 'MB',
    ],
    'Fibroblast': [
        'FN1', 'VIM', 'S100A4', 'POSTN', 'SPARC', 'FBLN1', 'FBLN2',
        'FBLN5', 'ELN', 'BGN', 'FMOD', 'OGN', 'ASPN', 'OMD',
        'PCOLCE', 'PCOLCE2', 'LOX', 'LOXL1', 'LOXL2', 'LOXL3',
        'COL3A1', 'COL5A1', 'COL5A2', 'COL6A1', 'COL6A2', 'COL6A3',
        'COL12A1', 'COL14A1', 'COL15A1', 'MMP2', 'MMP3', 'MMP14',
        'TIMP1', 'TIMP3', 'SERPINE1', 'SERPINE2', 'TGFBI', 'CTGF',
        'IGFBP3', 'IGFBP5', 'IGFBP7', 'PDGFRB', 'PDGFRL',
        'CXCL12', 'CXCL14', 'CCL2', 'IL6', 'WNT2', 'WNT5A',
        'TCF21', 'TWIST1', 'SNAI2',
    ],
    'Immune': [
        'PTPRC', 'ITGAM', 'CD86', 'IL1B', 'TNF', 'CCL3', 'CCL4',
        'CCL5', 'CXCL8', 'CXCL10', 'IL6', 'IL10', 'TGFB1',
        'FCGR1A', 'FCGR2A', 'FCGR2B', 'CD163', 'MRC1', 'MSR1',
        'SIGLEC1', 'CD209', 'CLEC7A', 'CLEC10A', 'TLR2', 'TLR4',
        'MYD88', 'IRAK4', 'TRAF6', 'IRF5', 'IRF7', 'IRF8',
        'SPI1', 'CEBPB', 'MAFB', 'BATF3', 'ID2', 'ZBTB46',
        'CD74', 'HLA-DRA', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1',
        'CTSS', 'CTSB', 'CTSL', 'LYZ', 'TYROBP', 'FCER1G',
        'SYK', 'ITGAX', 'ITGB2',
    ],
    'Smooth_Muscle': [
        'CALD1', 'DES', 'LMOD1', 'SYNPO2', 'SMTN', 'CARMN',
        'PLN', 'PRKG1', 'GUCY1A1', 'GUCY1B1', 'PDE5A', 'PDE3A',
        'KCNMA1', 'KCNMB1', 'CACNA1C', 'RGS5', 'NOTCH3', 'JAG1',
        'PDGFRB', 'ABCC9', 'KCNJ8', 'ANGPT1', 'ADIRF', 'RERGL',
        'FOXF1', 'FOXF2', 'MYOCD', 'SRF', 'TEAD1', 'TEAD3',
        'FLNA', 'FLNC', 'TPM1', 'TPM2', 'CFL2', 'FHOD3',
        'AEBP1', 'IGFBP2', 'GJA1', 'GJC1', 'PTN', 'MDK',
        'NTRK3', 'MEF2C', 'MEIS1', 'MEIS2', 'PBX1',
        'ACTN1', 'VCL', 'TLN1',
    ],
}

# Housekeeping / background genes
HOUSEKEEPING_GENES = [
    'ACTB', 'ACTG1', 'TUBA1A', 'TUBA1B', 'TUBB', 'TUBB4B',
    'B2M', 'UBC', 'UBB', 'UBA52',
    'HSP90AA1', 'HSP90AB1', 'HSPA1A', 'HSPA1B', 'HSPA5', 'HSPA8', 'HSPB1', 'HSPD1',
    'EEF1A1', 'EEF2', 'EIF4A1', 'EIF4G1',
    'CALM1', 'CALM2', 'CALM3',
    'PPIA', 'PPIB', 'PPID',
    'HMGB1', 'HMGB2',
    'YWHAB', 'YWHAE', 'YWHAG', 'YWHAH', 'YWHAZ',
    'RAC1', 'CDC42', 'RHOA',
    'RACK1', 'GNB1', 'GNB2',
    'FOS', 'JUN', 'JUNB', 'JUND', 'MYC',
    'TP53', 'RB1', 'CDKN1A', 'CDKN2A',
    'VEGFA', 'FGF2', 'PDGFB', 'TGFB1', 'EGF',
    'SOD1', 'SOD2', 'CAT', 'GPX1', 'GPX4',
    'MMP2', 'MMP9', 'MMP14', 'TIMP1', 'TIMP2',
    'CASP3', 'CASP9', 'BCL2', 'BAX', 'BID',
    'STAT1', 'STAT3', 'JAK1', 'JAK2',
    'NFKB1', 'RELA', 'IKBKB',
    'AKT1', 'AKT2', 'MTOR', 'PIK3CA',
    'MAPK1', 'MAPK3', 'MAPK14', 'MAP2K1',
    'NOTCH1', 'NOTCH2', 'HES1', 'HEY1',
    'WNT3A', 'WNT5A', 'CTNNB1', 'APC',
    'SHH', 'PTCH1', 'SMO', 'GLI1',
]

# Ribosomal protein genes (RPL/RPS families)
_RPL_GENES = [f'RPL{i}' for i in list(range(3, 42)) + ['3A', '7A', '7L1', '10A', '13A', '22L1', '26L1', '27A', '29', '30', '31', '32', '34', '35', '35A', '36', '36A', '37', '37A', '38', '39', '39L', '40', '41']]
_RPS_GENES = [f'RPS{i}' for i in list(range(2, 30)) + ['3A', '4X', '4Y1', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '15A', '16', '17', '18', '19', '20', '21', '23', '24', '25', '26', '27', '27A', '28', '29']]


# Dataset generation

# Cell type distribution per patient (400 cells total)
_CELL_TYPE_COUNTS = {
    'Endothelial': 100,
    'Fibroblast': 80,
    'Immune': 80,
    'Smooth_Muscle': 70,
    'Cardiomyocyte': 70,
}

# Disease fold changes for specific pathways
_DISEASE_FOLD_CHANGES = {
    'Glycolysis_Gluconeogenesis': 1.8,
    'Oxidative_Phosphorylation': 1.8,
    'Fatty_Acid_Oxidation': 0.6,
}


def generate_tutorial_dataset():
    """Generate a realistic dummy scRNA-seq AnnData object.

    Returns
    -------
    anndata.AnnData
        2400 cells x ~1600 genes with obs columns:
        sample, condition, cell_type
    """
    import anndata as ad
    import pandas as pd

    rng = np.random.default_rng(42)

    # --- Build gene list (deduplicated, ordered) ---
    all_genes = []
    seen = set()
    # Pathway genes first
    for pathway_genes in METABOLIC_GENES.values():
        for g in pathway_genes:
            if g not in seen:
                all_genes.append(g)
                seen.add(g)
    # Marker genes
    for marker_genes in MARKER_GENES.values():
        for g in marker_genes:
            if g not in seen:
                all_genes.append(g)
                seen.add(g)
    # Cell type program genes
    for program_genes in CELL_TYPE_PROGRAMS.values():
        for g in program_genes:
            if g not in seen:
                all_genes.append(g)
                seen.add(g)
    # Housekeeping
    for g in HOUSEKEEPING_GENES:
        if g not in seen:
            all_genes.append(g)
            seen.add(g)
    # Ribosomal — deduplicate and sort
    for g in sorted(set(_RPL_GENES + _RPS_GENES)):
        if g not in seen:
            all_genes.append(g)
            seen.add(g)

    # Pad to ~1600 with additional background genes
    extra_genes = [
        # Collagen family
        *[f'COL{t}A{n}' for t in range(1, 18) for n in range(1, 4)],
        # SLC transporters
        *[f'SLC{fam}A{m}' for fam in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 17, 22, 25, 26, 27, 29, 30, 35, 38, 39, 44] for m in range(1, 10)],
        # Zinc finger proteins
        *[f'ZNF{i}' for i in range(1, 200)],
        # Interleukins
        *[f'IL{i}' for i in range(1, 40)],
        # Keratins
        *[f'KRT{i}' for i in range(1, 30)],
        # HOX genes
        *[f'HOX{c}{n}' for c in 'ABCD' for n in range(1, 14)],
        # ADAM/ADAMTS
        *[f'ADAM{i}' for i in range(1, 34)],
        *[f'ADAMTS{i}' for i in range(1, 21)],
        # Claudins
        *[f'CLDN{i}' for i in range(1, 26)],
        # Cadherins
        *[f'CDH{i}' for i in range(1, 24)],
        # G-protein coupled receptors
        *[f'GPR{i}' for i in range(1, 90)],
        # Solute carriers additional
        *[f'SLC{fam}A{m}' for fam in [11, 13, 14, 15, 18, 19, 20, 23, 24, 28, 31, 33, 34, 36, 37, 40, 41, 43, 45, 46, 47, 48, 51, 52] for m in range(1, 6)],
        # Cytochrome P450
        *[f'CYP{f}{s}{m}' for f in ['1', '2', '3', '4', '7', '11', '17', '19', '24', '26', '27'] for s in ['A', 'B', 'C', 'D'] for m in range(1, 4)],
        # Olfactory receptors (to fill remaining)
        *[f'OR{i}A{j}' for i in range(1, 15) for j in range(1, 10)],
    ]
    for g in extra_genes:
        if g not in seen:
            all_genes.append(g)
            seen.add(g)
        if len(all_genes) >= 1600:
            break

    n_genes = len(all_genes)
    gene_to_idx = {g: i for i, g in enumerate(all_genes)}

    # --- Build cell metadata ---
    patients = [
        ('Patient_1', 'control'), ('Patient_2', 'control'), ('Patient_3', 'control'),
        ('Patient_4', 'disease'), ('Patient_5', 'disease'), ('Patient_6', 'disease'),
    ]

    obs_records = []
    for sample, condition in patients:
        for ct, count in _CELL_TYPE_COUNTS.items():
            for _ in range(count):
                obs_records.append({'sample': sample, 'condition': condition, 'cell_type': ct})

    n_cells = len(obs_records)  # 2400
    obs = pd.DataFrame(obs_records)
    obs.index = [f'cell_{i}' for i in range(n_cells)]

    # --- Build expression matrix ---
    # Base: negative binomial counts with ~50 % zeros
    base = rng.negative_binomial(n=2, p=0.33, size=(n_cells, n_genes)).astype(np.float32)
    # Introduce sparsity
    dropout_mask = rng.random((n_cells, n_genes)) < 0.50
    base[dropout_mask] = 0.0

    # --- Apply broader cell-type-specific gene programs (~2x enrichment) ---
    for ct, program_genes in CELL_TYPE_PROGRAMS.items():
        cell_mask = (obs['cell_type'] == ct).values
        n_ct_cells = int(cell_mask.sum())
        for gene in program_genes:
            if gene in gene_to_idx:
                gi = gene_to_idx[gene]
                # Moderate enrichment: ~2x expression with reduced dropout (30%)
                prog_expr = rng.negative_binomial(n=3, p=0.33, size=n_ct_cells).astype(np.float32)
                prog_dropout = rng.random(n_ct_cells) < 0.30
                prog_expr[prog_dropout] = 0.0
                base[cell_mask, gi] = prog_expr

    # --- Apply cell type marker enrichment (strong signal) ---
    for ct, markers in MARKER_GENES.items():
        cell_mask = (obs['cell_type'] == ct).values
        for gene in markers:
            if gene in gene_to_idx:
                gi = gene_to_idx[gene]
                # Strong marker expression with very low dropout
                marker_expr = rng.negative_binomial(n=5, p=0.25, size=int(cell_mask.sum())).astype(np.float32)
                marker_dropout = rng.random(int(cell_mask.sum())) < 0.10
                marker_expr[marker_dropout] = 0.0
                base[cell_mask, gi] = marker_expr

    # --- Apply disease fold changes for specific pathways ---
    disease_mask = (obs['condition'] == 'disease').values
    for pathway, fc in _DISEASE_FOLD_CHANGES.items():
        if pathway in METABOLIC_GENES:
            for gene in METABOLIC_GENES[pathway]:
                if gene in gene_to_idx:
                    gi = gene_to_idx[gene]
                    base[disease_mask, gi] = (base[disease_mask, gi] * fc).astype(np.float32)

    # --- Add patient-level variation (batch effect) ---
    for sample, _ in patients:
        sample_mask = (obs['sample'] == sample).values
        patient_scale = rng.normal(1.0, 0.1)
        base[sample_mask] = (base[sample_mask] * max(patient_scale, 0.5)).astype(np.float32)

    # Ensure non-negative integers for raw counts
    base = np.maximum(base, 0).astype(np.float32)

    # Convert to sparse CSR
    X = sp.csr_matrix(base)

    # --- Build AnnData ---
    var = pd.DataFrame(index=all_genes)
    var.index.name = 'gene'

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs['sample'] = adata.obs['sample'].astype('category')
    adata.obs['condition'] = adata.obs['condition'].astype('category')
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

    return adata


def get_tutorial_dataset_path() -> Path:
    """Return the default path for the tutorial dataset."""
    return Path(str(files("kosmic.reference.tutorial") / "tutorial_dataset.h5ad"))


# Bump this version whenever the dataset generation logic changes
# to force regeneration of cached datasets.
_DATASET_VERSION = 2

def ensure_tutorial_dataset() -> Path:
    """Generate the tutorial dataset if it doesn't exist (or is outdated) and return its path."""
    path = get_tutorial_dataset_path()
    version_file = path.with_suffix('.version')

    # Regenerate if missing or version mismatch
    needs_regen = not path.exists()
    if not needs_regen and version_file.exists():
        try:
            stored_version = int(version_file.read_text().strip())
            if stored_version < _DATASET_VERSION:
                needs_regen = True
        except (ValueError, OSError):
            needs_regen = True
    elif not needs_regen and not version_file.exists():
        needs_regen = True

    if needs_regen:
        path.parent.mkdir(parents=True, exist_ok=True)
        adata = generate_tutorial_dataset()
        adata.write_h5ad(path)
        version_file.write_text(str(_DATASET_VERSION))

    return path
