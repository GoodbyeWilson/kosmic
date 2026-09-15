# Third-party software and data

KOSMIC is licensed under the GNU General Public License v3.0 (see
`LICENSE`). This file records the licences of the software it depends
on and the terms attached to the reference data bundled in
`kosmic/reference/`.

## Why GPL-3.0

KOSMIC's GUI is built on **PyQt6**, which is available under the GPL-3.0
or a commercial licence from Riverbank Computing, and its clustering
uses **leidenalg** and **igraph**, which are GPL. A distributed build of
KOSMIC that includes these libraries is therefore a GPL-licensed work
whatever licence the KOSMIC code itself carries. Licensing KOSMIC under
GPL-3.0 states that plainly. Anyone wishing to ship KOSMIC inside a
proprietary product would need a commercial PyQt licence and would need
to replace the GPL clustering libraries; that is a constraint of the
dependencies, not a choice of this project.

## Software dependencies

Installed from PyPI; licences as declared by each package.

| Package | Licence |
|---|---|
| PyQt6 | GPL-3.0 / commercial |
| leidenalg | GPL-3.0-or-later |
| python-igraph | GPL-2.0 |
| scanpy, anndata | BSD-3-Clause |
| pydeseq2 | MIT |
| harmonypy | MIT |
| gseapy | BSD-3-Clause |
| goatools | BSD |
| numpy, scipy, pandas, matplotlib, scikit-learn, statsmodels, h5py, pyqtgraph, umap-learn, pynndescent, torch | BSD / MIT |
| celltypist (optional) | MIT |
| decontx (optional) | see package |

Ambient-RNA correction in `kosmic/scrna/qc/soupx.py` is an independent
re-implementation of the SoupX method (Young & Behjati, *GigaScience*
2020) and contains no SoupX code.

## Bundled reference data

Every file under `kosmic/reference/` is listed here with its source and
the terms that apply to it. Redistribution here is under the source's
own terms, which take precedence over the GPL for that data.

### Gene nomenclature and annotation

| File | Source | Terms |
|---|---|---|
| `hgnc/hgnc_symbols.tsv`, `hgnc/hgnc_locus_groups.tsv.gz` | HGNC (genenames.org) | CC0 |
| `go/go-basic.obo` | Gene Ontology Consortium | CC BY 4.0 |
| `go/gene2go_human.gz`, `go/hgnc_entrez.tsv` | NCBI Gene | public domain (US government work) |
| `orthologs/mgi_human_mouse_1to1.tsv` | Mouse Genome Informatics (MGI) | CC BY 4.0 |
| `gwas/gene_positions_GRCh37.tsv` | UCSC refFlat (hg19); built by `dev/scripts/generate_gene_positions.py` | free for all uses |

### Cell-type markers and atlases

| File | Source | Terms |
|---|---|---|
| `marker_dbs/PanglaoDB_markers.tsv` | PanglaoDB (Franzén, Gan & Björkegren, *Database* 2019) | The site states the data are public but publishes no licence. Academic use is the stated purpose; contact the authors before commercial redistribution. |
| `marker_dbs/Cell_marker_Seq.xlsx` | CellMarker 2.0 (Hu et al., *Nucleic Acids Research* 2023) | Freely downloadable; no licence text published. Same caution as PanglaoDB. |
| `atlases/heartmap_lv_broad.json.gz` | Per-cell-type mean expression profiles (13 types, 2,000 genes) derived from HeartMap (Datar et al., *Nat Cardiovasc Res* 2026; Broad Single Cell Portal SCP3689) | Derived summary, not the data. The SCP deposit is public; check the study's terms before commercial redistribution. |
| `atlases/gao_lv_broad.json.gz` | Per-cell-type mean expression profiles (13 types, 2,000 genes) derived from the Gao/Wu human LV atlas (Gao et al., *Genome Biol* 2026; GEO GSE290367) | Derived summary, not the data. GEO deposits are public. |
| `cell_type_focus/cardiac.json`, `cell_type_focus/stroke.json` | Lists of cell-type names selected from the two marker databases above | Names only. |
| `contamination/cardiomyocyte.json` | Derived from GSE183852 (Koenig et al., *Nature Cardiovascular Research* 2022): 34 cardiomyocyte-restricted genes chosen by expression and specificity thresholds recorded in the file | Derived gene list; the GEO deposit is public. |

Reference atlases for label transfer beyond the bundled heart atlas are
built on demand by `dev/scripts/build_all_references.py` from Tabula
Sapiens and Tabula Muris data (CC BY 4.0) and are not shipped.

### Gene sets and pathways

| File | Source | Terms |
|---|---|---|
| `gene_sets/ion_channels.gmt` | IUPHAR/BPS Guide to Pharmacology (GtoPdb 2026.2) families, mapped to HGNC symbols | GtoPdb: CC BY-SA 4.0 |
| `gene_sets/ion_channels_hgnc_original.gmt`, `gene_sets/gpcr.gmt`, `gene_sets/nicotinic.gmt` | HGNC gene groups (groups 177, 139 and the nicotinic receptor subunits) | CC0 |
| `gene_sets/mechanosensing.gmt` | Gene Ontology terms, hand-selected | CC BY 4.0 |
| `gene_sets/etc_complex.gmt` | Electron-transport-chain complexes I–V, curated by the authors | Part of KOSMIC (GPL-3.0) |
| `gene_sets/metabolic_kirk.gmt`, `gene_sets/metabolic_comprehensive*.gmt` | Curated by the authors for cardiac energy metabolism | Part of KOSMIC (GPL-3.0) |
| `gene_sets/metabolic_kegg_scmetabolism.gmt` | 85 KEGG metabolic pathways as redistributed by scMetabolism (Wu et al., *Cancer Discovery* 2022; GPL-3.0 R package) | **KEGG content.** KEGG is copyright Kanehisa Laboratories; academic use is free, non-academic use requires a KEGG licence. |
| `gene_sets/calcium_signaling.gmt` | KEGG hsa04020 (Calcium signaling pathway), fetched through Enrichr's KEGG_2021_Human library on 2026-09-03 | **KEGG content**; same terms as above. |

The two KEGG-derived files are included for the convenience of academic
users. If you use KOSMIC in a non-academic setting, KEGG's terms apply
to those two files independently of KOSMIC's licence.

### Proteomics panels

| File | Source | Terms |
|---|---|---|
| `olink_panels/Olink_CVD_panels.json`, `olink_panels/Olink_Explore_panels.json` | Protein lists of the Olink Target 96 CVD II/III and Olink Explore 3072 panels, from Olink's published panel contents | Lists of assay targets (facts), transcribed for the validation step. Olink® is a trademark of Olink Proteomics AB. |
| `olink_panels/Tromp2018_CVD_II_panel_genes.json` | Tromp et al., *European Journal of Heart Failure* 2018: the CVD II panel genes as used in that study | Gene list from a published paper. Not yet used by the application; kept for planned validation work. |
| `olink_panels/Wang2024_DCM_Olink_S1.xlsx` | Wang et al. 2024, supplementary table S1 (DCM plasma proteomics) | Published supplementary data. Not yet used by the application; kept for planned validation work. Cite the paper if you use the comparison. |

### Tutorial data

The tutorial dataset is **synthetic**: `kosmic/gui/help/tutorial/data.py`
generates it (six simulated donors, five cell types) on first use into
`tutorial/`. It contains no measured data, carries no third-party terms,
and is not stored in the repository.
