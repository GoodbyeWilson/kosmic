# scRNA-seq analysis utilities, organised by GUI workflow tab.
#
# Each subpackage mirrors a workflow tab:
#
# - 'kosmic.scrna.load'     -- Load Data tab (download, format conversion, RDS)
# - 'kosmic.scrna.inspect'  -- Inspect tab (column detection, roles, gene harmonisation)
# - 'kosmic.scrna.qc'       -- QC tab (filter pipeline, normalisation, Scrublet)
# - 'kosmic.scrna.cluster'  -- Cluster tab (HVG, Harmony, neighbors, Leiden, t-SNE)
# - 'kosmic.scrna.annotate' -- Annotate sub-tab (marker scoring, reference annotation)
#
# Each subpackage's '__init__' re-exports its public functions so
# callers can write 'from kosmic.scrna.cluster import run_harmony' rather
# than reaching into the per-method file.
