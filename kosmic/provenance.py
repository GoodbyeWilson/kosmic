# Per-study provenance sidecar.
#
# An append-only JSON record of the settings each analysis stage ran with, plus
# a cheap structural fingerprint of the data for staleness detection. Lives next
# to the study's h5ad as 'provenance.json'.
#
# An absent sidecar means "no provenance": callers fall back to inference and
# nothing is flagged as stale. The sidecar is independent of the h5ad -- writing
# it never touches the data file.
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

SCHEMA_VERSION = 1
PROVENANCE_FILENAME = 'provenance.json'

# Structural fields that define the data DE / meta consume. The fingerprint
# token hashes only these, so re-saving identical data (or re-clustering, which
# DE ignores) does not falsely flag staleness, while re-QC / role changes /
# decontX do. It deliberately excludes the sample/condition column *names* (a
# DE-side choice, not data state) so scRNA and DE compute the same token; a
# changed role split is captured via 'role_counts'.
_TOKEN_FIELDS = ('n_obs', 'n_vars', 'n_raw_vars', 'layers', 'role_counts')


def provenance_path(study_dir) -> Path:
    """Sidecar path for a study directory."""
    return Path(study_dir) / PROVENANCE_FILENAME


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def _token(fingerprint: dict) -> str:
    payload = {k: fingerprint.get(k) for k in _TOKEN_FIELDS}
    blob = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha1(blob.encode('utf-8')).hexdigest()[:12]


def compute_fingerprint(adata, sample_col: Optional[str] = None,
                        condition_col: Optional[str] = None,
                        h5ad_path=None) -> dict:
    """Cheap structural fingerprint of the data DE / meta depend on.

    No content hash -- just dimensions, layers, the role/condition columns and
    role class counts, plus (informational only) the h5ad size/mtime. 'token' is
    a short hash of the structural fields; two states with the same structure and
    role split share a token.
    """
    obs_cols = sorted(
        c for c in (sample_col, condition_col, '_role')
        if c and c in adata.obs.columns
    )
    fp = {
        'n_obs': int(adata.n_obs),
        'n_vars': int(adata.n_vars),
        'n_raw_vars': int(adata.raw.n_vars) if adata.raw is not None else None,
        'layers': sorted(adata.layers.keys()),
        'obs_cols': obs_cols,
    }
    if '_role' in adata.obs.columns:
        vc = adata.obs['_role'].astype(str).value_counts().to_dict()
        fp['role_counts'] = {str(k): int(v) for k, v in sorted(vc.items())}

    if h5ad_path:
        p = Path(h5ad_path)
        if p.is_file():
            stat = p.stat()
            fp['h5ad_size'] = int(stat.st_size)
            fp['h5ad_mtime'] = int(stat.st_mtime)

    fp['token'] = _token(fp)
    return fp


def load(study_dir) -> Optional[dict]:
    """Read a study's provenance record, or None when absent/unreadable."""
    path = provenance_path(study_dir)
    if not path.is_file():
        return None
    try:
        with open(path, encoding='utf-8') as f:
            rec = json.load(f)
        return rec if isinstance(rec, dict) else None
    except (OSError, json.JSONDecodeError):
        return None


def _save(study_dir, record: dict) -> None:
    path = provenance_path(study_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(record, f, indent=2, default=str)


def record_stage(study_dir, study: str, stage: str, params: dict, *,
                 fingerprint: Optional[dict] = None,
                 consumed_token: Optional[str] = None,
                 source: Optional[str] = None,
                 replace: bool = False) -> dict:
    """Append a stage entry to the study's sidecar (creating it if absent).

    'fingerprint' (from compute_fingerprint) updates the record's current data
    state and stamps this stage's 'produced_token' -- pass it for data-changing
    stages (qc / normalize / cluster / decontx / roles). 'consumed_token' records
    the upstream token a downstream stage (gene_de / pathway_de / meta) was
    computed against, so later staleness checks can compare. 'source' is the
    path to the parent h5ad this study was derived from (e.g. a cell-type
    subset's full dataset); it is stored once so the lineage can be followed
    back to the upstream processing.

    'replace' drops any earlier entries of the same stage name before
    adding this one, so the record keeps only the latest. Append is the
    default and is right for a study, where each stage transforms the
    h5ad and the chain is a real lineage. It is wrong for an analysis
    that emits a file: re-running a meta-analysis overwrites its output
    but was appending another entry, so the record grew a click history
    describing results that no longer existed.
    """
    rec = load(study_dir) or {
        'schema_version': SCHEMA_VERSION,
        'study': study,
        'fingerprint': None,
        'stages': [],
    }
    if source is not None:
        rec['source'] = str(source)
    entry = {'stage': stage, 'timestamp': _now(), 'params': params}
    if fingerprint is not None:
        rec['fingerprint'] = fingerprint
        entry['produced_token'] = fingerprint.get('token')
    if consumed_token is not None:
        entry['consumed_token'] = consumed_token
    stages = rec.setdefault('stages', [])
    if replace:
        stages[:] = [st for st in stages if st.get('stage') != stage]
    stages.append(entry)
    _save(study_dir, rec)
    write_methods(study_dir)  # keep the readable companion in sync
    return rec


def load_lineage(study_dir) -> list:
    """Records along the derivation chain, oldest ancestor first.

    Follows each record's 'source' pointer (a parent h5ad path) back to the
    study it was derived from, so a cell-type subset resolves to the full
    dataset's processing. Returns a list of (study_dir, record) pairs; a study
    with no sidecar yields (study_dir, None). Cycle-guarded.
    """
    chain = []
    seen = set()
    cur = Path(study_dir)
    while cur is not None:
        key = str(cur)
        if key in seen:
            break
        seen.add(key)
        rec = load(cur)
        chain.append((cur, rec))
        src = (rec or {}).get('source')
        cur = Path(src).parent if src else None
    chain.reverse()
    return chain


def current_token(study_dir) -> Optional[str]:
    """The data token the sidecar currently records, or None."""
    rec = load(study_dir)
    if not rec:
        return None
    fp = rec.get('fingerprint') or {}
    return fp.get('token')


def last_stage(record_or_dir, stage_name: str) -> Optional[dict]:
    """The most recent entry for 'stage_name', or None. Accepts a loaded record
    or a study_dir."""
    rec = record_or_dir if isinstance(record_or_dir, dict) else load(record_or_dir)
    if not rec:
        return None
    for entry in reversed(rec.get('stages', [])):
        if entry.get('stage') == stage_name:
            return entry
    return None


def staleness(study_dir, adata, sample_col: Optional[str] = None,
              condition_col: Optional[str] = None, h5ad_path=None) -> str:
    """Compare the recorded fingerprint to the current data.

    Returns 'none' (no sidecar -> today's behaviour), 'match' (data unchanged),
    or 'stale' (data changed since the record was written).
    """
    recorded = current_token(study_dir)
    if recorded is None:
        return 'none'
    now = compute_fingerprint(adata, sample_col, condition_col, h5ad_path)
    return 'match' if now.get('token') == recorded else 'stale'


# --- Human-readable rendering --------------------------------------------

_STAGE_TITLES = {
    'load': 'Load', 'qc': 'Quality control', 'normalize': 'Normalisation',
    'cluster': 'Clustering', 'decontx': 'DecontX decontamination',
    'embedding': 'HVG / PCA / batch correction',
    'subset': 'Cell-type subset',
    'filter_dataset': 'Pre-analysis cell filter',
    'gene_names': 'Gene-name harmonisation',
    'setup': 'Sample / role setup', 'roles': 'Role assignment',
    'sample_conditions': 'Per-sample condition labelling',
    'gene_de': 'Gene differential expression',
    'pathway_de': 'Pathway differential expression',
    'meta_gene': 'Gene meta-analysis', 'meta_pathway': 'Pathway meta-analysis',
    'meta_enrichment': 'Pathway enrichment of the pooled genes',
    'meta_reproducibility': 'Cross-dataset reproducibility',
    'meta_loo': 'Leave-one-study-out validation',
}


def _fmt_value(v) -> str:
    if isinstance(v, float):
        return f"{v:.4g}"
    if isinstance(v, dict):
        return ", ".join(f"{k}={_fmt_value(x)}" for k, x in v.items()) or "(none)"
    if isinstance(v, (list, tuple)):
        return ", ".join(str(x) for x in v) if v else "(none)"
    return str(v)


# A methods document reports the settings that produced the result. The
# record keeps everything; these two helpers decide what is worth reading.

_VALUE_LIMIT = 200


def _abbrev_value(key: str, value, limit: int = _VALUE_LIMIT) -> str:
    """Render a parameter, summarising ones too long to read.

    A gene-name harmonisation records every rename it made -- 47,746
    characters on one line for one study, which was 78% of the whole
    document. The count is the useful part; the list belongs in
    'provenance.json', which is where the reader is sent.
    """
    text = _fmt_value(value)
    if len(text) <= limit:
        return text
    n = len(value) if isinstance(value, (list, tuple, dict, set)) else None
    head = text[:limit].rsplit(', ', 1)[0]
    if n is not None:
        return f"{head}, ... ({n:,} total -- full list in provenance.json)"
    return f"{head}, ... ({len(text):,} chars -- full value in provenance.json)"


def _final_stages(stages: list) -> list:
    """Keep the last run of each stage, noting how many there were.

    Re-running clustering or DE overwrites what came before, so listing
    every attempt as an equal makes the document describe settings that
    produced nothing. The count is kept so nothing is hidden silently.
    """
    counts = {}
    for st in stages:
        counts[st.get('stage')] = counts.get(st.get('stage'), 0) + 1
    latest = {}
    for st in stages:
        latest[st.get('stage')] = st
    out = []
    for st in stages:
        if latest.get(st.get('stage')) is st:
            entry = dict(st)
            entry['_n_runs'] = counts.get(st.get('stage'), 1)
            out.append(entry)
    return out


def render_methods(record: Optional[dict], *, full: bool = False) -> str:
    """Render a provenance record as readable, ordered methods text.

    Returns '' for an absent record so callers can fall back to their own
    (e.g. inferred) rendering.

    By default this reports the *final* configuration of each stage and
    summarises parameters too long to read. 'full=True' renders every
    stage and every value verbatim, for when the audit trail is what you
    want rather than a methods section.
    """
    if not record or not record.get('stages'):
        return ''
    lines = []
    study = record.get('study')
    if study:
        lines.append(f"Study: {study}")
    src = record.get('source')
    if src:
        lines.append(f"Derived from: {Path(src).name}")
    fp = record.get('fingerprint') or {}
    if fp:
        bits = []
        if fp.get('n_obs') is not None:
            bits.append(f"{fp['n_obs']:,} cells")
        if fp.get('n_vars') is not None:
            bits.append(f"{fp['n_vars']:,} genes")
        if fp.get('role_counts'):
            bits.append("roles " + _fmt_value(fp['role_counts']))
        if bits:
            lines.append("Data: " + " | ".join(bits))
        if fp.get('token'):
            lines.append(f"Data token: {fp['token']}")
    lines.append("")

    stages = (record['stages'] if full
              else _final_stages(record['stages']))
    for i, stage in enumerate(stages, 1):
        title = _STAGE_TITLES.get(stage.get('stage'), stage.get('stage', 'stage'))
        ts = stage.get('timestamp', '')
        n_runs = stage.get('_n_runs', 1)
        rerun = f"  [last of {n_runs} runs]" if n_runs > 1 else ""
        lines.append(f"{i}. {title}  ({ts}){rerun}")
        for key, val in (stage.get('params') or {}).items():
            rendered = (_fmt_value(val) if full
                        else _abbrev_value(key, val))
            lines.append(f"     {key}: {rendered}")
        for tok_key, label in (('produced_token', 'produced token'),
                               ('consumed_token', 'ran against token')):
            if stage.get(tok_key):
                lines.append(f"     {label}: {stage[tok_key]}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def render_combined(study_records, meta_record: Optional[dict] = None, *,
                    full: bool = False) -> str:
    """Combined methods document across studies plus the meta-analysis step.

    'study_records' is an iterable of (label, record) pairs, where record is a
    loaded provenance dict or None. 'meta_record' is the meta-analysis sidecar
    (its stages describe the pooling). Returns '' when there is nothing to show.
    """
    blocks = []
    for label, rec in study_records:
        body = render_methods(rec, full=full)
        head = f"===== Study: {label} =====\n"
        blocks.append(head + (body or "(no provenance recorded for this study)\n"))
    meta_body = render_methods(meta_record, full=full)
    if meta_body:
        blocks.append("===== Meta-analysis =====\n" + meta_body)
    doc = "\n".join(blocks).strip()
    return doc + "\n" if doc else ""


def write_methods(study_dir) -> Optional[Path]:
    """Render the study's provenance to a readable 'methods.md' next to the
    sidecar. No-op (returns None) when there is no provenance yet."""
    rec = load(study_dir)
    text = render_methods(rec)
    if not text:
        return None
    path = Path(study_dir) / 'methods.md'
    try:
        path.write_text(text, encoding='utf-8')
        return path
    except OSError:
        return None
