# GEO Download & Archive Extraction
# Query GEO for supplementary files, download with progress, extract tar archives,
# and load multi-sample 10X data from extracted tarballs.

import re
import tarfile
import urllib.error
import urllib.request
from html.parser import HTMLParser
from pathlib import Path
from typing import Callable, Dict, List, Optional

from kosmic.scrna.load.converters import detect_file_format, load_10x_folder


# GEO query

class _GEOFileParser(HTMLParser):
    """Parse GEO FTP directory listing for filenames."""

    def __init__(self):
        super().__init__()
        self.files: List[str] = []

    def handle_data(self, data):
        data = data.strip()
        if data and ('.' in data) and not data.startswith('Parent'):
            self.files.append(data)


_DATA_EXTENSIONS = (
    '.rds', '.rds.gz', '.robj.gz', '.robj', '.h5ad', '.h5ad.gz',
    '.tar', '.tar.gz', '.h5', '.mtx.gz', '.csv.gz', '.txt.gz', '.tsv.gz',
)


def query_geo_files(gse: str) -> List[Dict]:
    """
    Query GEO for supplementary files for a given GSE accession.

    Parameters
    ----------
    gse : str
        GEO Series accession (e.g. 'GSE12345').

    Returns
    -------
    list of dict
        Each dict has keys: name, url, size, format.

    Raises
    ------
    urllib.error.HTTPError
        If the GSE number is invalid or server returns an error.
    """
    gse = gse.upper().strip()
    prefix = gse[:6]
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{gse}/suppl/"

    with urllib.request.urlopen(url, timeout=15) as response:
        html = response.read().decode('utf-8')

    # Parse HTML directory listing
    parser = _GEOFileParser()
    parser.feed(html)

    # Also extract from href tags
    files = parser.files.copy()
    for match in re.findall(r'href="([^"]+)"', html):
        if '.' in match and not match.startswith('?') and match not in files:
            files.append(match)

    # Filter to data files and get sizes
    data_files = []
    for f in files:
        f = f.strip()
        if any(f.lower().endswith(ext) for ext in _DATA_EXTENSIONS):
            size = 0
            try:
                file_url = url + f
                req = urllib.request.Request(file_url, method='HEAD')
                with urllib.request.urlopen(req, timeout=10) as resp:
                    size = int(resp.headers.get('content-length', 0))
            except (urllib.error.URLError, OSError, ValueError):
                pass

            data_files.append({
                'name': f,
                'url': url + f,
                'size': size,
                'format': detect_file_format(f),
            })

    return data_files


# HTTP download

def download_file(
    url: str,
    output_path: str,
    progress_callback: Optional[Callable[[int, int], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
    expected_size: int = 0,
) -> str:
    """
    Download a file from a URL with progress reporting.

    Parameters
    ----------
    url : str
        URL to download.
    output_path : str
        Local file path to save to.
    progress_callback : callable, optional
        Called with (downloaded_bytes, total_bytes) after each chunk.
    cancel_check : callable, optional
        Called before each chunk read; if returns True, download is cancelled.
    expected_size : int, optional
        File size in bytes obtained out-of-band (typically from a HEAD
        request before the GET). Used as a safety net when the GET
        response's 'content-length' is missing, zero, or implausibly
        small relative to 'expected_size' -- which happens with some
        GEO redirects / CDN handoffs and would otherwise cause the
        progress bar to flash to 100% within the first chunk.

    Returns
    -------
    str
        Path to the downloaded file.

    Raises
    ------
    Exception
        If download fails or is cancelled.
    """
    with urllib.request.urlopen(url) as response:
        total = int(response.headers.get('content-length', 0))
        # Prefer expected_size when the GET response advertises a size
        # that is missing or much smaller than the HEAD-reported size.
        if expected_size > 0 and (total <= 0 or total < expected_size * 0.9):
            total = expected_size
        downloaded = 0

        with open(output_path, 'wb') as f:
            while True:
                if cancel_check and cancel_check():
                    raise RuntimeError("Download cancelled")
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if progress_callback:
                    # Cap reported 'total' to at least the bytes already
                    # received so progress never exceeds 100% if the
                    # actual download outruns the advertised size.
                    progress_callback(downloaded, max(total, downloaded))

    return str(output_path)


# Archive extraction

def extract_archive(
    archive_path: str,
    output_dir: str,
    progress_callback: Optional[Callable[[str], None]] = None,
) -> List[str]:
    """
    Extract a tar archive (gz, bz2, xz, or plain).

    Parameters
    ----------
    archive_path : str
        Path to the archive file.
    output_dir : str
        Directory to extract into.
    progress_callback : callable, optional
        Called with status string after each member is extracted.

    Returns
    -------
    list of str
        Names of extracted members.
    """
    extracted = []
    with tarfile.open(archive_path, 'r:*') as tar:
        members = tar.getmembers()
        for i, member in enumerate(members):
            tar.extract(member, output_dir)
            extracted.append(member.name)
            if progress_callback:
                progress_callback(f"Extracting: {i + 1}/{len(members)}")
    return extracted


# Multi-sample 10X loading from tarballs

def find_sample_tarballs(raw_data_dir: str) -> List[Path]:
    """
    Find sample-level tar.gz files in a directory (e.g. GSM*_filtered_feature_bc_matrix.tar.gz).

    Filters out main RAW.tar and validates that tarballs contain 10X data.

    Parameters
    ----------
    raw_data_dir : str or Path
        Directory containing extracted tarballs.

    Returns
    -------
    list of Path
        Sorted list of valid tarball paths.
    """
    raw_data_dir = Path(raw_data_dir)

    tarballs = set()
    for pattern in [
        '*_filtered_feature_bc_matrix.tar.gz',
        '*_raw_feature_bc_matrix.tar.gz',
        'GSM*.tar.gz',
        '*_matrix.tar.gz',
    ]:
        tarballs.update(raw_data_dir.glob(pattern))

    # Filter out the main RAW.tar
    tarballs = [t for t in tarballs if '_RAW.tar' not in t.name and t.name != 'RAW.tar']

    # Validate: must contain matrix.mtx + barcodes.tsv
    valid = []
    for tb in tarballs:
        try:
            with tarfile.open(tb, 'r:gz') as tar:
                members = tar.getnames()
                has_matrix = any('matrix.mtx' in m for m in members)
                has_barcodes = any('barcodes.tsv' in m for m in members)
                if has_matrix and has_barcodes:
                    valid.append(tb)
        except (tarfile.TarError, OSError):
            continue

    return sorted(valid)


def extract_tarball(tarball_path: str, extract_dir: str) -> Path:
    """
    Extract a single sample tar.gz and return the extracted folder path.

    Handles both tarballs with subdirectories and flat tarballs.

    Parameters
    ----------
    tarball_path : str or Path
        Path to the tar.gz file.
    extract_dir : str or Path
        Directory to extract into.

    Returns
    -------
    Path
        Path to the extracted folder containing 10X files.
    """
    tarball_path = Path(tarball_path)
    extract_dir = Path(extract_dir)

    with tarfile.open(tarball_path, 'r:gz') as tar:
        members = tar.getnames()
        cleaned = [m.lstrip('./') for m in members if m.lstrip('./')]

        has_subdirectory = any('/' in m for m in cleaned)

        if has_subdirectory:
            top_dir = cleaned[0].split('/')[0]
            tar.extractall(extract_dir)
            return extract_dir / top_dir
        else:
            sample_name = tarball_path.stem.replace('.tar', '')
            sample_folder = extract_dir / sample_name
            sample_folder.mkdir(parents=True, exist_ok=True)
            tar.extractall(sample_folder)
            return sample_folder


def load_multi_sample_10x(
    raw_data_dir: str,
    output_path: str,
    progress_callback: Optional[Callable[[str], None]] = None,
):
    """
    Find, extract, and load multi-sample 10X tarballs into a combined AnnData.

    Parameters
    ----------
    raw_data_dir : str or Path
        Directory containing sample tar.gz files.
    output_path : str or Path
        Path to save the combined h5ad file.
    progress_callback : callable, optional
        Called with status messages.

    Returns
    -------
    anndata.AnnData
        Combined AnnData with 'sample' column in obs.
    """
    import scanpy as sc

    raw_data_dir = Path(raw_data_dir)
    output_path = Path(output_path)

    tarballs = find_sample_tarballs(raw_data_dir)
    if not tarballs:
        raise FileNotFoundError(
            "No sample tar.gz files found. Expected files like "
            "GSM*_filtered_feature_bc_matrix.tar.gz"
        )

    if progress_callback:
        progress_callback(f"Found {len(tarballs)} sample tar.gz file(s)")

    # Extract each tarball
    extracted_samples = []
    for tb in tarballs:
        if progress_callback:
            progress_callback(f"Extracting {tb.name}...")
        try:
            folder = extract_tarball(tb, raw_data_dir)
            if folder.exists():
                sample_name = tb.stem.replace('.tar', '')
                for suffix in ['_filtered_feature_bc_matrix', '_raw_feature_bc_matrix', '_matrix']:
                    sample_name = sample_name.replace(suffix, '')
                extracted_samples.append((folder, sample_name))
        except (tarfile.TarError, OSError, RuntimeError):
            continue

    if not extracted_samples:
        raise RuntimeError("Failed to extract any tar.gz files")

    # Load each sample
    adata_list = []
    sample_names = []
    for folder, sample_name in extracted_samples:
        if progress_callback:
            progress_callback(f"Loading {sample_name}...")
        try:
            adata = load_10x_folder(folder)
            adata.obs['sample'] = sample_name
            adata_list.append(adata)
            sample_names.append(sample_name)
        except (OSError, ValueError, KeyError):
            continue

    if not adata_list:
        raise RuntimeError("Failed to load any samples from extracted folders")

    # Make sample names unique
    if len(sample_names) != len(set(sample_names)):
        seen = {}
        unique_names = []
        for name in sample_names:
            if name in seen:
                seen[name] += 1
                unique_names.append(f"{name}_{seen[name]}")
            else:
                seen[name] = 0
                unique_names.append(name)
        sample_names = unique_names
        for adata, new_name in zip(adata_list, sample_names):
            adata.obs['sample'] = new_name

    # Combine
    if len(adata_list) > 1:
        adata_combined = sc.concat(adata_list, label='sample', keys=sample_names)
        adata_combined.obs_names_make_unique()
    else:
        adata_combined = adata_list[0]

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata_combined.write_h5ad(output_path)

    return adata_combined
