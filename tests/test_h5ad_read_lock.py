"""Stage 2b of ADR-014: Figure Export reads a study's file in a worker
thread while a step may be saving it.

Without coordination a save fails while a reader has the file open (on
Windows a replace fails too), and a reader that opens the file during a
save reads half of it. The readers and 'save_h5ad' share a per-file lock:
readers share it, a save has it alone, and a reader that finds the file
saved since it began stops with 'StudyFileChanged'.
"""
from __future__ import annotations

import os
import sys
import threading
import time

import anndata as ad
import h5py
import numpy as np
import pytest

from kosmic.scrna.load import h5ad_meta
from kosmic.scrna.load.h5ad_meta import (
    StudyFileChanged, file_lock, read_figure_meta, read_genes, save_h5ad,
)
from tests.study_folder import make_adata

SMALL_CHUNK = 1500


@pytest.fixture
def path(tmp_path):
    p = tmp_path / 's.h5ad'
    make_adata('counts').write_h5ad(p)
    return p


def test_without_the_lock_a_save_fails_while_a_reader_has_the_file_open(path):
    """The failure the lock exists for (HDF5 refuses to truncate an open file)."""
    with h5py.File(path, 'r'):
        with pytest.raises(OSError, match='already open'):
            make_adata('counts').write_h5ad(path)


@pytest.mark.skipif(sys.platform != 'win32', reason='replacing an open file is allowed elsewhere')
def test_without_the_lock_a_replace_fails_while_a_reader_has_the_file_open(path):
    tmp = path.with_name('s.h5ad.tmp')
    make_adata('counts').write_h5ad(tmp)
    with h5py.File(path, 'r'):
        with pytest.raises(PermissionError):
            os.replace(tmp, path)


@pytest.mark.parametrize('atomic', [False, True])
def test_a_save_during_a_chunked_read_waits_and_the_read_stops(path, monkeypatch, atomic):
    """The save waits for the chunk being read and succeeds; the reader's
    next chunk sees a new file and stops instead of mixing two versions."""
    first_chunk_read = threading.Event()
    save_done = threading.Event()
    real_block = h5ad_meta._Matrix._csr_block

    def _block(self, r0, r1):
        out = real_block(self, r0, r1)
        first_chunk_read.set()
        save_done.wait(timeout=10)             # the next chunk only after the save
        return out

    monkeypatch.setattr(h5ad_meta._Matrix, '_csr_block', _block)
    outcome = {}

    def _read():
        try:
            read_genes(path, ['VWF', 'G1'], source='log', chunk_bytes=SMALL_CHUNK)
            outcome['read'] = 'complete'
        except StudyFileChanged:
            outcome['read'] = 'changed'

    reader = threading.Thread(target=_read)
    reader.start()
    assert first_chunk_read.wait(timeout=10)
    new = make_adata('counts', seed=1)
    save_h5ad(new, path, atomic=atomic)            # no error, reader or not
    save_done.set()
    reader.join(timeout=10)
    assert outcome['read'] == 'changed'
    back = ad.read_h5ad(path)
    np.testing.assert_array_equal(back.X.toarray(), new.X.toarray())


def test_a_save_waits_for_a_reader_and_readers_share(path):
    lock = file_lock(path)
    events = []

    def _other_reader():
        with lock.reading():
            events.append('read')

    def _save():
        save_h5ad(make_adata('counts', seed=2), path)
        events.append('saved')

    with lock.reading():
        other_reader = threading.Thread(target=_other_reader)
        other_reader.start()
        other_reader.join(timeout=5)
        assert events == ['read']                  # readers do not wait for readers
        saver = threading.Thread(target=_save)
        saver.start()
        time.sleep(0.3)
        assert events == ['read']                  # the save waits for the reader
    saver.join(timeout=10)
    assert events == ['read', 'saved']


def test_a_waiting_save_goes_before_later_readers(path):
    lock = file_lock(path)
    order = []
    first = lock.reading()
    first.__enter__()
    saver = threading.Thread(
        target=lambda: (save_h5ad(make_adata('counts', seed=3), path), order.append('save')))
    saver.start()
    time.sleep(0.2)                                 # the save is now waiting

    def _late_reader():
        with lock.reading():
            order.append('late read')

    late = threading.Thread(target=_late_reader)
    late.start()
    time.sleep(0.2)
    assert order == []
    first.__exit__(None, None, None)
    saver.join(timeout=10)
    late.join(timeout=10)
    assert order == ['save', 'late read']


def test_the_lock_is_one_per_file_whatever_the_spelling(path):
    assert file_lock(path) is file_lock(str(path))
    assert file_lock(path) is file_lock(path.parent / '.' / path.name)
    if sys.platform == 'win32':
        assert file_lock(path) is file_lock(str(path).upper())


def test_an_unreadable_file_raises_rather_than_returning_partial_data(tmp_path, path):
    cut = tmp_path / 'cut.h5ad'
    cut.write_bytes(path.read_bytes()[: path.stat().st_size // 3])
    with pytest.raises(OSError):
        read_figure_meta(cut)


def test_a_file_without_an_embedding_reads_without_it(tmp_path):
    a = make_adata('counts')
    del a.obsm['X_umap']
    p = tmp_path / 'no_umap.h5ad'
    a.write_h5ad(p)
    assert 'X_umap' not in read_figure_meta(p).obsm
