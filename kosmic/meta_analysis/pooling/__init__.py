# Gene-level meta-analysis pooling methods.
#
# Each method has its own file with a single public '<method>_fast'
# entry point. Shared matrix helpers live in '.pooling_core' and are
# imported by every per-method file via relative imports.
#
# Re-exports below let callers write::
#
#     from kosmic.meta_analysis.pooling import dl_fast, reml_fast, sumrank_fast
#
# instead of spelling out the sub-module path. Both forms work.
from kosmic.meta_analysis.pooling.dl import dl_fast
from kosmic.meta_analysis.pooling.reml import reml_fast
from kosmic.meta_analysis.pooling.sumrank import sumrank_fast
from kosmic.meta_analysis.pooling.rop import rop_fast
from kosmic.meta_analysis.pooling.fisher import fisher_fast
from kosmic.meta_analysis.pooling.stouffer import stouffer_fast
from kosmic.meta_analysis.pooling.wop import wop_fast
from kosmic.meta_analysis.pooling.gwop import gwop_fast
from kosmic.meta_analysis.pooling.sign_test import pool_sign_test

__all__ = [
    'dl_fast',
    'reml_fast',
    'sumrank_fast',
    'rop_fast',
    'fisher_fast',
    'stouffer_fast',
    'wop_fast',
    'gwop_fast',
    'pool_sign_test',
]
