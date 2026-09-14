# Shared GUI primitives used by every workspace.
#
# Layout:
#   shared/theme.py           - palette + stylesheet + fonts + NoScroll* widgets
#   shared/icon_provider.py   - SVG icon loader with LRU cache
#   shared/widgets/           - CopyableTableWidget + future shared QWidgets
#   shared/plots/             - interactive pyqtgraph plot widgets
#   shared/tutorial/          - guided tutorial state machine + overlay + data
#   shared/run_worker.py      - uniform helper for wiring + starting a BaseWorker

from kosmic.gui.shared.borderless import borderless
from kosmic.gui.shared.run_worker import run_worker

__all__ = ["borderless", "run_worker"]
