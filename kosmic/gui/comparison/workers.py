# Thread wrappers for the Methods Comparison workspace (ADR-017). The work
# is in kosmic.comparison.run; these run it off the GUI thread.

from kosmic.gui.shared.widgets import BaseWorker


class CheckWorker(BaseWorker):
    """Resolve the merged design and find the saved merged run.
    Emits 'finished_ok' with a ``kosmic.comparison.run.Check``."""

    def __init__(self, project, selection):
        super().__init__()
        self.project = project
        self.selection = selection

    def _run(self):
        from kosmic.comparison.run import check
        return check(self.project, self.selection)


class ComparisonWorker(BaseWorker):
    """Run and save a comparison. Emits 'finished_ok' with the
    ``ComparisonBundle``."""

    def __init__(self, project, selection, methods, rule, merged_action, checked):
        super().__init__()
        self.project = project
        self.selection = selection
        self.methods = list(methods)
        self.rule = rule
        self.merged_action = merged_action
        self.checked = checked

    def _run(self):
        from kosmic.comparison.run import run_comparison
        return run_comparison(self.project, self.selection, self.methods, self.rule,
                              merged_action=self.merged_action, checked=self.checked,
                              progress=self.progress.emit)
