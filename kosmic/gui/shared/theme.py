# Centralised dark / light theming via QPalette + a generated global QSS.
from enum import Enum
from PyQt6.QtCore import QSettings
from PyQt6.QtGui import QPalette, QColor
from PyQt6.QtWidgets import (
    QApplication, QComboBox, QDoubleSpinBox, QSizePolicy, QSpinBox,
)
import re

from kosmic import DEFAULT_BASE_FONT_SIZE, DEFAULT_EXPORT_DPI, UI_FONT_SCALE


class ThemeMode(Enum):
    DARK = "dark"
    LIGHT = "light"


# Color palettes

DARK_COLORS = {
    # Backgrounds
    'bg_primary': '#1e1e1e',
    'bg_secondary': '#252526',
    'bg_tertiary': '#2d2d30',
    'bg_hover': '#2a2d2e',
    'bg_active': '#37373d',
    # Recessed, not raised: at #3c3c3c inputs read as grey slabs
    # against the #252526 chrome, 20% lighter than every other dark
    # surface. Sitting at bg_primary lets the border do the work.
    'bg_input': '#1e1e1e',
    'bg_card': '#2d2d30',

    # Foregrounds
    'fg_primary': '#cccccc',
    'fg_secondary': '#858585',
    'fg_tertiary': '#6a6a6a',
    'fg_inverse': '#ffffff',

    # Accents
    'accent_primary': '#007acc',
    'accent_hover': '#1c8ad6',
    'accent_active': '#0e639c',

    # Semantic
    'success': '#4ec9b0',
    'warning': '#ce9178',
    'error': '#f48771',
    'info': '#4fc1ff',

    # Borders
    'border': '#3e3e42',
    'border_focus': '#007acc',
    'divider': '#2b2b2b',

    # Plot colors (semantic, used across all charts)
    'plot_upregulated': '#e05252',
    'plot_downregulated': '#5288e0',
    'plot_insignificant': '#666666',
    'plot_control': '#4A9BD9',
    'plot_disease': '#D95A5A',
    'plot_bg': '#1e1e1e',
    'plot_fg': '#cccccc',
    'plot_dot': 'rgba(255, 255, 255, 180)',
}

LIGHT_COLORS = {
    # Backgrounds  (aligned with Adobe Illustrator lightest)
    'bg_primary': '#ffffff',
    'bg_secondary': '#f0f0f0',
    'bg_tertiary': '#e4e4e4',
    'bg_hover': '#e0e0e0',
    'bg_active': '#d4d4d4',
    'bg_input': '#ffffff',
    'bg_card': '#ffffff',

    # Foregrounds
    'fg_primary': '#333333',
    'fg_secondary': '#505050',
    'fg_tertiary': '#6e6e6e',
    'fg_inverse': '#ffffff',

    # Accents
    'accent_primary': '#0078d4',
    'accent_hover': '#106ebe',
    'accent_active': '#005a9e',

    # Semantic
    'success': '#107c10',
    'warning': '#ca5010',
    'error': '#d13438',
    'info': '#0078d4',

    # Borders
    'border': '#c9c9c9',
    'border_focus': '#0078d4',
    'divider': '#c9c9c9',

    # Plot colors (semantic, used across all charts)
    'plot_upregulated': '#cc3333',
    'plot_downregulated': '#3366cc',
    'plot_insignificant': '#888888',
    'plot_control': '#3498DB',
    'plot_disease': '#E74C3C',
    'plot_bg': '#ffffff',
    'plot_fg': '#333333',
    'plot_dot': 'rgba(80, 80, 80, 150)',
}

_current_mode: ThemeMode = ThemeMode.DARK


def get_current_mode() -> ThemeMode:
    """Return the current theme mode."""
    return _current_mode


def get_color(key: str) -> str:
    """Return a single colour value for the current theme."""
    colors = DARK_COLORS if _current_mode == ThemeMode.DARK else LIGHT_COLORS
    return colors.get(key, '#000000')


def get_colors() -> dict:
    """Return the full colour dictionary for the current theme."""
    return DARK_COLORS if _current_mode == ThemeMode.DARK else LIGHT_COLORS


def style_pg_plot(plot_widget, title: str = '', left_label: str = '',
                  bottom_label: str = ''):
    """Apply the current theme to a pyqtgraph PlotWidget.

    Call this after creating a pg.PlotWidget to match the app theme.
    Reads font sizes from plot settings so the user's base font is respected.
    All labels are linked to View → Plot Settings base font size.
    """
    from PyQt6.QtGui import QFont

    bg = get_color('bg_primary')
    fg = get_color('fg_primary')
    fg2 = get_color('fg_secondary')
    border = get_color('border')
    fs = get_font_sizes()

    plot_widget.setBackground(bg)

    # Axis styling with themed tick font
    tick_font = QFont()
    tick_font.setPointSizeF(fs['tick'])
    for axis_name in ('bottom', 'left', 'top', 'right'):
        ax = plot_widget.getAxis(axis_name)
        ax.setPen(border)
        ax.setTextPen(fg2)
        ax.setTickFont(tick_font)

    title_size = f'{fs["title"]:.0f}pt'
    label_size = f'{fs["axis_label"]:.0f}pt'
    if title:
        plot_widget.setTitle(title, color=fg, size=title_size)
    if left_label:
        plot_widget.setLabel('left', left_label, color=fg2, size=label_size)
    if bottom_label:
        plot_widget.setLabel('bottom', bottom_label, color=fg2, size=label_size)

    # Title colour and size for future setTitle calls
    plot_widget.getPlotItem().titleLabel.opts['color'] = fg
    plot_widget.getPlotItem().titleLabel.opts['size'] = title_size
    # Store label size for future setLabel calls
    plot_widget._kosmic_label_size = label_size


def style_mpl_plot(fig, ax, title: str = '', left_label: str = '',
                   bottom_label: str = ''):
    """Apply the current theme to a matplotlib Figure + Axes pair.

    Mirror of 'style_pg_plot' for matplotlib-backed plots. Sets
    background, spine and tick colours, and label fonts from the user's
    plot-settings base font size. Pass either a single ax or an iterable.
    """
    bg = get_color('bg_primary')
    fg = get_color('fg_primary')
    fg2 = get_color('fg_secondary')
    border = get_color('border')
    fs = get_font_sizes()

    fig.patch.set_facecolor(bg)

    axes = ax if hasattr(ax, '__iter__') else [ax]
    for axis in axes:
        axis.set_facecolor(bg)
        for spine in axis.spines.values():
            spine.set_edgecolor(border)
        axis.tick_params(colors=fg2, labelsize=fs['tick'])
        axis.xaxis.label.set_color(fg2)
        axis.yaxis.label.set_color(fg2)
        if title:
            axis.set_title(title, color=fg, fontsize=fs['title'])
        if left_label:
            axis.set_ylabel(left_label, color=fg2, fontsize=fs['axis_label'])
        if bottom_label:
            axis.set_xlabel(bottom_label, color=fg2,
                            fontsize=fs['axis_label'])


def get_font_sizes() -> dict:
    """Return a dict of named font sizes based on the user's base font setting.

    Keys: title, suptitle, axis_label, tick, legend, annotation,
    colorbar_label, colorbar_tick, base.  Uses the same ratios as
    'kosmic.visualisation.default_font_sizes()'.
    """
    from kosmic.visualisation import default_font_sizes
    base = int(QSettings("KOSMIC", "KOSMIC").value(
        "plot/base_font_size", DEFAULT_BASE_FONT_SIZE))
    return default_font_sizes(base)


def get_plot_settings() -> dict:
    """Read global plot settings from QSettings.

    Returns dict with all plot configuration including font_sizes,
    export_dpi, and per-plot-type settings (heatmap, volcano, bar, dotplot).
    """
    s = QSettings("KOSMIC", "KOSMIC")
    return {
        'export_dpi': int(s.value("plot/export_dpi", DEFAULT_EXPORT_DPI)),
        'font_sizes': get_font_sizes(),
        # Heatmap
        'heatmap_cmap': s.value("plot/heatmap_cmap", "viridis"),
        'heatmap_vmin': float(s.value("plot/heatmap_vmin", 5)),
        'heatmap_vmax': float(s.value("plot/heatmap_vmax", 95)),
        # Volcano
        'volcano_color_scheme': int(s.value("plot/volcano_color_scheme", 0)),
        'volcano_point_size': int(s.value("plot/volcano_point_size", 9)),
        'volcano_max_labels': int(s.value("plot/volcano_max_labels", 10)),
    }


def get_export_dpi() -> int:
    """Return the user-configured export DPI (default 150)."""
    return int(QSettings("KOSMIC", "KOSMIC").value("plot/export_dpi", DEFAULT_EXPORT_DPI))




def get_palette() -> QPalette:
    """Build a QPalette for the current theme."""
    c = get_colors()
    palette = QPalette()

    palette.setColor(QPalette.ColorRole.Window, QColor(c['bg_secondary']))
    palette.setColor(QPalette.ColorRole.WindowText, QColor(c['fg_primary']))
    palette.setColor(QPalette.ColorRole.Base, QColor(c['bg_input']))
    palette.setColor(QPalette.ColorRole.AlternateBase, QColor(c['bg_tertiary']))
    palette.setColor(QPalette.ColorRole.ToolTipBase, QColor(c['bg_primary']))
    palette.setColor(QPalette.ColorRole.ToolTipText, QColor(c['fg_primary']))
    palette.setColor(QPalette.ColorRole.Text, QColor(c['fg_primary']))
    palette.setColor(QPalette.ColorRole.Button, QColor(c['bg_tertiary']))
    palette.setColor(QPalette.ColorRole.ButtonText, QColor(c['fg_primary']))
    palette.setColor(QPalette.ColorRole.BrightText, QColor(c['error']))
    palette.setColor(QPalette.ColorRole.Link, QColor(c['accent_primary']))
    palette.setColor(QPalette.ColorRole.Highlight, QColor(c['accent_primary']))
    palette.setColor(QPalette.ColorRole.HighlightedText, QColor(c['fg_inverse']))

    palette.setColor(QPalette.ColorGroup.Disabled,
                     QPalette.ColorRole.WindowText, QColor(c['fg_tertiary']))
    palette.setColor(QPalette.ColorGroup.Disabled,
                     QPalette.ColorRole.Text, QColor(c['fg_tertiary']))
    palette.setColor(QPalette.ColorGroup.Disabled,
                     QPalette.ColorRole.ButtonText, QColor(c['fg_tertiary']))

    return palette


#: Explicit type ladder: source size -> shipped size, in px.
#:
#: The stylesheet grew ten distinct sizes (8, 10, 11, 12, 13, 14, 15, 16,
#: 21, 22) doing overlapping jobs -- 11px and 12px were both "small", 12px
#: and 13px were both "body", 14/15/16 were all "heading". Multiplying them
#: by a scale factor preserved that muddle and rounding made it worse, so
#: sizes one pixel apart ended up two or three apart.
#:
#: This collapses them onto six deliberate steps. Even 2px gaps at the
#: bottom where the eye notices small differences, wider jumps at the top.
_FONT_LADDER_PX = {
    8: 12, 10: 12, 11: 12,      # caption, secondary, footnotes
    12: 14, 13: 14,             # body -- the bulk of the app
    14: 16, 15: 16,             # emphasis, section + inline headings
    16: 18,                     # page headers
    21: 22, 22: 24,             # page titles, stat-tile numbers
}
#: No pt sizes remain in the sheet; kept so a stray one is not
#: silently left off the ladder.
_FONT_LADDER_PT: dict[int, int] = {}


def normalise_font_sizes(qss: str, scale: float = None) -> str:
    """Map every ``font-size`` in a stylesheet onto :data:`_FONT_LADDER_PX`.

    ``scale`` (default :data:`UI_FONT_SCALE`) is applied afterwards, so the
    ladder fixes the *relationships* and the config knob moves the whole
    thing up or down together. A size not on the ladder is scaled but
    otherwise left alone, so a new rule cannot silently vanish.

    Sizes never fall below 8px / 6pt.
    """
    if scale is None:
        scale = UI_FONT_SCALE

    def _remap(match: re.Match) -> str:
        value, unit = float(match.group(1)), match.group(2)
        ladder = _FONT_LADDER_PX if unit == 'px' else _FONT_LADDER_PT
        mapped = ladder.get(int(value), value) if value.is_integer() else value
        floor = 8.0 if unit == 'px' else 6.0
        return f"font-size: {max(round(mapped * scale), floor):g}{unit}"

    return re.sub(r"font-size:\s*([0-9.]+)(px|pt)", _remap, qss)


def get_stylesheet() -> str:
    """Generate the global QSS stylesheet for the current theme."""
    c = get_colors()

    return normalise_font_sizes(f"""
    /* === Base ============================================= */

    QWidget {{
        background-color: {c['bg_secondary']};
        color: {c['fg_primary']};
        font-size: 13px;
    }}

    /* === Scroll Areas ===================================== */

    QScrollArea {{
        background: transparent;
        border: none;
    }}

    QScrollArea > QWidget > QWidget {{
        background: transparent;
    }}

    /* === Group Boxes ====================================== */

    QGroupBox {{
        background-color: {c['bg_secondary']};
        border: 1px solid {c['border']};
        border-radius: 6px;
        margin-top: 8px;
        padding: 8px 8px 6px 8px;
        font-weight: 600;
    }}

    QGroupBox::title {{
        subcontrol-origin: margin;
        subcontrol-position: top left;
        left: 12px;
        top: 0px;
        padding: 0 6px;
        background-color: {c['bg_secondary']};
        color: {c['fg_primary']};
        font-size: 13px;
        font-weight: 600;
    }}

    /* === Buttons ========================================== */

    QPushButton {{
        background-color: {c['bg_tertiary']};
        color: {c['fg_primary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        padding: 3px 10px;
        font-weight: 500;
    }}

    QPushButton:hover {{
        background-color: {c['bg_active']};
        border-color: {c['accent_primary']};
    }}

    QPushButton:pressed {{
        background-color: {c['bg_active']};
    }}

    QPushButton:disabled {{
        background-color: {c['bg_secondary']};
        color: {c['fg_tertiary']};
        border-color: {c['border']};
    }}

    QPushButton:checked {{
        background-color: {c['accent_primary']};
        color: {c['fg_inverse']};
        border-color: {c['accent_primary']};
    }}

    /* Primary action buttons */
    QPushButton#primary_button {{
        background-color: {c['accent_primary']};
        color: {c['fg_inverse']};
        border-color: {c['accent_primary']};
        font-weight: 600;
    }}

    QPushButton#primary_button:hover {{
        background-color: {c['accent_hover']};
        border-color: {c['accent_hover']};
    }}

    QPushButton#primary_button:pressed {{
        background-color: {c['accent_active']};
        border-color: {c['accent_active']};
    }}

    /* The ID selector above outranks the generic
       'QPushButton:disabled', so without this a disabled primary
       button keeps its accent colour and reads as clickable. */
    QPushButton#primary_button:disabled {{
        background-color: {c['bg_secondary']};
        color: {c['fg_tertiary']};
        border-color: {c['border']};
    }}

    /* === Inputs =========================================== */

    QLineEdit, QTextEdit, QPlainTextEdit {{
        background-color: {c['bg_input']};
        color: {c['fg_primary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        padding: 3px 8px;
        selection-background-color: {c['accent_primary']};
    }}

    QLineEdit:focus, QTextEdit:focus, QPlainTextEdit:focus {{
        border: 1px solid {c['border_focus']};
    }}

    QLineEdit:disabled, QTextEdit:disabled, QPlainTextEdit:disabled {{
        background-color: {c['bg_secondary']};
        color: {c['fg_tertiary']};
    }}

    QSpinBox, QDoubleSpinBox {{
        padding: 3px 4px;
    }}

    /* === ComboBox ========================================= */

    QComboBox {{
        background-color: {c['bg_input']};
        color: {c['fg_primary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        padding: 5px 8px;
        min-height: 20px;
    }}

    QComboBox:hover {{
        border-color: {c['accent_primary']};
    }}

    QComboBox:focus {{
        border-color: {c['border_focus']};
    }}

    QComboBox QAbstractItemView {{
        background-color: {c['bg_input']};
        color: {c['fg_primary']};
        border: 1px solid {c['border']};
        selection-background-color: {c['accent_primary']};
        selection-color: {c['fg_inverse']};
        outline: none;
    }}

    /* === Lists & Tables ================================== */

    QListWidget, QTableWidget {{
        background-color: {c['bg_input']};
        color: {c['fg_primary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        alternate-background-color: {c['bg_secondary']};
        selection-background-color: {c['accent_primary']};
        selection-color: {c['fg_inverse']};
        outline: none;
    }}

    QListWidget::item, QTableWidget::item {{
        padding: 4px;
    }}

    QListWidget::item:hover, QTableWidget::item:hover {{
        background-color: {c['bg_hover']};
    }}

    QListWidget::item:selected, QTableWidget::item:selected {{
        background-color: {c['accent_primary']};
        color: {c['fg_inverse']};
    }}

    QHeaderView::section {{
        background-color: {c['bg_tertiary']};
        color: {c['fg_primary']};
        border: none;
        border-right: 1px solid {c['border']};
        border-bottom: 1px solid {c['border']};
        padding: 6px 8px;
        font-weight: 600;
    }}

    /* === Progress Bar ===================================== */

    QProgressBar {{
        background-color: {c['bg_secondary']};
        border: 1px solid {c['border']};
        border-radius: 3px;
        text-align: center;
        color: {c['fg_primary']};
        font-size: 11px;
        max-height: 18px;
    }}

    QProgressBar::chunk {{
        background-color: {c['accent_primary']};
        border-radius: 2px;
    }}

    /* === Scrollbars ======================================= */

    QScrollBar:vertical {{
        background: {c['bg_primary']};
        width: 12px;
        margin: 0;
        border-radius: 6px;
    }}

    QScrollBar::handle:vertical {{
        background: {c['bg_active']};
        min-height: 30px;
        border-radius: 6px;
    }}

    QScrollBar::handle:vertical:hover {{
        background: {c['fg_tertiary']};
    }}

    QScrollBar:horizontal {{
        background: {c['bg_primary']};
        height: 12px;
        margin: 0;
        border-radius: 6px;
    }}

    QScrollBar::handle:horizontal {{
        background: {c['bg_active']};
        min-width: 30px;
        border-radius: 6px;
    }}

    QScrollBar::handle:horizontal:hover {{
        background: {c['fg_tertiary']};
    }}

    QScrollBar::add-line, QScrollBar::sub-line,
    QScrollBar::add-page, QScrollBar::sub-page {{
        background: none;
        border: none;
        height: 0px;
        width: 0px;
    }}

    /* === Splitter ========================================= */

    QSplitter::handle {{
        background-color: {c['divider']};
    }}

    QSplitter::handle:horizontal {{
        width: 1px;
    }}

    QSplitter::handle:vertical {{
        height: 1px;
    }}

    /* === Checkboxes & Radios ============================== */

    QCheckBox, QRadioButton {{
        spacing: 6px;
        color: {c['fg_primary']};
    }}

    QCheckBox::indicator, QRadioButton::indicator {{
        width: 16px;
        height: 16px;
        border: 1px solid {c['border']};
        background-color: {c['bg_input']};
    }}

    QCheckBox::indicator {{
        border-radius: 3px;
    }}

    QRadioButton::indicator {{
        border-radius: 8px;
    }}

    QCheckBox::indicator:checked, QRadioButton::indicator:checked {{
        background-color: {c['accent_primary']};
        border-color: {c['accent_primary']};
    }}

    QCheckBox::indicator:hover, QRadioButton::indicator:hover {{
        border-color: {c['accent_primary']};
    }}

    /* === Tab Widget (internal tabs, e.g. DE Analysis) ===== */

    QTabWidget::pane {{
        border: 1px solid {c['border']};
        background-color: {c['bg_primary']};
        border-radius: 4px;
    }}

    QTabBar::tab {{
        background-color: {c['bg_tertiary']};
        color: {c['fg_secondary']};
        border: 1px solid {c['border']};
        border-bottom: none;
        border-top-left-radius: 4px;
        border-top-right-radius: 4px;
        padding: 6px 16px;
        margin-right: 2px;
    }}

    QTabBar::tab:selected {{
        background-color: {c['bg_primary']};
        color: {c['fg_primary']};
        font-weight: 600;
        border-bottom: 1px solid {c['bg_primary']};
    }}

    QTabBar::tab:hover:!selected {{
        background-color: {c['bg_hover']};
    }}

    /* === Menu Bar ========================================= */

    QMenuBar {{
        background-color: {c['bg_secondary']};
        color: {c['fg_primary']};
        border-bottom: 1px solid {c['border']};
        padding: 1px 0;
    }}

    QMenuBar::item {{
        padding: 4px 8px;
        background: transparent;
    }}

    QMenuBar::item:selected {{
        background-color: {c['bg_hover']};
    }}

    /* === Status Bar ======================================= */

    QStatusBar {{
        background-color: {c['bg_secondary']};
        color: {c['fg_secondary']};
        border-top: 1px solid {c['border']};
        font-size: 12px;
        padding: 2px 8px;
    }}

    /* === Tooltips ========================================= */

    QToolTip {{
        background-color: {c['bg_tertiary']};
        color: {c['fg_primary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        padding: 4px 8px;
    }}

    /* === Semantic labels ================================== */

    QLabel[role="header"] {{
        font-size: 15px;
        font-weight: 600;
        color: {c['fg_primary']};
        background-color: transparent;
    }}

    QLabel[role="secondary"] {{
        font-size: 11px;
        color: {c['fg_secondary']};
        background-color: transparent;
    }}

    /* === Section headers (uppercase sidebar labels) =============== */

    QLabel[role="section_header"] {{
        color: {c['fg_secondary']};
        font-family: "Segoe UI";
        font-size: 11px;
        font-weight: 600;
        padding-top: 4px;
        background-color: transparent;
    }}

    /* === Cards (bordered content panels) ========================= */

    QFrame[role="card"] {{
        background-color: {c['bg_secondary']};
        border: 2px solid {c['border']};
        border-radius: 8px;
        padding: 20px;
    }}

    /* Child labels inside cards inherit no border so inline
       'border: none' overrides are unnecessary. */
    QFrame[role="card"] QLabel {{
        background-color: transparent;
        border: none;
    }}

    /* Compact bordered panel (smaller than card).  Used for inline
       config sections / dataset lists / loaded-state boxes.  Padding
       is intentionally modest; callers may override via contentsMargins. */
    QFrame[role="panel"] {{
        background-color: {c['bg_secondary']};
        border: 1px solid {c['border']};
        border-radius: 6px;
    }}

    QFrame[role="panel"] QLabel {{
        background-color: transparent;
        border: none;
    }}

    /* === 1px horizontal dividers ================================== */

    QWidget[role="divider"] {{
        background-color: {c['border']};
        border: none;
    }}

    /* === File tree (main.py explorer pane) ======================= */

    QTreeView {{
        background: {c['bg_secondary']};
        color: {c['fg_primary']};
        border: none;
        selection-background-color: {c['bg_active']};
        selection-color: {c['fg_primary']};
    }}

    QTreeView::item {{
        border: none;
    }}

    QTreeView::item:hover {{
        background: {c['bg_hover']};
        color: {c['fg_primary']};
    }}

    QTreeView::item:selected,
    QTreeView::item:selected:active {{
        background: {c['bg_active']};
        color: {c['fg_primary']};
    }}

    QTreeView::branch {{
        background: {c['bg_secondary']};
    }}

    QTreeView::branch:selected {{
        background: {c['bg_active']};
    }}

    /* === Hint labels (small italic secondary text) =============== */

    QLabel[role="hint"] {{
        color: {c['fg_secondary']};
        font-size: 10px;
        font-style: italic;
        background-color: transparent;
    }}

    /* === Caption / value labels (KPI grids, field/value pairs) ===== */

    QLabel[role="caption"] {{
        color: {c['fg_secondary']};
        font-size: 11px;
        background-color: transparent;
        border: none;
    }}

    QLabel[role="value"] {{
        color: {c['fg_primary']};
        font-size: 12px;
        font-weight: bold;
        background-color: transparent;
        border: none;
    }}

    /* Status-bar permanent metric (CPU/RAM monitor in main window). */
    QLabel[role="status_metric"] {{
        padding: 0 10px;
    }}

    /* === Status labels (state-coloured inline text) ================ */

    QLabel[role="status_info"] {{
        color: {c['info']};
        font-size: 12px;
        background-color: transparent;
    }}

    QLabel[role="status_success"] {{
        color: {c['success']};
        font-size: 12px;
        background-color: transparent;
    }}

    QLabel[role="status_warning"] {{
        color: {c['warning']};
        font-size: 12px;
        background-color: transparent;
    }}

    QLabel[role="status_error"] {{
        color: {c['error']};
        font-size: 12px;
        background-color: transparent;
    }}

    /* === Secondary / danger / tab buttons ========================== */

    QPushButton#secondary_button {{
        background-color: {c['bg_tertiary']};
        color: {c['fg_primary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        padding: 4px 20px;
        font-size: 12px;
    }}

    QPushButton#secondary_button:hover {{
        background-color: {c['bg_active']};
        border-color: {c['fg_tertiary']};
    }}

    QPushButton#danger_button {{
        background-color: {c['error']};
        color: {c['fg_inverse']};
        border: 1px solid {c['error']};
        border-radius: 4px;
        padding: 4px 24px;
        font-size: 13px;
        font-weight: bold;
    }}

    QPushButton#danger_button:hover {{
        background-color: {c['warning']};
        border-color: {c['warning']};
    }}

    QPushButton#tab_button {{
        background-color: transparent;
        color: {c['fg_secondary']};
        border: none;
        border-bottom: 2px solid transparent;
        font-size: 12px;
        padding: 6px 16px;
    }}

    QPushButton#tab_button:hover {{
        color: {c['fg_primary']};
    }}

    QPushButton#tab_button[active="true"] {{
        color: {c['fg_primary']};
        border-bottom: 2px solid {c['accent_primary']};
        font-weight: bold;
    }}

    /* === Ghost button (transparent bg + thin border + hover) ========= */

    QPushButton#ghost_button {{
        background-color: {c['bg_primary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        font-size: 12px;
        padding: 6px 12px;
    }}

    QPushButton#ghost_button:hover {{
        background-color: {c['bg_hover']};
    }}

    /* === Flat text button (tiny transparent "link"-style buttons) === */

    QPushButton[role="flat_text"] {{
        background-color: transparent;
        color: {c['fg_secondary']};
        border: none;
        font-size: 11px;
    }}

    QPushButton[role="flat_text"]:hover {{
        color: {c['fg_primary']};
    }}

    /* === Log output text area (OutputPanel) ===================== */

    QTextEdit[role="log_output"] {{
        background-color: {c['bg_primary']};
        color: {c['fg_primary']};
        padding: 4px 6px;
        selection-background-color: {c['bg_hover']};
    }}

    QFrame[role="workflow_card"] {{
        background-color: {c['bg_secondary']};
        border: 1px solid {c['border']};
        border-radius: 6px;
    }}
    QLabel[role="workflow_card_title"] {{
        color: {c['fg_primary']};
    }}

    QFrame[role="walkthrough_card"] {{
        background-color: {c['bg_secondary']};
        border: 1px solid {c['border']};
        border-radius: 6px;
    }}
    QToolButton[role="walkthrough_header"] {{
        background-color: transparent;
        border: none;
        text-align: left;
        padding: 4px 12px;
        color: {c['fg_primary']};
        font-weight: 600;
    }}
    QToolButton[role="walkthrough_header"]:hover {{
        background-color: {c['bg_hover']};
    }}
    QToolButton[role="walkthrough_header"]:checked {{
        background-color: {c['bg_tertiary']};
        border-bottom: 1px solid {c['border']};
    }}

    QLabel[role="active_badge"] {{
        background-color: {c['accent_primary']};
        color: white;
        font-size: 10px;
        font-weight: 600;
        padding: 2px 8px;
        border-radius: 8px;
    }}

    /* Flat in-table action buttons -- transparent so the row colour
       shows through; only background changes on hover. */
    QPushButton[role="row_action"] {{
        background-color: transparent;
        color: {c['accent_primary']};
        border: none;
        padding: 2px 10px;
        font-size: 11px;
        font-weight: 500;
    }}
    QPushButton[role="row_action"]:hover {{
        background-color: {c['bg_hover']};
        border-radius: 3px;
    }}

    /* Hamburger-icon menu trigger inside table rows. */
    QToolButton[role="row_kebab"] {{
        background-color: transparent;
        border: none;
        padding: 0px;
    }}
    QToolButton[role="row_kebab"]:hover {{
        background-color: {c['bg_hover']};
        border-radius: 3px;
    }}
    QToolButton[role="row_kebab"]::menu-indicator {{
        image: none;
    }}

    /* === Compact list-item button used by setup/download dataset cards
           -- same shape but slightly roomier padding, different active state */

    QPushButton[role="dataset_card"] {{
        background-color: {c['bg_primary']};
        color: {c['fg_secondary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        padding: 8px 10px;
        font-size: 13px;
        text-align: left;
    }}

    QPushButton[role="dataset_card"]:hover {{
        color: {c['fg_primary']};
        border-color: {c['accent_primary']};
    }}

    QPushButton[role="dataset_card"][active="true"] {{
        color: {c['fg_primary']};
        border: 2px solid {c['success']};
    }}

    /* === Overview primitives (widgets/overview.py) =================== */

    /* Hero card for the active dataset (scRNA Dataset screen) */
    QFrame#dataset_hero_card {{
        background-color: {c['bg_card']};
        border: 1px solid {c['border']};
        border-radius: 8px;
    }}

    /* Card children paint no background of their own (the base QWidget
       rule would draw bg_secondary boxes on the lighter card). */
    QFrame#dataset_hero_card QLabel,
    QFrame#info_panel QLabel,
    StatBlock {{
        background-color: transparent;
    }}

    /* Re-asserted after the transparent sweep above (higher-specificity
       descendant selectors would otherwise win). */
    QFrame#dataset_hero_card QLabel[role="active_badge"] {{
        background-color: {c['accent_primary']};
    }}

    /* IconTile: rounded file-type / source tile */
    QFrame#icon_tile {{
        background-color: {c['bg_tertiary']};
        border: 1px solid {c['border']};
        border-radius: 10px;
    }}

    QFrame#icon_tile QLabel {{
        background-color: transparent;
        border: none;
    }}

    QLabel[role="icon_tile_caption"] {{
        color: {c['accent_primary']};
        font-size: 8px;
        font-weight: 700;
        letter-spacing: 0.5px;
    }}

    QLabel[role="icon_tile_text"] {{
        color: {c['accent_primary']};
        font-size: 13px;
        font-weight: 700;
        letter-spacing: 1px;
    }}

    /* Selectable card row (GEO file list, assemble-from-parts rows) */
    QFrame[role="file_row"] {{
        background-color: {c['bg_card']};
        border: 1px solid {c['border']};
        border-radius: 8px;
    }}

    QFrame[role="file_row"]:hover {{
        border-color: {c['accent_hover']};
    }}

    QFrame[role="file_row"][selected="true"] {{
        border: 1px solid {c['accent_primary']};
    }}

    QFrame[role="file_row"] QLabel {{
        background-color: transparent;
    }}

    QLabel[role="file_row_name"] {{
        font-size: 13px;
        font-weight: 600;
        color: {c['fg_primary']};
    }}

    /* Tiny neutral tag pill (Required / Optional) */
    QLabel[role="tag"] {{
        background-color: {c['bg_tertiary']};
        color: {c['fg_secondary']};
        border: 1px solid {c['border']};
        border-radius: 8px;
        padding: 1px 8px;
        font-size: 10px;
    }}

    /* === Stage-flow primitives (widgets/stage_flow.py) =============== */

    QFrame#stage_card {{
        background-color: {c['bg_card']};
        border: 1px solid {c['border']};
        border-radius: 8px;
    }}

    QFrame#stage_card:hover {{
        border-color: {c['accent_hover']};
    }}

    QFrame#stage_card QLabel {{
        background-color: transparent;
    }}

    /* Same reason as the QLabel rule above: once any stylesheet rule
       matches a widget, Qt paints its palette background rather than
       leaving it clear, so a bare checkbox or radio inside a card drew
       a bg_secondary band across the card's lighter face. Input widgets
       (combo, spin, line edit) are left alone -- they are meant to
       carry their own fill. */
    QFrame#stage_card QCheckBox,
    QFrame#stage_card QRadioButton {{
        background-color: transparent;
    }}

    /* Selected-gene readout floating over the volcano. Opaque, because
       scatter points run underneath it. */
    QLabel#volcano_gene_legend {{
        background-color: {c['bg_card']};
        border: 1px solid {c['border']};
        border-radius: 6px;
        padding: 6px 8px;
        color: {c['fg_primary']};
    }}

    QLabel[role="stage_card_title"] {{
        font-size: 14px;
        font-weight: 600;
        color: {c['fg_primary']};
    }}

    /* Thin horizontal rule under flat section headers */
    QFrame[role="section_rule"] {{
        background-color: {c['border']};
        border: none;
        max-height: 1px;
        min-height: 1px;
    }}

    /* Numbered step badge on section cards. The card-scoped selector
       re-asserts the fill over the info_panel transparent-QLabel rule. */
    QLabel[role="step_badge"],
    QFrame#info_panel QLabel[role="step_badge"] {{
        background-color: {c['accent_primary']};
        color: #ffffff;
        border-radius: 11px;
        font-size: 11px;
        font-weight: 600;
    }}

    /* Completion banners (Setup page) */
    QFrame[role="banner_success"] {{
        background-color: transparent;
        border: 1px solid {c['success']};
        border-radius: 8px;
    }}

    QFrame[role="banner_success"] QLabel {{
        color: {c['success']};
        background-color: transparent;
        font-size: 12px;
    }}

    QFrame[role="banner_info"] {{
        background-color: transparent;
        border: 1px solid {c['border']};
        border-radius: 8px;
    }}

    QFrame[role="banner_info"] QLabel {{
        color: {c['fg_secondary']};
        background-color: transparent;
        font-size: 12px;
    }}

    /* Clickable expander bar (View all metadata columns) */
    QFrame[role="expander_bar"] {{
        background-color: {c['bg_card']};
        border: 1px solid {c['border']};
        border-radius: 8px;
    }}

    QFrame[role="expander_bar"]:hover {{
        border-color: {c['accent_hover']};
    }}

    QFrame[role="expander_bar"] QLabel {{
        background-color: transparent;
    }}


    /* Quiet accent link button (Edit, Close, View elbow plot) */
    QPushButton[role="link_accent"] {{
        background-color: transparent;
        border: none;
        color: {c['accent_primary']};
        font-size: 12px;
        padding: 2px 4px;
    }}

    QPushButton[role="link_accent"]:hover {{
        color: {c['accent_hover']};
        text-decoration: underline;
    }}

    /* Dashed drag-and-drop target */
    QFrame[role="drop_zone"] {{
        background-color: transparent;
        border: 1px dashed {c['border']};
        border-radius: 8px;
    }}

    QFrame[role="drop_zone"] QLabel {{
        background-color: transparent;
    }}

    QLabel#hero_dataset_name {{
        font-size: 21px;
        font-weight: 600;
        color: {c['fg_primary']};
    }}

    /* StatBlock: big number over a small caption */
    QLabel[role="stat_value"] {{
        font-size: 22px;
        font-weight: 600;
        color: {c['fg_primary']};
    }}

    QFrame[role="stat_divider"] {{
        background-color: {c['border']};
    }}

    /* InfoPanel: titled key/value panel */
    QFrame#info_panel {{
        background-color: {c['bg_card']};
        border: 1px solid {c['border']};
        border-radius: 8px;
    }}

    QLabel[role="panel_title"] {{
        color: {c['fg_secondary']};
        font-size: 11px;
        font-weight: 600;
        letter-spacing: 1px;
    }}

    QLabel[role="kv_key"] {{
        color: {c['fg_secondary']};
        font-size: 12px;
    }}

    QLabel[role="kv_value"] {{
        color: {c['fg_primary']};
        font-size: 12px;
        font-weight: 500;
    }}

    /* StatusChip: rounded page-state pill */
    QLabel[role="chip_success"] {{
        background-color: {c['bg_tertiary']};
        color: {c['success']};
        border: 1px solid {c['border']};
        border-radius: 13px;
        padding: 5px 14px;
        font-size: 12px;
        font-weight: 600;
    }}

    QLabel[role="chip_info"] {{
        background-color: {c['bg_tertiary']};
        color: {c['info']};
        border: 1px solid {c['border']};
        border-radius: 13px;
        padding: 5px 14px;
        font-size: 12px;
        font-weight: 600;
    }}

    QLabel[role="chip_warning"] {{
        background-color: {c['bg_tertiary']};
        color: {c['warning']};
        border: 1px solid {c['border']};
        border-radius: 13px;
        padding: 5px 14px;
        font-size: 12px;
        font-weight: 600;
    }}

    /* Full-width informational strip under overview panels */
    QLabel[role="info_banner"] {{
        background-color: {c['bg_card']};
        color: {c['fg_secondary']};
        border: 1px solid {c['border']};
        border-radius: 6px;
        padding: 10px 14px;
        font-size: 12px;
    }}

    /* === Mode-selection cards (choose_analysis, choose_consensus) ==== */

    _ModeCard {{
        background-color: {c['bg_secondary']};
        border: 2px solid {c['border']};
        border-radius: 10px;
    }}

    _ModeCard[state="hover"], _ModeCard[state="selected"] {{
        border-color: {c['accent_primary']};
    }}

    _ModeCard QLabel {{
        background-color: transparent;
        border: none;
    }}

    /* Circular step-number badge inside _ModeCard */
    _ModeCard QLabel#step_badge {{
        background-color: {c['accent_primary']};
        color: {c['fg_inverse']};
        border-radius: 11px;
        border: none;
        font-weight: bold;
    }}

    /* Workflow "steps" section header sub-title (letter-spaced) */
    _ModeCard QLabel#mode_card_steps_title {{
        color: {c['fg_secondary']};
        letter-spacing: 1.5px;
        background-color: transparent;
        border: none;
    }}

    /* === Tutorial overlay card ====================================== */

    QFrame#tutorial_card {{
        background-color: {c['bg_primary']};
        border: 1px solid {c['border']};
        border-left: 4px solid {c['accent_primary']};
        border-radius: 6px;
    }}

    QFrame#tutorial_card QLabel#tutorial_title {{
        color: {c['fg_primary']};
        font-size: 16px;
        font-weight: bold;
        background-color: transparent;
        border: none;
    }}

    QFrame#tutorial_card QLabel#tutorial_body {{
        color: {c['fg_primary']};
        font-size: 13px;
        background-color: transparent;
        border: none;
    }}

    QFrame#tutorial_card QLabel#tutorial_step_counter {{
        color: {c['fg_secondary']};
        font-size: 11px;
        background-color: transparent;
        border: none;
    }}

    QFrame#tutorial_card QPushButton#tutorial_btn {{
        background-color: transparent;
        color: {c['fg_secondary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        padding: 5px 12px;
        font-size: 12px;
    }}

    QFrame#tutorial_card QPushButton#tutorial_btn_primary {{
        background-color: {c['accent_primary']};
        color: {c['fg_inverse']};
        border: 1px solid {c['accent_primary']};
        border-radius: 4px;
        padding: 5px 16px;
        font-size: 12px;
        font-weight: 600;
    }}

    QFrame#tutorial_card QPushButton#tutorial_btn_auto {{
        background-color: transparent;
        color: {c['accent_primary']};
        border: 1px solid {c['accent_primary']};
        border-radius: 4px;
        padding: 5px 12px;
        font-size: 12px;
        font-weight: 500;
    }}

    /* === Success-bordered "loaded state" panel (setup, de workspace) == */

    QFrame#loaded_frame {{
        background-color: {c['bg_secondary']};
        border: 1px solid {c['success']};
        border-radius: 6px;
    }}

    QFrame#loaded_frame QLabel {{
        background-color: transparent;
        border: none;
    }}

    QLabel#loaded_title {{
        color: {c['success']};
        font-size: 14px;
        font-weight: bold;
    }}

    QLabel#loaded_title[state="error"] {{
        color: {c['error']};
    }}

    /* === Inline image preview (interactive_plots.ImagePreviewWidget) == */

    ImagePreviewWidget {{
        background-color: {c['bg_primary']};
        border: 1px solid {c['border']};
    }}

    /* === Dialog / simulation-benchmark progress stage labels ========= */

    QLabel[role="stage_idle"] {{ color: {c['fg_tertiary']}; }}
    QLabel[role="stage_active"] {{ color: {c['fg_primary']}; }}
    QLabel[role="stage_done_icon"] {{
        color: {c['success']};
        font-weight: bold;
    }}

    /* === Bottom status strip (download_tab tab_status) ============== */

    QLabel#tab_status_strip {{
        padding: 4px 10px;
        border-top: 1px solid {c['border']};
    }}

    /* === Content scroll area (borderless, uses primary bg) ========== */

    QScrollArea[role="content_scroll"] {{
        border: none;
        background-color: {c['bg_primary']};
    }}

    /* === Data-table styling used by result-area tables ============== */

    QTableWidget[role="result_table"] {{
        background-color: {c['bg_primary']};
        gridline-color: {c['border']};
        color: {c['fg_primary']};
        font-size: 12px;
    }}

    QTableWidget[role="result_table"] QHeaderView::section {{
        background-color: {c['bg_secondary']};
        color: {c['fg_primary']};
        border: 1px solid {c['border']};
        padding: 4px;
        font-size: 11px;
        font-weight: bold;
    }}

    /* === Emphasised inline headlines (bold + padding) =============== */

    QLabel[role="headline"] {{
        font-weight: bold;
        padding: 4px;
    }}

    /* === Narrow thumb scrollbars (filter_tab) ======================== */

    QScrollArea[role="thin_scroll"] QScrollBar:vertical {{
        width: 6px;
    }}

    QScrollArea[role="thin_scroll"] QScrollBar::handle:vertical {{
        border-radius: 3px;
    }}

    /* === Large centered italic empty-placeholder ===================== */

    QLabel[role="empty_placeholder"] {{
        color: {c['fg_secondary']};
        font-size: 12px;
        font-style: italic;
        padding: 40px;
        background-color: transparent;
    }}

    /* === Page title ================================================== */

    QLabel[role="page_title"] {{
        font-size: 21px;
        font-weight: bold;
    }}

    /* === Sub-section title (13px bold, used inside cards) ============ */

    QLabel[role="subsection_title"] {{
        font-size: 13px;
        font-weight: bold;
        background-color: transparent;
        border: none;
    }}

    /* === Validation headline (padded rounded bg_secondary banner) ==== */

    QLabel[role="validation_headline"] {{
        padding: 12px;
        border-radius: 6px;
        background-color: {c['bg_secondary']};
    }}

    /* === Success banner (green bold with padding) ==================== */

    QLabel[role="success_banner"] {{
        color: {c['success']};
        font-size: 12px;
        font-weight: bold;
        padding: 6px;
        background-color: transparent;
    }}

    /* === Bold italic text (hint_label in characterise dialog) ======== */

    QLabel[role="italic_hint"] {{
        font-style: italic;
        background-color: transparent;
    }}

    /* === Borderless list widget (dataset picker) ===================== */

    QListWidget[role="borderless"] {{
        border: none;
    }}

    /* === Bold-only title (no colour change) ========================== */

    QLabel[role="bold_label"] {{
        font-weight: bold;
    }}

    /* === Cluster-tab per-colour legend overlay ======================= */

    QFrame#legend_overlay {{
        background-color: {c['bg_secondary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
    }}

    QFrame#legend_overlay QLabel[role="legend_text"] {{
        font-size: 10px;
        background-color: transparent;
        border: none;
    }}

    /* === Cluster-tab compact toggle buttons (color/label/key) ======== */

    QPushButton[role="toggle"] {{
        background-color: {c['bg_tertiary']};
        color: {c['fg_secondary']};
        border: none;
        border-radius: 4px;
        padding: 4px;
        font-size: 11px;
    }}

    QPushButton[role="toggle"]:hover {{
        background-color: {c['bg_secondary']};
    }}

    QPushButton[role="toggle"][state="active"] {{
        background-color: {c['accent_primary']};
        color: {c['fg_inverse']};
        font-weight: bold;
    }}

    QPushButton[role="toggle"][state="warning"] {{
        background-color: {c['warning']};
        color: {c['fg_inverse']};
        font-weight: bold;
    }}

    /* === Transient flash ring on edited spinboxes ==================== */

    QDoubleSpinBox[flash="true"], QSpinBox[flash="true"] {{
        border: 2px solid {c['success']};
    }}

    /* === Cluster-tab embedding-plot container ======================== */

    QWidget[role="embedding_plot"] {{
        background-color: {c['bg_secondary']};
        border: 1px solid {c['border']};
    }}

    QWidget[role="embedding_plot"][hover="true"] {{
        border: 2px solid {c['accent_primary']};
    }}

    /* === Volcano overlay filter button (gene_ma_page) ============== */

    QPushButton#volcano_filter_btn {{
        background-color: rgba(40, 40, 40, 180);
        color: #ddd;
        border: 1px solid #666;
        border-radius: 3px;
        padding: 3px 8px;
        font-size: 12px;
    }}

    QPushButton#volcano_filter_btn:hover {{
        background-color: rgba(60, 60, 60, 200);
    }}

    /* === Compact padded status label (gene_ma_page panel/gwas) ===== */

    QLabel[role="padded_status"] {{
        padding: 4px;
        background-color: transparent;
    }}

    /* === Venn bar plot top spacer ============================== */

    QWidget#venn_top_spacer {{
        background-color: {c['bg_primary']};
    }}

    /* === Tree widget with no outer border (meta workspace) =========== */

    QTreeWidget[role="borderless"] {{
        border: none;
        background-color: transparent;
    }}

    /* === Cluster-tab resolution-sweep thumbnail (default + selected) == */

    _SweepThumbnail {{
        border: 1px solid {c['border']};
        background-color: {c['bg_secondary']};
    }}

    _SweepThumbnail[selected="true"] {{
        border: 2px solid {c['accent_primary']};
    }}

    /* === Compact bordered button with small padding (reload / toggles) */

    QPushButton[role="compact_secondary"] {{
        border: 1px solid {c['border']};
        border-radius: 4px;
        padding: 1px 6px;
        font-size: 11px;
    }}

    QPushButton[role="compact_secondary"]:hover {{
        border-color: {c['accent_primary']};
    }}

    /* Toggleable compact button (qc_tab metric selector) */
    QPushButton[role="metric_toggle"] {{
        border: 1px solid {c['border']};
        border-radius: 3px;
        padding: 1px 6px;
        font-size: 11px;
        background-color: {c['bg_secondary']};
        color: {c['fg_secondary']};
    }}

    QPushButton[role="metric_toggle"]:checked {{
        background-color: {c['accent_primary']};
        color: {c['fg_inverse']};
        border-color: {c['accent_primary']};
    }}

    /* Filter chips: checkable pills that sit above a results table and
       carry the count of rows in each category. Larger than
       'metric_toggle', which is a compact 11px inline switch. */
    QPushButton[role="filter_chip"] {{
        border: 1px solid {c['border']};
        border-radius: 13px;
        padding: 5px 14px;
        font-size: 12px;
        font-weight: 600;
        background-color: {c['bg_tertiary']};
        color: {c['fg_secondary']};
    }}

    QPushButton[role="filter_chip"]:hover {{
        border-color: {c['accent_primary']};
        color: {c['fg_primary']};
    }}

    QPushButton[role="filter_chip"]:checked {{
        background-color: {c['accent_primary']};
        border-color: {c['accent_primary']};
        color: {c['fg_inverse']};
    }}

    QPushButton[role="filter_chip"]:disabled {{
        color: {c['fg_tertiary']};
        border-color: {c['border']};
    }}

    /* === 16px page-level titles ====================================== */

    QLabel[role="page_header"] {{
        font-size: 16px;
        font-weight: bold;
    }}

    /* === Selection header (bold 13px, used in filter_tab) ============ */

    QLabel[role="selection_header"] {{
        font-size: 13px;
        font-weight: bold;
    }}
    """)


def apply_theme(app: QApplication, mode: ThemeMode) -> None:
    """Apply palette + stylesheet to a QApplication."""
    global _current_mode
    _current_mode = mode
    app.setPalette(get_palette())
    app.setStyleSheet(get_stylesheet())


# No-scroll widget variants

class NoScrollComboBox(QComboBox):
    """ComboBox that ignores scroll wheel events, and does not widen its
    parent.

    Qt sizes a combo to its widest entry by default. These sit in the
    fixed-width (320px) sidebars, so one long option -- "Median absolute
    deviation (MAD)" -- pushes the widget past the panel edge. Asking for a
    minimum content length instead lets it take the width it is given and
    clip, rather than dictating the layout.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        # 6 characters, not the widest entry: these live in 320px
        # sidebars beside a label, so anything larger pushes the row
        # past the panel edge. Expanding lets it use whatever width
        # the row actually has.
        self.setMinimumContentsLength(6)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Fixed)

    def wheelEvent(self, event):
        event.ignore()


class NoScrollSpinBox(QSpinBox):
    """SpinBox that ignores scroll wheel events to prevent accidental changes."""
    def wheelEvent(self, event):
        event.ignore()


class NoScrollDoubleSpinBox(QDoubleSpinBox):
    """DoubleSpinBox that ignores scroll wheel events to prevent accidental changes."""
    def wheelEvent(self, event):
        event.ignore()
