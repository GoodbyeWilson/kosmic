"""The app's type ladder.

The stylesheet had grown ten distinct font sizes across two units doing
overlapping jobs -- 11px and 12px both "small", 12px and 13px both "body",
9pt and 8pt duplicating them again. These pin the collapsed ladder so a new
rule cannot quietly reintroduce a size nobody chose.
"""
import re

from kosmic.gui.shared.theme import (
    _FONT_LADDER_PX,
    get_stylesheet,
    normalise_font_sizes,
)

#: Every size the shipped sheet is allowed to use.
ALLOWED_PX = {12, 14, 16, 18, 22, 24}


def _sizes(qss):
    return re.findall(r"font-size:\s*([0-9.]+)(px|pt)", qss)


def test_the_shipped_sheet_uses_only_ladder_sizes():
    sizes = _sizes(get_stylesheet())
    assert sizes, "stylesheet declares no font sizes at all"
    off_ladder = {f"{v}{u}" for v, u in sizes
                  if u != "px" or int(float(v)) not in ALLOWED_PX}
    assert not off_ladder, f"sizes outside the ladder: {sorted(off_ladder)}"


def test_the_whole_sheet_is_one_unit():
    """Mixing px and pt was half the inconsistency: 9pt and 11px are the
    same size but scale differently."""
    assert {u for _v, u in _sizes(get_stylesheet())} == {"px"}


def test_neighbouring_source_sizes_collapse_together():
    """11px and 12px were both 'small'; they must not stay one apart."""
    assert _FONT_LADDER_PX[11] == _FONT_LADDER_PX[10]
    assert _FONT_LADDER_PX[12] == _FONT_LADDER_PX[13]
    assert _FONT_LADDER_PX[14] == _FONT_LADDER_PX[15]


def test_the_ladder_is_monotonic():
    """A caption must never come out larger than a heading."""
    pairs = sorted(_FONT_LADDER_PX.items())
    mapped = [new for _old, new in pairs]
    assert mapped == sorted(mapped)


def test_scale_multiplies_the_ladder_not_the_source():
    out = normalise_font_sizes("a{font-size: 12px;}", 2.0)
    assert "font-size: 28px" in out          # 12 -> ladder 14 -> x2


def test_a_size_not_on_the_ladder_survives():
    """An unmapped size is passed through, never dropped."""
    assert "font-size: 37px" in normalise_font_sizes("a{font-size: 37px;}", 1.0)


def test_sizes_never_fall_below_the_readable_floor():
    assert "font-size: 8px" in normalise_font_sizes("a{font-size: 12px;}", 0.01)


def test_whitespace_variants_are_matched():
    assert "font-size: 14px" in normalise_font_sizes("a{font-size:13px;}")
