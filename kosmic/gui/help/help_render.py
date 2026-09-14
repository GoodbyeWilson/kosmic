# Markdown -> styled HTML rendering for the help system. Admonitions
# use blockquote syntax with a leading '**Note:**' / '**Tip:**' /
# '**Warning:**' / '**Danger:**' -- coloured left border by kind.
from __future__ import annotations

import re
from functools import lru_cache

from markdown_it import MarkdownIt

from kosmic.gui.shared.theme import get_color


_md = MarkdownIt("commonmark", {"breaks": False, "linkify": True}).enable("table")


def _stylesheet() -> str:
    """Return the CSS for help content. Pulls colours from the active theme."""
    fg = get_color("fg_primary")
    fg2 = get_color("fg_secondary")
    accent = get_color("accent_primary")
    code_bg = get_color("bg_secondary")  # distinct from QTextBrowser bg_primary
    border = get_color("border")
    return f"""
    body {{
        font-family: "Segoe UI", "Helvetica Neue", sans-serif;
        font-size: 10pt;
        color: {fg};
        line-height: 1.4;
    }}

    h1 {{
        font-size: 13pt;
        font-weight: 600;
        color: {fg};
        margin: 0 0 6px 0;
    }}
    h2 {{
        font-size: 11.5pt;
        font-weight: 600;
        color: {fg};
        margin: 16px 0 4px 0;
    }}
    h3 {{
        font-size: 10.5pt;
        font-weight: 600;
        color: {fg};
        margin: 12px 0 3px 0;
    }}
    h4 {{
        font-size: 9.5pt;
        font-weight: 600;
        color: {fg2};
        margin: 10px 0 2px 0;
        text-transform: uppercase;
        letter-spacing: 0.6px;
    }}

    p {{
        margin: 0 0 8px 0;
    }}

    a {{
        color: {accent};
        text-decoration: none;
    }}
    a:hover {{
        text-decoration: underline;
    }}

    ul, ol {{
        margin: 0 0 8px 0;
        padding-left: 20px;
    }}
    li {{
        margin: 0 0 2px 0;
    }}

    code {{
        font-family: "Consolas", "Courier New", monospace;
        background-color: {code_bg};
        color: {accent};
        padding: 1px 5px;
        border-radius: 3px;
        font-size: 9.5pt;
    }}

    pre {{
        background-color: {code_bg};
        border: 1px solid {border};
        border-radius: 4px;
        padding: 10px 12px;
        margin: 0 0 14px 0;
        font-family: "Consolas", "Courier New", monospace;
        font-size: 9.5pt;
        color: {fg};
    }}
    pre code {{
        background-color: transparent;
        color: {fg};
        padding: 0;
        border: none;
    }}

    blockquote {{
        margin: 14px 0;
        padding: 10px 14px;
        border-left: 3px solid {accent};
        background-color: {code_bg};
        color: {fg};
    }}
    blockquote p {{
        margin: 0 0 6px 0;
    }}
    blockquote p:last-child {{
        margin-bottom: 0;
    }}

    table {{
        border-collapse: collapse;
        margin: 14px 0;
        font-size: 10pt;
    }}
    th, td {{
        border: 1px solid {border};
        padding: 6px 12px;
        text-align: left;
    }}
    th {{
        background-color: {code_bg};
        color: {fg};
        font-weight: 600;
    }}

    hr {{
        border: none;
        border-top: 1px solid {border};
        margin: 20px 0;
    }}

    strong {{
        font-weight: 600;
        color: {fg};
    }}
    em {{
        font-style: italic;
        color: {fg2};
    }}
    """


_ADMON_KIND_COLOURS = {
    "tip":     "#3aa563",
    "note":    "#3a7bd5",
    "warning": "#d59333",
    "danger":  "#d54545",
}


def _wrap_admonitions_html(html: str) -> str:
    """Tag rendered '<blockquote>'s by their leading keyword.

    Looks at the first '<strong>...:</strong>' inside each '<blockquote>'
    and re-emits it with a 'class="admon-<kind>"' attribute so the
    stylesheet can colour the left border accordingly. Cheap regex pass --
    Qt's stylesheet can't do attribute selectors anyway, so we inline a
    coloured border directly.
    """
    pattern = re.compile(
        r"<blockquote>\s*<p>\s*<strong>(Tip|Note|Warning|Danger):</strong>",
        re.IGNORECASE,
    )

    def repl(m: re.Match) -> str:
        kind = m.group(1).lower()
        colour = _ADMON_KIND_COLOURS.get(kind, "#888")
        return (
            f'<blockquote style="border-left: 3px solid {colour};">'
            f'<p><strong style="color:{colour};">{m.group(1)}:</strong>'
        )

    return pattern.sub(repl, html)


_IMG_RE = re.compile(r"<img src=" + chr(34))  # '<img src="'


def _size_images(html: str, width: int) -> str:
    """Give every image an explicit width.

    Qt's rich text has no max-width, so a 1600px screenshot in a 400px
    side panel would scroll sideways. Width scales height with it."""
    return _IMG_RE.sub(f'<img width="{width}" src=' + chr(34), html)


def render(markdown_text: str, image_width: int | None = None) -> str:
    """Render markdown to a complete HTML document with embedded stylesheet.

    The leading H1 is stripped because the panel header already shows it.
    'image_width' pins embedded images to that many pixels (the viewer's
    usable width); None leaves them at natural size.
    """
    text = markdown_text or ""
    stripped = text.lstrip()
    if stripped.startswith("# "):
        nl = stripped.index("\n")
        text = stripped[nl + 1:].lstrip("\n")
    body = _md.render(text)
    body = _wrap_admonitions_html(body)
    if image_width:
        body = _size_images(body, image_width)
    return (
        "<html><head><style>"
        + _stylesheet()
        + "</style></head><body>"
        + body
        + "</body></html>"
    )


@lru_cache(maxsize=128)
def render_cached(markdown_text: str, image_width: int | None = None) -> str:
    """Cached renderer for content that may be displayed repeatedly."""
    return render(markdown_text, image_width)
