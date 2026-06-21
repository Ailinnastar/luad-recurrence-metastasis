#!/usr/bin/env python3
"""Build public GitHub Pages site from the final report (no student contributions)."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT_SRC = ROOT / "THE FINAL REPORT" / "report.html"
DOCS = ROOT / "docs"


def strip_student_contributions(html: str) -> str:
    html = re.sub(
        r'\s*<li><a href="#student-contributions"[^>]*>.*?</a></li>\s*',
        "",
        html,
        flags=re.DOTALL,
    )
    html = re.sub(
        r'<section id="student-contributions" class="level1"[^>]*>.*?</section>\s*',
        "",
        html,
        flags=re.DOTALL,
    )
    html = re.sub(
        r'<meta name="author" content="[^"]*">\s*',
        "",
        html,
    )
    html = re.sub(
        r'\s*<div>\s*<div class="quarto-title-meta-heading">Author</div>\s*'
        r'<div class="quarto-title-meta-contents">\s*<p>[^<]*</p>\s*</div>\s*</div>\s*',
        "",
        html,
        flags=re.DOTALL,
    )
    html = re.sub(
        r'<span class="co"> &quot;550167750, 530616715, 540094367, 530598798, 540673139&quot;</span></span>\s*',
        "",
        html,
    )
    html = re.sub(
        r'<span id="cb12-5">.*?</span>\s*',
        "",
        html,
        count=1,
        flags=re.DOTALL,
    )
    html = re.sub(
        r'<span id="cb12-647">.*?<span id="cb12-682">',
        r'<span id="cb12-682">',
        html,
        count=1,
        flags=re.DOTALL,
    )
    # Renumber top-level sections after removal
    renumber = [
        ("references", "7", "6"),
        ("ai-declaration", "8", "7"),
        ("appendix", "9", "8"),
    ]
    for section_id, old_num, new_num in renumber:
        html = html.replace(
            f'id="toc-{section_id}"><span class="header-section-number">{old_num}</span>',
            f'id="toc-{section_id}"><span class="header-section-number">{new_num}</span>',
        )
        html = html.replace(
            f'id="{section_id}" class="level1" data-number="{old_num}"',
            f'id="{section_id}" class="level1" data-number="{new_num}"',
        )
        html = html.replace(
            f'data-number="{old_num}"><span class="header-section-number">{old_num}</span>',
            f'data-number="{new_num}"><span class="header-section-number">{new_num}</span>',
            1,
        )
    # Appendix subsection numbers: 9.x -> 8.x
    html = re.sub(
        r'(<section id="appendix"[\s\S]*?</section>\s*</main>)',
        lambda m: m.group(1).replace("header-section-number\">9.", "header-section-number\">8."),
        html,
        count=1,
    )
    return html


def main() -> None:
    DOCS.mkdir(exist_ok=True)
    html = REPORT_SRC.read_text(encoding="utf-8")
    (DOCS / "report.html").write_text(strip_student_contributions(html), encoding="utf-8")
    print(f"Wrote {DOCS / 'report.html'}")


if __name__ == "__main__":
    main()
