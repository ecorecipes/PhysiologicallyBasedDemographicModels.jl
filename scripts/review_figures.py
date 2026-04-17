#!/usr/bin/env python3
"""
Self-consistency figure review via a local VLM (Granite 3.2 Vision).

For each rendered vignette `*.md` in the vignettes tree, locate inline
figure references of the form `![caption](path/to/fig.png)`, extract the
caption plus surrounding narrative, and ask a vision-language model
whether the figure is consistent with its caption and with the text
immediately before and after it.

This is an **advisory** tool, not a test. It writes a Markdown report
to `vignettes/_reports/figure-review.md`. Per the design doc
(`docs/design/vlm_figure_review.md`), we do not attempt pixel-level or
reference-image comparison here because reference paper figures are
under publisher copyright.

Usage (from inside the PBDM package):
    python3 scripts/review_figures.py
    python3 scripts/review_figures.py --model granite3.2-vision
    python3 scripts/review_figures.py --vignette 16_aedes_albopictus
    python3 scripts/review_figures.py --limit 3
"""

from __future__ import annotations

import argparse
import base64
import json
import re
import sys
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path


OLLAMA_URL = "http://localhost:11434/api/generate"

FIGURE_RE = re.compile(r"!\[(?P<caption>[^\]]*)\]\((?P<path>[^)]+\.(?:png|jpg|jpeg|svg))\)")
IMG_HTML_RE = re.compile(
    r'<img\s+[^>]*src="(?P<path>[^"]+\.(?:png|jpg|jpeg|svg))"[^>]*'
    r'(?:alt="(?P<alt>[^"]*)")?[^>]*/?>',
    re.IGNORECASE,
)
FIGCAPTION_RE = re.compile(
    r"<figcaption[^>]*>(?P<cap>.*?)</figcaption>", re.IGNORECASE | re.S
)

PROMPT_TEMPLATE = """You are reviewing a figure from a scientific vignette for self-consistency.

Caption:
{caption}

Narrative before the figure:
{before}

Narrative after the figure:
{after}

Look at the attached figure and answer in four short lines:

1. DESCRIPTION: one sentence describing what the figure actually shows.
2. CAPTION_MATCH: PASS / FAIL / UNSURE — does the figure match its caption?
3. NARRATIVE_MATCH: PASS / FAIL / UNSURE — is the figure consistent with the surrounding narrative?
4. NOTES: one short sentence flagging any obvious problem (wrong variable, empty plot, swapped axes, degenerate output). Write NONE if nothing looks wrong.

Be concise. Ignore styling differences (fonts, colors). Do not invent details that are not visible in the figure."""


@dataclass
class FigureRef:
    vignette: str
    md_path: Path
    image_path: Path
    caption: str
    before: str
    after: str


def strip_md(text: str) -> str:
    """Flatten markdown to plaintext-ish narrative for the prompt."""
    # Drop fenced code blocks.
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    # Drop HTML comments.
    text = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    # Drop image/link wrappers, keep link text.
    text = re.sub(r"!\[[^\]]*\]\([^)]+\)", "", text)
    text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text)
    # Collapse whitespace.
    text = re.sub(r"\s+", " ", text).strip()
    return text


def truncate_words(text: str, n: int) -> str:
    words = text.split()
    if len(words) <= n:
        return text
    return " ".join(words[-n:]) if n > 0 else ""


def truncate_words_head(text: str, n: int) -> str:
    words = text.split()
    if len(words) <= n:
        return text
    return " ".join(words[:n])


def extract_figures(md_path: Path, vignette: str) -> list[FigureRef]:
    content = md_path.read_text(encoding="utf-8", errors="replace")
    figures: list[FigureRef] = []

    matches: list[tuple[int, int, str, str]] = []
    for m in FIGURE_RE.finditer(content):
        cap = (m.group("caption") or "").strip() or "(no caption)"
        matches.append((m.start(), m.end(), m.group("path"), cap))
    for m in IMG_HTML_RE.finditer(content):
        alt = (m.group("alt") or "").strip()
        # Look for a <figcaption> immediately after the img (same div wrapper).
        window = content[m.end(): m.end() + 600]
        fc = FIGCAPTION_RE.search(window)
        cap = strip_md(fc.group("cap")) if fc else alt
        cap = cap or "(no caption)"
        matches.append((m.start(), m.end(), m.group("path"), cap))

    matches.sort()
    for start, end, rel, caption in matches:
        img = (md_path.parent / rel.strip()).resolve()
        if not img.exists():
            continue
        before = strip_md(content[:start])
        after = strip_md(content[end:])
        figures.append(
            FigureRef(
                vignette=vignette,
                md_path=md_path,
                image_path=img,
                caption=caption,
                before=truncate_words(before, 150),
                after=truncate_words_head(after, 80),
            )
        )
    return figures


def call_vlm(model: str, prompt: str, image: Path, timeout: int = 300) -> str:
    img_b64 = base64.b64encode(image.read_bytes()).decode("ascii")
    body = json.dumps(
        {
            "model": model,
            "prompt": prompt,
            "images": [img_b64],
            "stream": False,
            "options": {"temperature": 0.0},
        }
    ).encode("utf-8")
    req = urllib.request.Request(
        OLLAMA_URL, data=body, headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
    return payload.get("response", "").strip()


def collect(vignettes_dir: Path, only: str | None) -> list[FigureRef]:
    figures: list[FigureRef] = []
    for d in sorted(vignettes_dir.iterdir()):
        if not d.is_dir() or d.name.startswith((".", "_")):
            continue
        if only and d.name != only:
            continue
        md = d / f"{d.name}.md"
        if not md.exists():
            continue
        figures.extend(extract_figures(md, d.name))
    return figures


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", default="granite3.2-vision")
    ap.add_argument("--vignette", default=None,
                    help="Restrict to a single vignette directory name")
    ap.add_argument("--limit", type=int, default=0,
                    help="Limit number of figures reviewed (0 = all)")
    ap.add_argument("--report",
                    default="vignettes/_reports/figure-review.md")
    ap.add_argument("--timeout", type=int, default=300)
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent
    vignettes_dir = root / "vignettes"
    report_path = root / args.report
    report_path.parent.mkdir(parents=True, exist_ok=True)

    figures = collect(vignettes_dir, args.vignette)
    if args.limit > 0:
        figures = figures[: args.limit]

    print(f"[review] model={args.model} figures={len(figures)}")

    rows: list[str] = []
    rows.append("# Vignette figure self-consistency review\n")
    rows.append(f"- Model: `{args.model}`\n")
    rows.append(f"- Figures reviewed: {len(figures)}\n")
    rows.append("- Scope: self-consistency only (caption + surrounding "
                "narrative). No reference-paper comparison.\n")

    detail_rows: list[str] = ["\n## Results\n"]

    cap_counts = {"PASS": 0, "FAIL": 0, "UNSURE": 0, "?": 0}
    nar_counts = {"PASS": 0, "FAIL": 0, "UNSURE": 0, "?": 0}
    flagged_fails: list[str] = []
    flagged_notes: list[tuple[str, str]] = []

    t0 = time.time()
    for i, fig in enumerate(figures, 1):
        prompt = PROMPT_TEMPLATE.format(
            caption=fig.caption,
            before=fig.before or "(no narrative)",
            after=fig.after or "(no narrative)",
        )
        print(f"[review] {i}/{len(figures)} {fig.vignette} "
              f"{fig.image_path.name} ...", flush=True)
        try:
            answer = call_vlm(args.model, prompt, fig.image_path, args.timeout)
        except Exception as e:
            answer = f"ERROR: {e}"
        try:
            rel_img = fig.image_path.relative_to(vignettes_dir)
        except ValueError:
            rel_img = fig.image_path
        header = f"{fig.vignette} — `{rel_img}`"
        detail_rows.append(f"\n### {header}\n")
        detail_rows.append(f"**Caption:** {fig.caption}\n\n")
        detail_rows.append("```\n")
        detail_rows.append(answer.strip() + "\n")
        detail_rows.append("```\n")

        cm = re.search(r"CAPTION_MATCH:\s*(PASS|FAIL|UNSURE)", answer, re.I)
        nm = re.search(r"NARRATIVE_MATCH:\s*(PASS|FAIL|UNSURE)", answer, re.I)
        cap_counts[cm.group(1).upper() if cm else "?"] += 1
        nar_counts[nm.group(1).upper() if nm else "?"] += 1
        if cm and cm.group(1).upper() == "FAIL":
            flagged_fails.append(f"{header} (caption)")
        if nm and nm.group(1).upper() == "FAIL":
            flagged_fails.append(f"{header} (narrative)")
        note_match = re.search(r"NOTES:\s*(.+?)(?:\n|$)", answer)
        if note_match:
            note = note_match.group(1).strip().rstrip("`").strip()
            if note and not note.upper().startswith("NONE") \
                    and note.upper() not in ("N/A", "-"):
                flagged_notes.append((header, note))

    rows.append("\n## Summary\n\n")
    rows.append(f"- Total figures reviewed: **{sum(cap_counts.values())}**\n")
    rows.append(
        f"- Caption match: PASS={cap_counts['PASS']} "
        f"FAIL={cap_counts['FAIL']} UNSURE={cap_counts['UNSURE']}\n"
    )
    rows.append(
        f"- Narrative match: PASS={nar_counts['PASS']} "
        f"FAIL={nar_counts['FAIL']} UNSURE={nar_counts['UNSURE']}\n"
    )
    rows.append(f"- Figures with non-NONE notes: **{len(flagged_notes)}**\n")
    rows.append(f"- Figures with any FAIL: **{len(flagged_fails)}**\n")
    if flagged_fails:
        rows.append("\n### Flagged FAILs\n\n")
        for line in flagged_fails:
            rows.append(f"- {line}\n")
    if flagged_notes:
        rows.append("\n### Flagged notes\n\n")
        for header, note in flagged_notes:
            rows.append(f"- {header}: {note}\n")
    rows.extend(detail_rows)

    elapsed = time.time() - t0
    rows.append(f"\n---\n_Elapsed: {elapsed:.1f}s_\n")
    report_path.write_text("".join(rows), encoding="utf-8")
    print(f"[review] wrote {report_path} ({elapsed:.1f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
