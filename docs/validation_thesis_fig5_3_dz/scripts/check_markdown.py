#!/usr/bin/env python3
"""Markdown and folder checks of the report (exit 1 on any problem):
- links and images resolve, and their targets are not git-ignored;
- every `#item-N` / `#anchor` target exists (explicit <a id> or a GitHub heading slug);
- table rows have the same number of cells as their header (escaped \\| not counted);
- <details>: a blank line after </summary> and before </details>; no <details> inside a list item;
- no `$` math; `runs/...` and `~/...` paths are not links;
- every file in the folder: no trailing whitespace (text files), final newline, <= 5 MB, not
  git-ignored, no __pycache__, binaries only PNG and .gz."""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

from common import DOC, REPO

MD = ["README.md", "stages.md", "appendix.md", "scripts/criteria.md", "scripts/STATUS.md", "scripts/NOTES.md",
      "scripts/review_log.md"]
problems = []


def slug(heading: str) -> str:
    s = heading.strip().lower()
    s = re.sub(r"<[^>]+>", "", s)
    s = re.sub(r"[^\w\- ]", "", s)
    return s.replace(" ", "-")


def anchors(path: Path):
    text = path.read_text()
    a = set(re.findall(r'<a id="([^"]+)"></a>', text))
    for line in text.splitlines():
        m = re.match(r"^#+\s+(.*)$", line)
        if m:
            a.add(slug(m.group(1)))
    return a


def ignored(path: Path) -> bool:
    r = subprocess.run(["git", "check-ignore", "-q", str(path)], cwd=REPO)
    return r.returncode == 0


def check_md(rel: str):
    path = DOC / rel
    if not path.exists():
        return
    text = path.read_text()
    lines = text.splitlines()
    in_code = False
    header_cells = None
    for i, line in enumerate(lines, 1):
        if line.startswith("```"):
            in_code = not in_code
            continue
        if in_code:
            continue
        if "$" in re.sub(r"`[^`]*`", "", line):
            problems.append(f"{rel}:{i}: '$' outside code")
        # table cell counts
        if line.startswith("|"):
            cells = len(re.findall(r"(?<!\\)\|", line)) - 1
            if header_cells is None:
                header_cells = cells
            elif cells != header_cells and not re.match(r"^\|[-|]+\|$", line):
                problems.append(f"{rel}:{i}: {cells} cells, header has {header_cells}")
        else:
            header_cells = None
        if line.strip() == "</summary>" or line.rstrip().endswith("</summary>"):
            if i < len(lines) and lines[i].strip() != "":
                problems.append(f"{rel}:{i}: no blank line after </summary>")
        if line.strip() == "</details>" and i >= 2 and lines[i - 2].strip() != "":
            problems.append(f"{rel}:{i}: no blank line before </details>")
        if re.match(r"^\s*[-*\d.]+\s+.*<details>", line):
            problems.append(f"{rel}:{i}: <details> inside a list")
        # links and images
        for m in re.finditer(r"!?\[[^\]]*\]\(([^)\s]+)\)", re.sub(r"`[^`]*`", "", line)):
            target = m.group(1)
            if target.startswith("http"):
                continue
            if target.startswith(("runs/", "~/", "/")):
                problems.append(f"{rel}:{i}: link to a local path {target}")
                continue
            file_part, _, frag = target.partition("#")
            dest = (path.parent / file_part).resolve() if file_part else path
            if not dest.exists():
                problems.append(f"{rel}:{i}: missing target {target}")
                continue
            if ignored(dest):
                problems.append(f"{rel}:{i}: target is git-ignored {target}")
            if frag and dest.suffix == ".md" and frag not in anchors(dest):
                problems.append(f"{rel}:{i}: missing anchor #{frag} in {dest.name}")


def check_files():
    for p in sorted(DOC.rglob("*")):
        if p.is_dir():
            if p.name == "__pycache__":
                problems.append(f"{p.relative_to(DOC)}: __pycache__")
            continue
        rel = p.relative_to(DOC)
        if p.stat().st_size > 5 * 1024 * 1024:
            problems.append(f"{rel}: larger than 5 MB")
        if ignored(p):
            problems.append(f"{rel}: git-ignored")
        if p.suffix in (".png", ".gz"):
            continue
        try:
            text = p.read_text()
        except UnicodeDecodeError:
            problems.append(f"{rel}: binary file")
            continue
        if text and not text.endswith("\n"):
            problems.append(f"{rel}: no final newline")
        if text.endswith("\n\n"):
            problems.append(f"{rel}: blank line at end of file")
        for i, line in enumerate(text.splitlines(), 1):
            if line != line.rstrip():
                problems.append(f"{rel}:{i}: trailing whitespace")
                break


def main():
    for rel in MD:
        check_md(rel)
    check_files()
    for p in problems[:40]:
        print(p)
    print(f"check_markdown: {len(problems)} problem(s)")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
