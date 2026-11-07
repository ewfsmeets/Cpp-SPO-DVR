#!/usr/bin/env python3
"""Item 17, step 1 (needs pymupdf; run as `uv run --no-project --with pymupdf python -I thesis_extract.py PDF`):
find Figure 5.3 of M. F. Somers' PhD thesis by its caption, dump the page's vector drawings to
JSON and render the figure region (for the check image). Output in
runs/thesis_fig5_3_dz_validation/thesis/ (git-ignored)."""
import json
import sys
from pathlib import Path

import pymupdf

OUT = Path(__file__).resolve().parents[3] / "runs" / "thesis_fig5_3_dz_validation" / "thesis"
CAPTION = "Figure 5.3: Reaction probabilities computed for the DZ PES"


def main():
    doc = pymupdf.open(sys.argv[1])
    pages = [i for i, p in enumerate(doc) if CAPTION in p.get_text().replace("\n", " ")]
    assert len(pages) == 1, pages
    page = doc[pages[0]]
    drawings = []
    for x in page.get_drawings():
        items = []
        for it in x["items"]:
            if it[0] == "re":
                r = it[1]
                pts = [[r.x0, r.y0], [r.x1, r.y1]]
            else:
                pts = [[p.x, p.y] for p in it[1:] if hasattr(p, "x")]
            items.append([it[0], pts])
        r = x["rect"]
        drawings.append({"type": x["type"], "rect": [r.x0, r.y0, r.x1, r.y1], "fill": x.get("fill"),
                         "color": x.get("color"), "width": x.get("width"), "items": items})
    OUT.mkdir(parents=True, exist_ok=True)
    clip = pymupdf.Rect(140, 285, 500, 670)
    dpi = 200
    page.get_pixmap(dpi=dpi, clip=clip).save(OUT / "fig5_3_crop.png")
    meta = {"pdf_page_index": pages[0], "printed_page": pages[0] + 1, "clip_pt": [clip.x0, clip.y0, clip.x1, clip.y1],
            "dpi": dpi, "pymupdf": pymupdf.__version__, "n_drawings": len(drawings)}
    json.dump({"meta": meta, "drawings": drawings}, open(OUT / "fig5_3_drawings.json", "w"))
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
