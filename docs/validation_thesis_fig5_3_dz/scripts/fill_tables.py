#!/usr/bin/env python3
"""Copy the generated table fragments (scripts/tables/NAME.md) into README.md, stages.md and
appendix.md between <!-- BEGIN table:NAME --> and <!-- END table:NAME --> (blank line kept on
both sides so GitHub renders the table)."""
from __future__ import annotations

import re
import sys

from common import DOC, TABLES

PATTERN = re.compile(r"(<!-- BEGIN table:(\S+) -->\n)(.*?)(<!-- END table:\2 -->)", re.S)


def main():
    missing = []
    for name in ("README.md", "stages.md", "appendix.md"):
        path = DOC / name
        if not path.exists():
            continue
        text = path.read_text()

        def repl(m):
            frag = TABLES / f"{m.group(2)}.md"
            if not frag.exists():
                missing.append(f"{name}: {m.group(2)}")
                return m.group(0)
            return m.group(1) + "\n" + frag.read_text().rstrip("\n") + "\n\n" + m.group(4)
        new = PATTERN.sub(repl, text)
        if new != text:
            path.write_text(new)
        print(name, "tables:", len(PATTERN.findall(new)))
    if missing:
        print("missing fragments:", missing)
        sys.exit(1)


if __name__ == "__main__":
    main()
