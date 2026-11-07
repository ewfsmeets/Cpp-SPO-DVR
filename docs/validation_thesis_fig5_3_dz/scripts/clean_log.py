#!/usr/bin/env python3
"""Clean a `script` log: CRLF and \\r progress redraws split into lines, ANSI escapes removed,
progress-bar fragments dropped, trailing whitespace stripped.

Usage: clean_log.py RAW_LOG CLEAN_OUT
"""
import re
import sys

ANSI = re.compile(r"\x1b\[[0-9;?]*[A-Za-z]")


def clean_lines(raw: str):
    for chunk in raw.split("\n"):
        for piece in chunk.split("\r"):
            line = ANSI.sub("", piece).rstrip()
            if "▕" in line:  # progress bar fragment (▕)
                continue
            yield line


def main():
    raw = open(sys.argv[1], encoding="utf-8", errors="replace").read()
    out = []
    previous_empty = False
    for line in clean_lines(raw):
        # \r splits create runs of empty lines: keep at most one
        if line == "":
            if previous_empty:
                continue
            previous_empty = True
        else:
            previous_empty = False
        out.append(line)
    with open(sys.argv[2], "w", encoding="utf-8") as f:
        f.write("\n".join(out).rstrip("\n") + "\n")


if __name__ == "__main__":
    main()
