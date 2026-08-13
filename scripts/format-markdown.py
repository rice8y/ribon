#!/usr/bin/env python3
"""Keep each Markdown prose paragraph on one physical source line."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


FENCE_RE = re.compile(r"^ {0,3}(`{3,}|~{3,})")
LIST_RE = re.compile(r"^(\s*)([-+*]|\d+[.)])(\s+)(.*)$")
QUOTE_RE = re.compile(r"^( {0,3}> ?)(.*)$")
EMPTY_IMAGE_ALT_RE = re.compile(r"!\[\s*\]\(")
HTML_IMAGE_WITHOUT_ALT_RE = re.compile(r"<img\b(?![^>]*\balt\s*=)[^>]*>", re.IGNORECASE)


def is_standalone(line: str) -> bool:
    stripped = line.strip()
    return bool(
        re.match(r"^#{1,6}(?:\s|$)", stripped)
        or re.match(r"^(?:={3,}|-{3,})$", stripped)
        or re.match(r"^(?:\*\s*){3,}$", stripped)
        or re.match(r"^\[[^]]+\]:\s*", stripped)
        or stripped.startswith("|")
        or stripped.endswith("|")
        or re.match(r"^:?-{3,}:?(?:\s*\|\s*:?-{3,}:?)+$", stripped)
        or re.match(r"^</?[A-Za-z][^>]*>", stripped)
        or stripped.startswith("<!--")
        or stripped.startswith("-->")
        or stripped.startswith("![")
        or line.startswith("    ")
        or line.startswith("\t")
    )


def format_markdown(source: str) -> str:
    lines = source.splitlines()
    output: list[str] = []
    paragraph: list[str] = []
    paragraph_prefix = ""
    paragraph_kind = "plain"
    fence_char: str | None = None
    fence_length = 0

    def flush() -> None:
        nonlocal paragraph, paragraph_prefix, paragraph_kind
        if paragraph:
            output.append(paragraph_prefix + " ".join(part.strip() for part in paragraph))
        paragraph = []
        paragraph_prefix = ""
        paragraph_kind = "plain"

    for line in lines:
        fence = FENCE_RE.match(line)
        if fence_char is not None:
            output.append(line)
            marker = fence.group(1) if fence else ""
            if marker.startswith(fence_char) and len(marker) >= fence_length:
                fence_char = None
                fence_length = 0
            continue

        if fence:
            flush()
            marker = fence.group(1)
            fence_char = marker[0]
            fence_length = len(marker)
            output.append(line)
            continue

        if not line.strip():
            flush()
            output.append("")
            continue

        quote = QUOTE_RE.match(line)
        if quote:
            prefix, body = quote.groups()
            if paragraph and (paragraph_kind != "quote" or paragraph_prefix != prefix):
                flush()
            paragraph_kind = "quote"
            paragraph_prefix = prefix
            paragraph.append(body)
            continue

        item = LIST_RE.match(line)
        if item:
            flush()
            indent, marker, spacing, body = item.groups()
            paragraph_kind = "list"
            paragraph_prefix = indent + marker + spacing
            paragraph.append(body)
            continue

        if is_standalone(line):
            flush()
            output.append(line)
            continue

        if paragraph_kind in ("plain", "list"):
            paragraph.append(line)
        else:
            flush()
            paragraph.append(line)

    flush()
    suffix = "\n" if source.endswith("\n") else ""
    return "\n".join(output) + suffix


def repository_markdown_files() -> list[Path]:
    result = subprocess.run(
        ["rg", "--files", "-g", "*.md"],
        check=True,
        capture_output=True,
        text=True,
    )
    return [Path(line) for line in result.stdout.splitlines() if line]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("paths", nargs="*", type=Path)
    args = parser.parse_args()
    paths = args.paths or repository_markdown_files()
    changed: list[Path] = []
    invalid_images: list[Path] = []

    for path in paths:
        source = path.read_text(encoding="utf-8")
        if path.name.lower() == "readme.md" and (
            EMPTY_IMAGE_ALT_RE.search(source) or HTML_IMAGE_WITHOUT_ALT_RE.search(source)
        ):
            invalid_images.append(path)
        formatted = format_markdown(source)
        if formatted == source:
            continue
        changed.append(path)
        if not args.check:
            path.write_text(formatted, encoding="utf-8")

    if args.check and (changed or invalid_images):
        for path in changed:
            print(f"paragraph wraps remain: {path}", file=sys.stderr)
        for path in invalid_images:
            print(f"README image without alt text: {path}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
