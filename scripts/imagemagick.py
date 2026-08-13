"""Portable ImageMagick command selection for render validators."""

from __future__ import annotations

import shutil


def convert_command(*arguments: str) -> list[str]:
    """Return an ImageMagick conversion command for versions 6 and 7."""
    if shutil.which("magick") is not None:
        return ["magick", *arguments]
    if shutil.which("convert") is not None:
        return ["convert", *arguments]
    raise FileNotFoundError(
        "ImageMagick is required: install either the 'magick' (v7) or "
        "'convert' (v6) command"
    )
