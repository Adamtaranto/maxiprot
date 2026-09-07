"""Shared I/O helpers for maxiprot subcommands."""

from __future__ import annotations

from collections.abc import Callable, Iterator
import sys
from typing import TextIO


def open_output(path: str | None) -> tuple[TextIO, bool]:
    """
    Open an output handle for text writing.

    Parameters
    ----------
    path : str or None
        Destination path. ``None``, ``''`` and ``'-'`` mean stdout.

    Returns
    -------
    tuple of (TextIO, bool)
        The handle and whether the caller must close it.
    """
    if path in (None, '-', ''):
        return sys.stdout, False
    return open(path, 'w', encoding='utf-8'), True


def read_lines(path: str | None) -> Iterator[str]:
    """
    Yield lines from a text file or stdin.

    Parameters
    ----------
    path : str or None
        Input path. ``None``, ``''`` and ``'-'`` mean stdin.

    Yields
    ------
    str
        One line at a time (newline preserved).
    """
    if path in (None, '-', ''):
        yield from sys.stdin
    else:
        with open(path, 'r', encoding='utf-8') as fh:
            yield from fh


def guard_broken_pipe(fn: Callable[[], int]) -> int:
    """
    Run ``fn`` and exit quietly on ``BrokenPipeError``.

    Allows piping output into ``head``/``tail`` without a traceback.

    Parameters
    ----------
    fn : callable
        Zero-argument callable returning an exit code.

    Returns
    -------
    int
        The exit code from ``fn``, or 1 if the pipe was broken.
    """
    try:
        return fn()
    except BrokenPipeError:
        for stream in (sys.stdout, sys.stderr):
            try:
                stream.close()
            except Exception:
                pass
        return 1
