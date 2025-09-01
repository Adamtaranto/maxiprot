#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
maxiprot CLI entrypoint.

Subcommands
-----------
- maxiprot filter   -> run the alignment filtering/selection tool (maxiprot.filter)
- maxiprot extract  -> run the sequence extraction tool (maxiprot.extract)

All arguments after the subcommand are forwarded unchanged to the corresponding
module's `main(argv)` function.

Examples
--------
    maxiprot filter --help
    maxiprot extract --help

    # Typical usage
    maxiprot filter input.gff3 --out-gff3 best.gff3 --out-tsv best.tsv
    maxiprot extract best.gff3 -g genome.fa --extract protein > best.faa
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional, Sequence

from maxiprot._version import __version__


def _load_filter_main():
    """
    Import and return maxiprot.filter.main.
    """
    try:
        from maxiprot.filter import main as filter_main  # type: ignore
    except Exception as e:  # broad to show helpful message
        raise ModuleNotFoundError(
            "Failed to import 'maxiprot.filter'. Ensure 'src/maxiprot/filter.py' "
            'is on the Python path and exposes a callable main(argv) -> int.'
        ) from e
    if not callable(filter_main):
        raise TypeError('maxiprot.filter.main is not callable')
    return filter_main


def _load_extract_main():
    """
    Import and return maxiprot.extract.main.
    """
    try:
        from maxiprot.extract import main as extract_main  # type: ignore
    except Exception as e:
        raise ModuleNotFoundError(
            "Failed to import 'maxiprot.extract'. Ensure 'src/maxiprot/extract.py' "
            'is on the Python path and exposes a callable main(argv) -> int.'
        ) from e
    if not callable(extract_main):
        raise TypeError('maxiprot.extract.main is not callable')
    return extract_main


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level CLI parser that selects a subcommand."""
    parser = argparse.ArgumentParser(
        prog='maxiprot',
        description='Unified CLI for maxiprot tools.',
        epilog="Use 'maxiprot filter --help' or 'maxiprot extract --help' for subcommand options.",
        add_help=True,
    )
    parser.add_argument(
        '--version', action='version', version=f'maxiprot {__version__}'
    )
    subparsers = parser.add_subparsers(
        dest='command',
        metavar='{filter,extract}',
        required=True,
        help='Subcommand to run',
    )

    # Minimal subparsers; actual options belong to the downstream tools.
    subparsers.add_parser(
        'filter',
        help='Select the best miniprot alignment per locus and emit GFF3/TSV.',
        add_help=False,  # let maxiprot.filter handle its own --help
    )
    subparsers.add_parser(
        'extract',
        help='Extract protein/CDS/gene sequences from maxiprot GFF3.',
        add_help=False,  # let maxiprot.extract handle its own --help
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Entrypoint. Parse the subcommand token and forward remaining args.
    """
    if argv is None:
        argv = sys.argv[1:]

    parser = build_parser()
    # Only parse the subcommand; forward the rest (including --help) to the sub-tool.
    args, remainder = parser.parse_known_args(argv)

    if args.command == 'filter':
        try:
            filter_main = _load_filter_main()
        except (ModuleNotFoundError, TypeError) as e:
            print(str(e), file=sys.stderr)
            return 2
        return int(filter_main(remainder))

    if args.command == 'extract':
        try:
            extract_main = _load_extract_main()
        except (ModuleNotFoundError, TypeError) as e:
            print(str(e), file=sys.stderr)
            return 2
        return int(extract_main(remainder))

    # Should not happen (subparsers.required=True), but keep a fallback:
    parser.print_usage(sys.stderr)
    return 2


if __name__ == '__main__':
    raise SystemExit(main())
