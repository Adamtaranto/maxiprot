#!/usr/bin/env python3
"""
maxiprot CLI entrypoint.

Subcommands
-----------
- maxiprot filter   -> select the best miniprot alignment per locus (maxiprot.filter)
- maxiprot extract  -> extract protein/CDS/gene sequences (maxiprot.extract)

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

from maxiprot import __version__, extract, filter
from maxiprot.ioutils import guard_broken_pipe


def build_parser() -> argparse.ArgumentParser:
    """
    Build the top-level CLI parser with real subcommand parsers.

    Returns
    -------
    argparse.ArgumentParser
        Parser whose subcommands dispatch via ``args.func(args)``.
    """
    parser = argparse.ArgumentParser(
        prog='maxiprot',
        description='Unified CLI for maxiprot tools.',
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

    sp_filter = subparsers.add_parser(
        'filter',
        help='Select the best miniprot alignment per locus and emit GFF3/TSV.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    filter.add_arguments(sp_filter)
    sp_filter.add_argument(
        '--version', action='version', version=f'maxiprot {__version__}'
    )
    sp_filter.set_defaults(func=filter.run)

    sp_extract = subparsers.add_parser(
        'extract',
        help='Extract protein/CDS/gene sequences from maxiprot GFF3.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    extract.add_arguments(sp_extract)
    sp_extract.add_argument(
        '--version', action='version', version=f'maxiprot {__version__}'
    )
    sp_extract.set_defaults(func=extract.run)

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Entrypoint: parse arguments and dispatch to the selected subcommand.

    Parameters
    ----------
    argv : Sequence[str] or None, optional
        Argument vector to parse instead of ``sys.argv[1:]``.

    Returns
    -------
    int
        Exit status code.
    """
    if argv is None:
        argv = sys.argv[1:]
    parser = build_parser()
    args = parser.parse_args(argv)
    return guard_broken_pipe(lambda: int(args.func(args)))


if __name__ == '__main__':
    raise SystemExit(main())
