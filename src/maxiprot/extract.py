#!/usr/bin/env python3
"""Extract sequences from GFF3 annotations produced by `maxiprot`.

Supports three extraction modes:
    - protein : translate spliced CDS to amino acids (default)
    - cds     : spliced CDS nucleotide sequence (5'->3' transcript orientation)
    - gene    : genomic span of the gene feature (strand-corrected; includes introns)

Reads a GFF3 (from file or stdin), fetches sequence from a genome FASTA
(plain or bgzip-compressed) using pyfaidx, and writes FASTA to stdout by default
(or to a file via --out-faa).

Key behaviors
-------------
- If no .fai index is present for the genome, one is created automatically.
- For CDS/protein: CDS parts are ordered by transcript orientation and the
  first part is trimmed by its `phase` (per GFF3 spec) before concatenation.
- Warn if concatenated CDS length is not divisible by 3 (remainder is dropped).
- Warn if a downstream CDS part's phase disagrees with the cumulative length
  (a sign of a frameshifted alignment).
- Translation and reverse complement use Biopython, so all NCBI translation
  tables are available and IUPAC ambiguity codes are preserved.
- `--exclude-pseudogenes` (protein mode only) drops translations with internal '*'.

Examples
--------
    # Extract proteins from a maxiprot GFF3 and write FASTA to stdout
    maxiprot filter ... | maxiprot extract - -g genome.fa > proteins.faa

    # Extract CDS nucleotides instead of proteins
    maxiprot extract annotations.gff3 -g genome.fa --extract cds > cds.fna

    # Extract gene genomic sequences (strand-corrected)
    maxiprot extract annotations.gff3 -g genome.fa --extract gene > genes.fna
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
import logging
from pathlib import Path
from typing import Dict, List, Optional, Sequence, TextIO, Tuple

from Bio.Data import CodonTable
from Bio.Seq import Seq
from pyfaidx import Fasta  # type: ignore

from maxiprot.gffio import GffFeature, iter_gff3
from maxiprot.ioutils import guard_broken_pipe, open_output, read_lines
from maxiprot.logs import LOG_LEVELS, init_logging

logger = logging.getLogger(__name__)

#: NCBI translation table IDs supported by Biopython.
TRANSL_TABLE_IDS = sorted(CodonTable.unambiguous_dna_by_id)


# ------------------------------- Data model ---------------------------------


@dataclass
class CdsPart:
    """
    A single CDS feature line.
    """

    seqid: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    strand: str  # '+' or '-'
    phase: int  # 0,1,2
    parent: str  # mRNA ID


@dataclass
class Mrna:
    """
    An mRNA with its CDS parts.
    """

    mrna_id: str
    seqid: str
    strand: str
    start: int
    end: int
    attrs: Dict[str, str]
    cds_parts: List[CdsPart] = field(default_factory=list)
    gene_parent: Optional[str] = None  # gene ID if available


@dataclass
class Gene:
    """
    A gene feature with genomic span.
    """

    gene_id: str
    seqid: str
    strand: str
    start: int
    end: int
    attrs: Dict[str, str]
    mrnas: List[str] = field(default_factory=list)


# ------------------------------- GFF parsing --------------------------------


def _parse_phase(phase_s: str, seqid: str, start: int, end: int) -> int:
    """Parse a CDS phase column, warning when it is missing/malformed."""
    if phase_s in {'0', '1', '2'}:
        return int(phase_s)
    logger.warning(
        'CDS at %s:%d-%d has missing/invalid phase %r; assuming 0',
        seqid,
        start,
        end,
        phase_s,
    )
    return 0


def parse_gff(
    gff_lines: Iterable[str],
) -> Tuple[Dict[str, Mrna], Dict[str, int], Dict[str, Gene], Dict[str, int]]:
    """
    Parse gene, mRNA, and CDS features from a GFF3 stream.

    Feature order does not matter: CDS lines may precede their mRNA, and
    mRNA lines may precede their gene.

    Parameters
    ----------
    gff_lines : Iterable[str]
        Lines of a GFF3 file.

    Returns
    -------
    tuple
        (mRNAs, contig_mrna_counts, genes, contig_gene_counts)
    """
    mrnas: Dict[str, Mrna] = {}
    genes: Dict[str, Gene] = {}

    for feat in iter_gff3(gff_lines):
        if not isinstance(feat, GffFeature):
            continue

        if feat.type == 'gene':
            gene_id = feat.attrs.get('ID')
            if not gene_id:
                logger.warning(
                    'gene feature missing ID on %s:%s:%s-%s',
                    feat.seqid,
                    feat.strand,
                    feat.start,
                    feat.end,
                )
                continue
            if gene_id in genes:
                logger.warning('Duplicate gene ID %s; keeping the last one', gene_id)
            genes[gene_id] = Gene(
                gene_id=gene_id,
                seqid=feat.seqid,
                strand=feat.strand,
                start=feat.start,
                end=feat.end,
                attrs=feat.attrs,
            )

        elif feat.type == 'mRNA':
            mrna_id = feat.attrs.get('ID')
            if not mrna_id:
                logger.warning(
                    'mRNA feature missing ID on %s:%s:%s-%s',
                    feat.seqid,
                    feat.strand,
                    feat.start,
                    feat.end,
                )
                continue
            parent_gene = feat.attrs.get('Parent')
            shell = mrnas.get(mrna_id)
            if shell is not None and shell.attrs == {}:
                # A shell created from CDS lines seen earlier: merge, keeping
                # the collected CDS parts.
                shell.seqid = feat.seqid
                shell.strand = feat.strand
                shell.start = feat.start
                shell.end = feat.end
                shell.attrs = feat.attrs
                shell.gene_parent = parent_gene
            else:
                if shell is not None:
                    logger.warning(
                        'Duplicate mRNA ID %s; keeping the last one', mrna_id
                    )
                mrnas[mrna_id] = Mrna(
                    mrna_id=mrna_id,
                    seqid=feat.seqid,
                    strand=feat.strand,
                    start=feat.start,
                    end=feat.end,
                    attrs=feat.attrs,
                    gene_parent=parent_gene,
                )

        elif feat.type == 'CDS':
            parent = feat.attrs.get('Parent')
            if not parent:
                continue
            part = CdsPart(
                seqid=feat.seqid,
                start=feat.start,
                end=feat.end,
                strand=feat.strand,
                phase=_parse_phase(feat.phase, feat.seqid, feat.start, feat.end),
                parent=parent,
            )
            if parent not in mrnas:
                # mRNA not seen yet: create a shell (merged if/when the mRNA
                # line arrives; attrs stays empty to mark it as a shell).
                mrnas[parent] = Mrna(
                    mrna_id=parent,
                    seqid=feat.seqid,
                    strand=feat.strand,
                    start=feat.start,
                    end=feat.end,
                    attrs={},
                    cds_parts=[part],
                )
            else:
                mrnas[parent].cds_parts.append(part)

    # Link genes to mRNAs regardless of feature order.
    for mrna_id, m in mrnas.items():
        if m.gene_parent and m.gene_parent in genes:
            gene = genes[m.gene_parent]
            if mrna_id not in gene.mrnas:
                gene.mrnas.append(mrna_id)

    contig_mrna_counts: Dict[str, int] = {}
    for m in mrnas.values():
        contig_mrna_counts[m.seqid] = contig_mrna_counts.get(m.seqid, 0) + 1
    contig_gene_counts: Dict[str, int] = {}
    for g in genes.values():
        contig_gene_counts[g.seqid] = contig_gene_counts.get(g.seqid, 0) + 1

    return mrnas, contig_mrna_counts, genes, contig_gene_counts


# ------------------------------ Sequence utils ------------------------------


def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement, preserving IUPAC ambiguity codes.

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    str
        Reverse-complemented sequence.
    """
    return str(Seq(seq).reverse_complement())


def translate_nt(seq_nt: str, table_id: int) -> str:
    """
    Translate a nucleotide CDS into amino acids using Biopython.

    Codons that cannot be resolved (e.g. containing 'N') yield 'X'.
    Trailing bases that don't complete a codon are ignored.

    Parameters
    ----------
    seq_nt : str
        Nucleotide CDS sequence (5'->3').
    table_id : int
        NCBI translation table ID.

    Returns
    -------
    str
        Amino-acid sequence (single-letter), '*' for stops.

    Raises
    ------
    ValueError
        If the translation table ID is unknown.
    """
    if table_id not in CodonTable.unambiguous_dna_by_id:
        raise ValueError(
            f'Unsupported translation table: {table_id}. '
            f'Supported: {", ".join(map(str, TRANSL_TABLE_IDS))}'
        )
    seq = seq_nt.upper().replace('U', 'T')
    seq = seq[: len(seq) - (len(seq) % 3)]
    try:
        return str(Seq(seq).translate(table=table_id))
    except CodonTable.TranslationError:
        # Fall back to per-codon translation, substituting 'X' for codons
        # Biopython cannot resolve.
        aa_chars: List[str] = []
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            try:
                aa_chars.append(str(Seq(codon).translate(table=table_id)))
            except CodonTable.TranslationError:
                aa_chars.append('X')
        return ''.join(aa_chars)


def fetch_interval(fa: Fasta, seqid: str, start_1based: int, end_1based: int) -> str:
    """
    Fetch genomic interval in FASTA (1-based inclusive coordinates).

    Raises
    ------
    KeyError
        If ``seqid`` is not present in the FASTA.
    """
    return str(fa[seqid][start_1based - 1 : end_1based]).upper()


def _expected_phase(cumulative_len: int) -> int:
    """Phase implied for the next CDS part after ``cumulative_len`` bases."""
    return (3 - (cumulative_len % 3)) % 3


def fetch_cds_sequence(fa: Fasta, cds_parts: List[CdsPart], strand: str) -> str:
    """
    Fetch and concatenate CDS nucleotide in transcript orientation.

    The GFF3 ``phase`` is applied **only to the first CDS part** in transcript
    order (5'->3' of the mRNA). The phase of each subsequent part is verified
    against the cumulative length; a mismatch (typical of frameshifted
    miniprot alignments) triggers a warning, since the translation downstream
    of the mismatch will be out of frame.

    Transcript order:
      - '+' strand: CDS parts sorted by ascending genomic start.
      - '-' strand: CDS parts sorted by descending genomic end; each piece is
        reverse-complemented to maintain 5'->3' transcript orientation.

    Parameters
    ----------
    fa : Fasta
        Open pyfaidx Fasta object.
    cds_parts : list of CdsPart
        CDS features for a single mRNA.
    strand : str
        '+' or '-'.

    Returns
    -------
    str
        Concatenated CDS nucleotide sequence (5'->3', transcript orientation).

    Raises
    ------
    ValueError
        If ``strand`` is not '+' or '-'.
    """
    if strand not in {'+', '-'}:
        raise ValueError(f'Invalid strand: {strand}')

    if strand == '-':
        parts = sorted(cds_parts, key=lambda p: p.end, reverse=True)
    else:
        parts = sorted(cds_parts, key=lambda p: p.start)

    chunks: List[str] = []
    cumulative = 0
    for idx, part in enumerate(parts):
        s = fetch_interval(fa, part.seqid, part.start, part.end)
        if strand == '-':
            s = reverse_complement(s)

        if idx == 0:
            # Apply phase ONLY to the first part in transcript order.
            if part.phase in (1, 2):
                if len(s) >= part.phase:
                    s = s[part.phase :]
                else:
                    logger.warning(
                        'First CDS part shorter than its phase; trimming to empty '
                        '(seqid=%s start=%d end=%d phase=%d)',
                        part.seqid,
                        part.start,
                        part.end,
                        part.phase,
                    )
                    s = ''
        else:
            expected = _expected_phase(cumulative)
            if part.phase != expected:
                logger.warning(
                    'CDS part at %s:%d-%d has phase %d but cumulative length '
                    'implies %d; the alignment is likely frameshifted and the '
                    'translation downstream may be incorrect.',
                    part.seqid,
                    part.start,
                    part.end,
                    part.phase,
                    expected,
                )

        cumulative += len(s)
        chunks.append(s)

    return ''.join(chunks)


def count_non_acgt(seq: str) -> int:
    """
    Count non-ACGT letters (after uppercasing and U->T).
    """
    seq = seq.upper().replace('U', 'T')
    return sum(1 for ch in seq if ch not in {'A', 'C', 'G', 'T'})


# --------------------------------- FASTA I/O --------------------------------


def write_fasta(
    records: Iterable[Tuple[str, str]], handle: TextIO, width: int = 60
) -> None:
    """
    Write FASTA records to a text handle.

    Parameters
    ----------
    records : iterable of (header, sequence)
        Header lines must NOT include the leading '>'.
    handle : TextIO
        Output stream.
    width : int
        Line wrap width for sequences (default 60).
    """
    for header, seq in records:
        handle.write(f'>{header}\n')
        for i in range(0, len(seq), width):
            handle.write(seq[i : i + width] + '\n')


def build_header(record_id: str, location: str, extra: Dict[str, object]) -> str:
    """
    Build a FASTA header from an ID, a location string, and key=value fields.

    Parameters
    ----------
    record_id : str
        Feature ID (first token of the header).
    location : str
        ``seqid:start-end(strand)`` location string.
    extra : dict
        Additional ``key=value`` fields; falsy values are skipped.

    Returns
    -------
    str
        Space-separated header (without the leading '>').
    """
    parts = [record_id, location]
    for key, value in extra.items():
        if value:
            parts.append(f'{key}={value}')
    return ' '.join(parts)


# ---------------------------------- CLI/Main --------------------------------


def add_arguments(ap: argparse.ArgumentParser) -> None:
    """
    Add ``maxiprot extract`` arguments to a parser.

    Parameters
    ----------
    ap : argparse.ArgumentParser
        Parser (top-level subparser or standalone) to populate.
    """
    ap.add_argument('gff', nargs='?', default='-', help="GFF3 path or '-' for stdin.")
    ap.add_argument(
        '-g',
        '--genome',
        required=True,
        help='Path to genome FASTA (plain or bgzip-compressed).',
    )
    ap.add_argument('--out-faa', default='-', help="Output FASTA path ('-' = stdout).")
    ap.add_argument(
        '--extract',
        choices=['protein', 'cds', 'gene'],
        default='protein',
        help='What to extract: translated protein, spliced CDS, or gene genomic span.',
    )
    ap.add_argument(
        '--transl-table',
        type=int,
        choices=TRANSL_TABLE_IDS,
        default=1,
        metavar='ID',
        help=(
            'NCBI translation table ID. '
            f'Supported: {", ".join(map(str, TRANSL_TABLE_IDS))}.'
        ),
    )
    ap.add_argument(
        '--max-annos-per-contig',
        type=int,
        default=0,
        help=(
            'Max number of annotations allowed per contig; '
            'if a contig exceeds this number, all sequences for that contig are '
            'skipped. Unit depends on --extract: mRNA for protein/cds; gene for '
            'gene mode. 0 or negative = unlimited.'
        ),
    )
    ap.add_argument(
        '--exclude-pseudogenes',
        action='store_true',
        help="(protein mode only) Exclude translations with internal '*'.",
    )
    ap.add_argument(
        '--log-level', default='INFO', choices=LOG_LEVELS, help='Logging level.'
    )


def build_arg_parser() -> argparse.ArgumentParser:
    """
    Build the standalone ``maxiprot extract`` parser.

    Returns
    -------
    argparse.ArgumentParser
        A configured argument parser instance.
    """
    ap = argparse.ArgumentParser(
        prog='maxiprot extract',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Extract protein/CDS/gene sequences from maxiprot GFF3.',
    )
    add_arguments(ap)
    return ap


@dataclass
class _Stats:
    """Counters shared by the record generators."""

    n_out: int = 0
    n_skipped_contig: int = 0
    n_skipped_pseudo: int = 0
    n_fetch_errors: int = 0


def _iter_mrna_records(
    mrnas: Dict[str, Mrna],
    fa: Fasta,
    mode: str,
    transl_table: int,
    exclude_pseudogenes: bool,
    disallowed_contigs: set[str],
    stats: _Stats,
) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence) records for protein/cds modes."""
    for mrna_id, m in mrnas.items():
        if m.seqid in disallowed_contigs:
            stats.n_skipped_contig += 1
            continue
        if not m.cds_parts:
            logger.warning('mRNA %s has no CDS parts; skipped', mrna_id)
            continue

        try:
            cds_nt = fetch_cds_sequence(fa, m.cds_parts, m.strand)
        except KeyError:
            logger.error(
                'Contig %s (mRNA %s) not found in genome FASTA; skipped',
                m.seqid,
                mrna_id,
            )
            stats.n_fetch_errors += 1
            continue

        non_acgt = count_non_acgt(cds_nt)
        if non_acgt > 0:
            logger.warning('CDS for %s contains %d non-ACGT bases', mrna_id, non_acgt)

        if len(cds_nt) % 3 != 0:
            logger.warning(
                'CDS length not divisible by 3 for %s (len=%d) - truncating remainder',
                mrna_id,
                len(cds_nt),
            )
            cds_nt = cds_nt[: len(cds_nt) - (len(cds_nt) % 3)]

        location = f'{m.seqid}:{m.start}-{m.end}({m.strand})'
        if mode == 'cds':
            header = build_header(
                mrna_id,
                location,
                {'feature': 'CDS', 'target': m.attrs.get('Target', '')},
            )
            stats.n_out += 1
            yield header, cds_nt
        else:  # protein
            prot = translate_nt(cds_nt, transl_table)
            internal_stop = '*' in prot[:-1] if prot else False
            if internal_stop and exclude_pseudogenes:
                stats.n_skipped_pseudo += 1
                continue
            header = build_header(
                mrna_id,
                location,
                {
                    'table': transl_table,
                    'feature': 'protein',
                    'target': m.attrs.get('Target', ''),
                    'pseudogene': 1 if internal_stop else None,
                },
            )
            stats.n_out += 1
            yield header, prot


def _iter_gene_records(
    genes: Dict[str, Gene],
    fa: Fasta,
    disallowed_contigs: set[str],
    stats: _Stats,
) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence) records for gene mode."""
    for gene_id, g in genes.items():
        if g.seqid in disallowed_contigs:
            stats.n_skipped_contig += 1
            continue
        try:
            gene_nt = fetch_interval(fa, g.seqid, g.start, g.end)
        except KeyError:
            logger.error(
                'Contig %s (gene %s) not found in genome FASTA; skipped',
                g.seqid,
                gene_id,
            )
            stats.n_fetch_errors += 1
            continue
        if g.strand == '-':
            gene_nt = reverse_complement(gene_nt)

        non_acgt = count_non_acgt(gene_nt)
        if non_acgt > 0:
            logger.warning(
                'Gene %s sequence contains %d non-ACGT bases', gene_id, non_acgt
            )

        header = build_header(
            gene_id,
            f'{g.seqid}:{g.start}-{g.end}({g.strand})',
            {
                'feature': 'gene',
                'name': g.attrs.get('Name'),
                'mrnas': len(g.mrnas) if g.mrnas else None,
            },
        )
        stats.n_out += 1
        yield header, gene_nt


def run(args: argparse.Namespace) -> int:
    """
    Run the extract pipeline for parsed arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments (see :func:`add_arguments`).

    Returns
    -------
    int
        Exit status: 0 success, 1 partial data errors, 2 fatal/usage error.
    """
    init_logging(loglevel=args.log_level)
    logger.info('Reading GFF: %s', args.gff)
    logger.info('Genome FASTA: %s', args.genome)
    logger.info('Extraction mode: %s', args.extract)

    mrnas, mrna_counts, genes, gene_counts = parse_gff(read_lines(args.gff))
    logger.info('Parsed %d mRNA(s) across %d contig(s)', len(mrnas), len(mrna_counts))
    logger.info('Parsed %d gene(s) across %d contig(s)', len(genes), len(gene_counts))

    disallowed_contigs: set[str] = set()
    if args.max_annos_per_contig and args.max_annos_per_contig > 0:
        counts = gene_counts if args.extract == 'gene' else mrna_counts
        for seqid, n in counts.items():
            if n > args.max_annos_per_contig:
                disallowed_contigs.add(seqid)
        if disallowed_contigs:
            logger.warning(
                'Skipping all sequences on %d contig(s) due to --max-annos-per-contig=%d: %s',
                len(disallowed_contigs),
                args.max_annos_per_contig,
                ', '.join(sorted(disallowed_contigs)),
            )

    if args.extract in {'cds', 'gene'} and args.exclude_pseudogenes:
        logger.warning(
            '--exclude-pseudogenes applies to protein mode only; ignoring for %s mode',
            args.extract,
        )
    if args.extract in {'cds', 'gene'} and args.transl_table != 1:
        logger.warning(
            '--transl-table applies to protein mode only; ignoring for %s mode',
            args.extract,
        )

    fai_path = Path(args.genome + '.fai')
    if not fai_path.exists():
        logger.info(
            'No FASTA index found (%s); creating with pyfaidx...', fai_path.name
        )
    try:
        fa = Fasta(args.genome, as_raw=True, sequence_always_upper=True)
    except (OSError, ValueError) as e:
        logger.error('Cannot open genome FASTA %s: %s', args.genome, e)
        return 2

    stats = _Stats()
    out_handle, close_out = open_output(args.out_faa)
    try:
        if args.extract in {'protein', 'cds'}:
            records: Iterator[Tuple[str, str]] = _iter_mrna_records(
                mrnas,
                fa,
                args.extract,
                args.transl_table,
                args.exclude_pseudogenes,
                disallowed_contigs,
                stats,
            )
        else:
            records = _iter_gene_records(genes, fa, disallowed_contigs, stats)
        write_fasta(records, out_handle)
    finally:
        if close_out:
            out_handle.close()
        fa.close()

    logger.info('Wrote %d record(s)', stats.n_out)
    if stats.n_skipped_contig:
        logger.info(
            'Skipped %d record(s) due to contig annotation limit',
            stats.n_skipped_contig,
        )
    if stats.n_skipped_pseudo:
        logger.info(
            'Excluded %d pseudogene translation(s) with internal stops',
            stats.n_skipped_pseudo,
        )
    if stats.n_fetch_errors:
        logger.error(
            '%d record(s) could not be fetched from the genome FASTA',
            stats.n_fetch_errors,
        )
        return 1

    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Run the ``maxiprot extract`` command-line interface.

    Parameters
    ----------
    argv : Sequence[str] or None, optional
        Argument vector to parse instead of ``sys.argv``. Useful for testing.

    Returns
    -------
    int
        Exit status code (0 on success).
    """
    args = build_arg_parser().parse_args(argv)
    return run(args)


if __name__ == '__main__':
    raise SystemExit(guard_broken_pipe(main))
