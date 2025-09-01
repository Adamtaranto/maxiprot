#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract sequences from GFF3 annotations produced by `maxiprot`.

Supports three extraction modes:
    - protein : translate spliced CDS to amino acids (default)
    - cds     : spliced CDS nucleotide sequence (5'→3' transcript orientation)
    - gene    : genomic span of the gene feature (strand-corrected; includes introns)

Reads a GFF3 (from file or stdin), fetches sequence from a genome FASTA
(plain or bgzip-compressed) using pyfaidx, and writes FASTA to stdout by default
(or to a file via --out-faa).

Key behaviors
-------------
- If no .fai index is present for the genome, one is created automatically.
- For CDS/protein: CDS parts are ordered by transcript orientation and trimmed
  by `phase` (per GFF3 spec) before concatenation.
- Warn if concatenated CDS length is not divisible by 3 (remainder is dropped).
- Warn if sequences contain non-ATGC letters (translation uses 'X' for ambiguous codons).
- `--exclude-pseudogenes` (protein mode only) drops translations with internal '*'.
  In CDS/gene modes this option is ignored (with a one-time warning).

Examples
--------
    # Extract proteins from a maxiprot GFF3 and write FASTA to stdout
    ./maxiprot ... | ./maxiprot_extract_proteins.py -g genome.fa > proteins.faa

    # Extract CDS nucleotides instead of proteins
    ./maxiprot_extract_proteins.py annotations.gff3 -g genome.fa --extract cds > cds.fna

    # Extract gene genomic sequences (strand-corrected)
    ./maxiprot_extract_proteins.py annotations.gff3 -g genome.fa --extract gene > genes.fna

    # Exclude pseudogenes (protein mode only) and use bacterial translation table (11)
    ./maxiprot_extract_proteins.py annotations.gff3 -g genome.fa \
        --exclude-pseudogenes --transl-table 11 > proteins.faa
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import logging
from pathlib import Path
import sys
from typing import Dict, Iterable, List, Optional, Sequence, TextIO, Tuple

from pyfaidx import Fasta  # type: ignore

from maxiprot._version import __version__

# --------------------------- Translation tables -----------------------------


def _codon_table_standard() -> Dict[str, str]:
    """NCBI translation table 1 (Standard Genetic Code)."""
    # fmt: off
    table = {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
        "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
        "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
        "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
        "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
        "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
        "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
        "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
        "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
        "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
        "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
        "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
        "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
        "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
        "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
    }
    # fmt: on
    return table


def _codon_table_bacterial() -> Dict[str, str]:
    """NCBI translation table 11 (Bacterial, Archaeal and Plant Plastid)."""
    return _codon_table_standard()


def _codon_table_mold_protozoan_mito() -> Dict[str, str]:
    """NCBI translation table 4 (Mold, Protozoan, etc. Mitochondrial Code)."""
    tbl = _codon_table_standard().copy()
    tbl['TGA'] = 'W'
    return tbl


def _codon_table_vert_mito() -> Dict[str, str]:
    """NCBI translation table 2 (Vertebrate Mitochondrial Code)."""
    tbl = _codon_table_standard().copy()
    tbl['AGA'] = '*'
    tbl['AGG'] = '*'
    tbl['TGA'] = 'W'
    return tbl


def get_codon_table(table_id: int) -> Dict[str, str]:
    """Return a codon→AA map for the requested NCBI translation table.

    Parameters
    ----------
    table_id : int
        NCBI genetic code table ID (1, 2, 4, 11).

    Returns
    -------
    dict
        Mapping from triplet codon (T) to single-letter amino acid.

    Raises
    ------
    ValueError
        If the table is not supported.
    """
    if table_id == 1:
        return _codon_table_standard()
    if table_id == 11:
        return _codon_table_bacterial()
    if table_id == 4:
        return _codon_table_mold_protozoan_mito()
    if table_id == 2:
        return _codon_table_vert_mito()
    raise ValueError(
        f'Unsupported translation table: {table_id}. Supported: 1, 2, 4, 11'
    )


def translate_nt(seq_nt: str, table_id: int) -> str:
    """Translate nucleotide sequence into amino acids.

    Any codon containing non-ACGT letters yields 'X'.
    Trailing bases that don't complete a codon are ignored.

    Parameters
    ----------
    seq_nt : str
        Nucleotide CDS sequence (5'→3').
    table_id : int
        NCBI translation table ID.

    Returns
    -------
    str
        Amino-acid sequence (single-letter), '*' for stops.
    """
    seq = seq_nt.upper().replace('U', 'T')
    tbl = get_codon_table(table_id)
    aa_chars: List[str] = []
    for i in range(0, len(seq) - (len(seq) % 3), 3):
        codon = seq[i : i + 3]
        if set(codon) <= {'A', 'C', 'G', 'T'}:
            aa_chars.append(tbl.get(codon, 'X'))
        else:
            aa_chars.append('X')
    return ''.join(aa_chars)


# ------------------------------- GFF parsing --------------------------------


@dataclass
class CdsPart:
    """A single CDS feature line."""

    seqid: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    strand: str  # '+' or '-'
    phase: int  # 0,1,2
    parent: str  # mRNA ID


@dataclass
class Mrna:
    """An mRNA with its CDS parts."""

    mrna_id: str
    seqid: str
    strand: str
    start: int
    end: int
    attrs: Dict[str, str]
    cds_parts: List[CdsPart]
    gene_parent: Optional[str] = None  # gene ID if available


@dataclass
class Gene:
    """A gene feature with genomic span."""

    gene_id: str
    seqid: str
    strand: str
    start: int
    end: int
    attrs: Dict[str, str]
    mrnas: List[str]


def parse_gff(
    gff_lines: Iterable[str],
) -> Tuple[Dict[str, Mrna], Dict[str, int], Dict[str, Gene], Dict[str, int]]:
    """Parse gene, mRNA, and CDS features from a GFF3 stream.

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
    contig_mrna_counts: Dict[str, int] = {}
    contig_gene_counts: Dict[str, int] = {}

    for ln in gff_lines:
        if not ln or ln.startswith('#'):
            continue
        cols = ln.rstrip('\n').split('\t')
        if len(cols) < 9:
            continue
        seqid, source, ftype, start_s, end_s, score, strand, phase_s, attrs_s = cols[:9]
        start = int(start_s)
        end = int(end_s)
        attrs: Dict[str, str] = {}
        for kv in attrs_s.split(';'):
            if '=' in kv:
                k, v = kv.split('=', 1)
                attrs[k] = v

        if ftype == 'gene':
            gene_id = attrs.get('ID')
            if not gene_id:
                logging.warning(
                    'gene feature missing ID on %s:%s:%s-%s', seqid, strand, start, end
                )
                continue
            genes[gene_id] = Gene(
                gene_id=gene_id,
                seqid=seqid,
                strand=strand,
                start=start,
                end=end,
                attrs=attrs,
                mrnas=[],
            )
            contig_gene_counts[seqid] = contig_gene_counts.get(seqid, 0) + 1

        elif ftype == 'mRNA':
            mrna_id = attrs.get('ID')
            if not mrna_id:
                logging.warning(
                    'mRNA feature missing ID on %s:%s:%s-%s', seqid, strand, start, end
                )
                continue
            parent_gene = attrs.get('Parent')
            mrnas[mrna_id] = Mrna(
                mrna_id=mrna_id,
                seqid=seqid,
                strand=strand,
                start=start,
                end=end,
                attrs=attrs,
                cds_parts=[],
                gene_parent=parent_gene,
            )
            contig_mrna_counts[seqid] = contig_mrna_counts.get(seqid, 0) + 1
            if parent_gene and parent_gene in genes:
                genes[parent_gene].mrnas.append(mrna_id)

        elif ftype == 'CDS':
            parent = attrs.get('Parent')
            if not parent:
                continue
            try:
                phase = int(phase_s) if phase_s in {'0', '1', '2'} else 0
            except Exception:
                phase = 0
            part = CdsPart(
                seqid=seqid,
                start=start,
                end=end,
                strand=strand,
                phase=phase,
                parent=parent,
            )
            if parent not in mrnas:
                # If mRNA not seen yet, create a shell (gene parent unknown)
                mrnas[parent] = Mrna(
                    mrna_id=parent,
                    seqid=seqid,
                    strand=strand,
                    start=start,
                    end=end,
                    attrs={},
                    cds_parts=[part],
                    gene_parent=None,
                )
                contig_mrna_counts[seqid] = contig_mrna_counts.get(seqid, 0) + 1
            else:
                mrnas[parent].cds_parts.append(part)

    return mrnas, contig_mrna_counts, genes, contig_gene_counts


# ------------------------------ Sequence utils ------------------------------


def reverse_complement(seq: str) -> str:
    """Return reverse complement for A/C/G/T/N sequences (others -> N)."""
    comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    seq_rc = seq.translate(comp)
    return ''.join(ch if ch in 'ACGTNacgtn' else 'N' for ch in seq_rc[::-1])


def fetch_interval(fa: Fasta, seqid: str, start_1based: int, end_1based: int) -> str:
    """Fetch genomic interval in FASTA (1-based inclusive coordinates)."""
    return str(fa[seqid][start_1based - 1 : end_1based]).upper()


def fetch_cds_sequence(fa: Fasta, cds_parts: List[CdsPart], strand: str) -> str:
    """Fetch and concatenate CDS nucleotide in transcript orientation.

    GFF3 uses 1-based inclusive; pyfaidx uses 0-based end-exclusive.
    Phase handling:
      - For '+' strand: remove `phase` from the LEFT of each part.
      - For '-' strand: reverse-complement each piece first, then remove its `phase` from the LEFT.

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
        Concatenated CDS nucleotide sequence (5'→3', transcript orientation).
    """
    if strand not in {'+', '-'}:
        raise ValueError(f'Invalid strand: {strand}')

    parts = sorted(cds_parts, key=lambda p: p.start, reverse=(strand == '-'))
    chunks: List[str] = []
    for part in parts:
        s = fetch_interval(fa, part.seqid, part.start, part.end)
        if strand == '-':
            s = reverse_complement(s)
        if part.phase in (1, 2):
            if len(s) >= part.phase:
                s = s[part.phase :]
            else:
                logging.warning(
                    'CDS part shorter than phase (seqid=%s start=%d end=%d phase=%d)',
                    part.seqid,
                    part.start,
                    part.end,
                    part.phase,
                )
                s = ''
        chunks.append(s)
    return ''.join(chunks)


def count_non_acgt(seq: str) -> int:
    """Count non-ACGT letters (after uppercasing and U->T)."""
    seq = seq.upper().replace('U', 'T')
    return sum(1 for ch in seq if ch not in {'A', 'C', 'G', 'T'})


# --------------------------------- FASTA I/O --------------------------------


def write_fasta(
    records: Iterable[Tuple[str, str]], handle: TextIO, width: int = 60
) -> None:
    """Write FASTA records to a text handle.

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


# ---------------------------------- CLI/Main --------------------------------


def configure_logging(level: str = 'INFO') -> None:
    """Configure root logger to stderr with a concise format."""
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format='%(asctime)s [%(levelname)s] %(message)s',
        stream=sys.stderr,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    """Build CLI parser."""
    ap = argparse.ArgumentParser(
        description=f'Extract sequences from maxiprot GFF3 (v{__version__})'
    )
    ap.add_argument(
        'gff', nargs='?', default='-', help="GFF3 path or '-' for stdin (default: '-')"
    )
    ap.add_argument(
        '-g',
        '--genome',
        required=True,
        help='Path to genome FASTA (plain or bgzip-compressed)',
    )
    ap.add_argument(
        '--out-faa', default='-', help="Output FASTA path (default: '-' = stdout)"
    )
    ap.add_argument(
        '--extract',
        choices=['protein', 'cds', 'gene'],
        default='protein',
        help='What to extract: translated protein, spliced CDS, or gene genomic span (default: protein)',
    )
    ap.add_argument(
        '--transl-table',
        type=int,
        default=1,
        help='NCBI translation table ID (default: 1)',
    )
    ap.add_argument(
        '--max-annos-per-contig',
        type=int,
        default=0,
        help=(
            'Max number of annotations allowed per contig; '
            'if a contig exceeds this number, all sequences for that contig are skipped. '
            'Unit depends on --extract: mRNA for protein/cds; gene for gene mode. '
            '0 or negative = unlimited (default: 0).'
        ),
    )
    ap.add_argument(
        '--exclude-pseudogenes',
        action='store_true',
        help="(protein mode only) Exclude translations with internal '*' (premature stop).",
    )
    ap.add_argument(
        '--log-level',
        default='INFO',
        help='Logging level (DEBUG, INFO, WARNING, ERROR)',
    )
    return ap


def open_gff_lines(path: str) -> Iterable[str]:
    """Yield lines from a GFF path or stdin."""
    if path == '-' or path is None:
        for ln in sys.stdin:
            yield ln
    else:
        with open(path, 'r', encoding='utf-8') as fh:
            for ln in fh:
                yield ln


def _open_out_handle(path: str) -> Tuple[TextIO, bool]:
    """Open output handle (stdout if '-' or empty)."""
    if path in (None, '', '-'):
        return sys.stdout, False
    return open(path, 'w', encoding='utf-8'), True


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entrypoint."""
    ap = build_arg_parser()
    args = ap.parse_args(argv)

    configure_logging(args.log_level)
    logging.info('maxiprot_extract_proteins version %s', __version__)
    logging.info('Reading GFF: %s', args.gff)
    logging.info('Genome FASTA: %s', args.genome)
    logging.info('Extraction mode: %s', args.extract)

    # Parse GFF and gather features
    mrnas, mrna_counts, genes, gene_counts = parse_gff(open_gff_lines(args.gff))
    logging.info('Parsed %d mRNA(s) across %d contig(s)', len(mrnas), len(mrna_counts))
    logging.info('Parsed %d gene(s) across %d contig(s)', len(genes), len(gene_counts))

    # Contig-level gating by mode
    disallowed_contigs: set[str] = set()
    if args.max_annos_per_contig and args.max_annos_per_contig > 0:
        counts = gene_counts if args.extract == 'gene' else mrna_counts
        for seqid, n in counts.items():
            if n > args.max_annos_per_contig:
                disallowed_contigs.add(seqid)
        if disallowed_contigs:
            logging.warning(
                'Skipping all sequences on %d contig(s) due to --max-annos-per-contig=%d: %s',
                len(disallowed_contigs),
                args.max_annos_per_contig,
                ', '.join(sorted(disallowed_contigs)),
            )

    # Warn about options unused in certain modes
    if args.extract in {'cds', 'gene'} and args.exclude_pseudogenes:
        logging.warning(
            '--exclude-pseudogenes applies to protein mode only; ignoring for %s mode',
            args.extract,
        )
    if args.extract in {'cds', 'gene'} and args.transl_table != 1:
        logging.warning(
            '--transl-table applies to protein mode only; ignoring for %s mode',
            args.extract,
        )

    # Open genome with pyfaidx (auto-creates .fai if missing)
    fai_path = Path(args.genome + '.fai')
    if not fai_path.exists():
        logging.info(
            'No FASTA index found (%s); creating with pyfaidx...', fai_path.name
        )
    fa = Fasta(
        args.genome, as_raw=True, sequence_always_upper=True
    )  # indexes if needed

    # Prepare output
    out_handle, close_out = _open_out_handle(args.out_faa)

    n_out = 0
    n_skipped_contig = 0
    n_skipped_pseudo = 0

    try:
        # Protein/CDS modes are per-mRNA; gene mode is per-gene
        records: List[Tuple[str, str]] = []

        if args.extract in {'protein', 'cds'}:
            for mrna_id, m in mrnas.items():
                if m.seqid in disallowed_contigs:
                    n_skipped_contig += 1
                    continue
                if not m.cds_parts:
                    logging.warning('mRNA %s has no CDS parts; skipped', mrna_id)
                    continue

                # Fetch spliced CDS (handles phase & strand)
                cds_nt = fetch_cds_sequence(fa, m.cds_parts, m.strand)

                # Sanity checks
                non_acgt = count_non_acgt(cds_nt)
                if non_acgt > 0:
                    logging.warning(
                        'CDS for %s contains %d non-ACGT bases', mrna_id, non_acgt
                    )

                if len(cds_nt) % 3 != 0:
                    logging.warning(
                        'CDS length not divisible by 3 for %s (len=%d) — truncating remainder',
                        mrna_id,
                        len(cds_nt),
                    )
                    cds_nt = cds_nt[: len(cds_nt) - (len(cds_nt) % 3)]

                # Decide output per mode
                if args.extract == 'cds':
                    header_parts = [
                        mrna_id,
                        f'{m.seqid}:{m.start}-{m.end}({m.strand})',
                        'feature=CDS',
                    ]
                    tgt = m.attrs.get('Target', '')
                    if tgt:
                        header_parts.append(f'target={tgt}')
                    records.append((' '.join(header_parts), cds_nt))
                    n_out += 1
                else:  # protein
                    try:
                        prot = translate_nt(cds_nt, args.transl_table)
                    except ValueError as e:
                        logging.error(str(e))
                        fa.close()
                        if close_out:
                            out_handle.close()
                        return 2

                    internal_stop = '*' in prot[:-1] if prot else False
                    if internal_stop and args.exclude_pseudogenes:
                        n_skipped_pseudo += 1
                        continue

                    header_parts = [
                        mrna_id,
                        f'{m.seqid}:{m.start}-{m.end}({m.strand})',
                        f'table={args.transl_table}',
                        'feature=protein',
                    ]
                    tgt = m.attrs.get('Target', '')
                    if tgt:
                        header_parts.append(f'target={tgt}')
                    if internal_stop:
                        header_parts.append('pseudogene=1')
                    records.append((' '.join(header_parts), prot))
                    n_out += 1

        else:  # gene mode
            for gene_id, g in genes.items():
                if g.seqid in disallowed_contigs:
                    n_skipped_contig += 1
                    continue
                # Fetch genomic span (includes introns)
                gene_nt = fetch_interval(fa, g.seqid, g.start, g.end)
                if g.strand == '-':
                    gene_nt = reverse_complement(gene_nt)

                non_acgt = count_non_acgt(gene_nt)
                if non_acgt > 0:
                    logging.warning(
                        'Gene %s sequence contains %d non-ACGT bases', gene_id, non_acgt
                    )

                header_parts = [
                    gene_id,
                    f'{g.seqid}:{g.start}-{g.end}({g.strand})',
                    'feature=gene',
                ]
                name = g.attrs.get('Name')
                if name:
                    header_parts.append(f'name={name}')
                # Optionally, mention #mRNAs
                if g.mrnas:
                    header_parts.append(f'mrnas={len(g.mrnas)}')
                records.append((' '.join(header_parts), gene_nt))
                n_out += 1

        # Write output FASTA
        write_fasta(records, out_handle)

    finally:
        if close_out:
            out_handle.close()
        fa.close()

    logging.info('Wrote %d record(s)', n_out)
    if n_skipped_contig:
        logging.info(
            'Skipped %d record(s) due to contig annotation limit', n_skipped_contig
        )
    if n_skipped_pseudo:
        logging.info(
            'Excluded %d pseudogene translation(s) with internal stops',
            n_skipped_pseudo,
        )

    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        # Allow piping to head/tail without traceback noise
        try:
            sys.stderr.close()
        except Exception:
            pass
        try:
            sys.stdout.close()
        except Exception:
            pass
        raise
