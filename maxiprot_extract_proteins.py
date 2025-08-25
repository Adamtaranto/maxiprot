#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract protein sequences from GFF3 CDS annotations produced by `maxiprot`.

Reads a GFF3 (from file or stdin), fetches CDS segments from a genome FASTA
(plain or bgzip-compressed) using pyfaidx, applies GFF3 phase to maintain frame,
translates with a chosen NCBI translation table, and writes proteins in FASTA
format to stdout (default) or to a file.

Key behaviors
-------------
- If no .fai index is present for the genome, one is created automatically.
- CDS segments are ordered by strand (5'→3' transcript order) and trimmed by
  `phase` (per GFF3 spec) before concatenation.
- Warn if concatenated CDS length is not divisible by 3; the trailing remainder
  is dropped for translation.
- Warn if concatenated CDS contains non-ATGC letters (translation uses 'X' for
  any codon with ambiguous letters).
- If `--exclude-pseudogenes` is set, any translation with premature '*' (i.e.,
  internal stop) is skipped. Otherwise, '*' is kept in output.

Examples
--------
    # Extract proteins from a maxiprot GFF3 and write FASTA to stdout
    ./maxiprot ... | ./maxiprot_extract_proteins.py -g genome.fa > proteins.faa

    # Read GFF from file, write FASTA to file, use bacterial table (11)
    ./maxiprot_extract_proteins.py annotations.gff3 -g genome.fa \
        --transl-table 11 --out-faa proteins.faa

    # Skip contigs with >= 2 annotations (require <=1 per contig)
    ./maxiprot_extract_proteins.py annotations.gff3 -g genome.fa \
        --max-annos-per-contig 1 > proteins.faa
"""
from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, TextIO

from pyfaidx import Fasta  # type: ignore

__version__ = "0.1.0"


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
    # Differs mainly in allowed start codons; coding table is the same as 1
    return _codon_table_standard()


def _codon_table_mold_protozoan_mito() -> Dict[str, str]:
    """NCBI translation table 4 (Mold, Protozoan, etc. Mitochondrial Code)."""
    tbl = _codon_table_standard().copy()
    # TGA encodes W
    tbl["TGA"] = "W"
    return tbl


def _codon_table_vert_mito() -> Dict[str, str]:
    """NCBI translation table 2 (Vertebrate Mitochondrial Code)."""
    tbl = _codon_table_standard().copy()
    # AGA/AGG are STOP; TGA is W
    tbl["AGA"] = "*"
    tbl["AGG"] = "*"
    tbl["TGA"] = "W"
    return tbl


def get_codon_table(table_id: int) -> Dict[str, str]:
    """Return a codon→AA dict for the requested NCBI translation table.

    Parameters
    ----------
    table_id : int
        NCBI genetic code table ID (e.g., 1, 2, 4, 11).

    Returns
    -------
    dict
        Mapping from triplet codon (T) to single-letter amino acid (incl. '*').

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
    raise ValueError(f"Unsupported translation table: {table_id}. Supported: 1, 2, 4, 11")


def translate_nt(seq_nt: str, table_id: int) -> str:
    """Translate nucleotide sequence into amino acids.

    Any codon containing non-ACGT letters yields 'X'.
    Trailing bases that don't complete a codon are ignored.

    Parameters
    ----------
    seq_nt : str
        Nucleotide CDS sequence (5'→3'), uppercase recommended.
    table_id : int
        NCBI translation table ID.

    Returns
    -------
    str
        Amino-acid sequence (single-letter), '*' for stops.
    """
    seq = seq_nt.upper().replace("U", "T")
    tbl = get_codon_table(table_id)
    aa_chars: List[str] = []
    for i in range(0, len(seq) - (len(seq) % 3), 3):
        codon = seq[i : i + 3]
        if set(codon) <= {"A", "C", "G", "T"}:
            aa_chars.append(tbl.get(codon, "X"))
        else:
            aa_chars.append("X")
    return "".join(aa_chars)


# ------------------------------- GFF parsing --------------------------------

@dataclass
class CdsPart:
    """A single CDS feature line."""
    seqid: str
    start: int      # 1-based inclusive
    end: int        # 1-based inclusive
    strand: str     # '+' or '-'
    phase: int      # 0,1,2
    parent: str     # mRNA ID


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


def parse_gff(gff_lines: Iterable[str]) -> Tuple[Dict[str, Mrna], Dict[str, int]]:
    """Parse mRNA and CDS features from a GFF3 stream.

    Parameters
    ----------
    gff_lines : Iterable[str]
        Lines of a GFF3 file.

    Returns
    -------
    tuple
        (mRNAs, contig_annotation_counts)
        mRNAs: map mRNA ID -> Mrna object with ordered CDS parts list (order not yet sorted).
        contig_annotation_counts: map seqid -> number of distinct mRNA annotations.
    """
    mrnas: Dict[str, Mrna] = {}
    contig_counts: Dict[str, int] = {}

    for ln in gff_lines:
        if not ln or ln.startswith("#"):
            continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        seqid, source, ftype, start_s, end_s, score, strand, phase_s, attrs_s = cols[:9]
        start = int(start_s)
        end = int(end_s)
        attrs: Dict[str, str] = {}
        for kv in attrs_s.split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                attrs[k] = v

        if ftype == "mRNA":
            mrna_id = attrs.get("ID")
            if not mrna_id:
                logging.warning("mRNA feature missing ID on %s:%s:%s-%s", seqid, strand, start, end)
                continue
            mrnas[mrna_id] = Mrna(
                mrna_id=mrna_id,
                seqid=seqid,
                strand=strand,
                start=start,
                end=end,
                attrs=attrs,
                cds_parts=[],
            )
            contig_counts[seqid] = contig_counts.get(seqid, 0) + 1

        elif ftype == "CDS":
            parent = attrs.get("Parent")
            if not parent:
                continue
            try:
                phase = int(phase_s) if phase_s in {"0", "1", "2"} else 0
            except Exception:
                phase = 0
            part = CdsPart(seqid=seqid, start=start, end=end, strand=strand, phase=phase, parent=parent)
            # It's possible the mRNA appears later; stash anyway
            if parent not in mrnas:
                mrnas[parent] = Mrna(
                    mrna_id=parent,
                    seqid=seqid,
                    strand=strand,
                    start=start,
                    end=end,
                    attrs={},
                    cds_parts=[part],
                )
                contig_counts[seqid] = contig_counts.get(seqid, 0) + 1
            else:
                mrnas[parent].cds_parts.append(part)

    return mrnas, contig_counts


# ------------------------------ Sequence utils ------------------------------

def reverse_complement(seq: str) -> str:
    """Return reverse complement for A/C/G/T/N sequences (others -> N)."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    # Any other symbol becomes 'N' to maintain length:
    seq_rc = seq.translate(comp)
    # Replace any untranslated chars (non-ACGTN) with 'N'
    return "".join(ch if ch in "ACGTNacgtn" else "N" for ch in seq_rc[::-1])


def fetch_cds_sequence(fa: Fasta, cds_parts: List[CdsPart], strand: str) -> str:
    """Fetch and concatenate CDS nucleotide sequence in transcript orientation.

    GFF3 uses 1-based inclusive coordinates. pyfaidx uses 0-based, end-exclusive
    slicing. We therefore slice [start-1 : end] for each part.

    Phase handling (GFF3 spec):
    - For '+' strand: remove `phase` bases from the **left** of each CDS part.
    - For '-' strand: reverse-complement each part first (to transcript 5'→3'),
      then remove `phase` bases from the **left** of that RC part.

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
    if strand not in {"+", "-"}:
        raise ValueError(f"Invalid strand: {strand}")

    # Sort to transcript order
    parts = sorted(cds_parts, key=lambda p: p.start, reverse=(strand == "-"))

    chunks: List[str] = []
    for part in parts:
        # Pyfaidx slice: start-1 .. end (end-exclusive)
        # pyfaidx returns upper/lower as in file; force upper for consistency.
        s = str(fa[part.seqid][part.start - 1 : part.end]).upper()

        if strand == "-":
            s = reverse_complement(s)

        # Apply phase as "trim from left of the transcript-oriented piece"
        if part.phase in (1, 2):
            if len(s) >= part.phase:
                s = s[part.phase :]
            else:
                logging.warning(
                    "CDS part shorter than phase (seqid=%s start=%d end=%d phase=%d)",
                    part.seqid, part.start, part.end, part.phase
                )
                s = ""

        chunks.append(s)

    return "".join(chunks)


def count_non_acgt(seq: str) -> int:
    """Count non-ACGT letters (after uppercasing and U->T)."""
    seq = seq.upper().replace("U", "T")
    return sum(1 for ch in seq if ch not in {"A", "C", "G", "T"})


# --------------------------------- FASTA I/O --------------------------------

def write_fasta(records: Iterable[Tuple[str, str]], handle: TextIO, width: int = 60) -> None:
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
        handle.write(f">{header}\n")
        for i in range(0, len(seq), width):
            handle.write(seq[i : i + width] + "\n")


# ---------------------------------- CLI/Main --------------------------------

def configure_logging(level: str = "INFO") -> None:
    """Configure root logger to stderr with a concise format."""
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        stream=sys.stderr,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    """Build CLI parser."""
    ap = argparse.ArgumentParser(
        description=f"Extract proteins from maxiprot GFF3 (v{__version__})"
    )
    ap.add_argument("gff", nargs="?", default="-", help="GFF3 path or '-' for stdin (default: '-')")
    ap.add_argument("-g", "--genome", required=True, help="Path to genome FASTA (plain or bgzip-compressed)")
    ap.add_argument("--out-faa", default="-", help="Output protein FASTA path (default: '-' = stdout)")
    ap.add_argument("--transl-table", type=int, default=1, help="NCBI translation table ID (default: 1)")
    ap.add_argument(
        "--max-annos-per-contig", type=int, default=0,
        help="Max number of mRNA annotations allowed per contig; "
             "if a contig exceeds this number, all its proteins are skipped. "
             "0 or negative = unlimited (default: 0)."
    )
    ap.add_argument(
        "--exclude-pseudogenes", action="store_true",
        help="Exclude translations that contain premature '*' (internal stop codon)."
    )
    ap.add_argument("--log-level", default="INFO", help="Logging level (DEBUG, INFO, WARNING, ERROR)")
    return ap


def open_gff_lines(path: str) -> Iterable[str]:
    """Yield lines from a GFF path or stdin."""
    if path == "-" or path is None:
        for ln in sys.stdin:
            yield ln
    else:
        with open(path, "r", encoding="utf-8") as fh:
            for ln in fh:
                yield ln


def _open_out_handle(path: str) -> Tuple[TextIO, bool]:
    """Open output handle (stdout if '-' or empty)."""
    if path in (None, "", "-"):
        return sys.stdout, False
    return open(path, "w", encoding="utf-8"), True


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entrypoint."""
    ap = build_arg_parser()
    args = ap.parse_args(argv)

    configure_logging(args.log_level)
    logging.info("maxiprot_extract_proteins version %s", __version__)
    logging.info("Reading GFF: %s", args.gff)
    logging.info("Genome FASTA: %s", args.genome)

    # Parse GFF and gather mRNAs/CDS
    mrnas, contig_counts = parse_gff(open_gff_lines(args.gff))
    logging.info("Parsed %d mRNA annotations across %d contigs", len(mrnas), len(contig_counts))

    # Contig filter based on max annotations per contig
    disallowed_contigs: set[str] = set()
    if args.max_annos_per_contig and args.max_annos_per_contig > 0:
        for seqid, n in contig_counts.items():
            if n > args.max_annos_per_contig:
                disallowed_contigs.add(seqid)
        if disallowed_contigs:
            logging.warning(
                "Skipping all proteins on %d contig(s) due to --max-annos-per-contig=%d: %s",
                len(disallowed_contigs), args.max_annos_per_contig, ", ".join(sorted(disallowed_contigs))
            )

    # Open genome with pyfaidx (auto-creates .fai if missing)
    fai_path = Path(args.genome + ".fai")
    if not fai_path.exists():
        logging.info("No FASTA index found (%s); creating with pyfaidx...", fai_path.name)
    fa = Fasta(args.genome, as_raw=True, sequence_always_upper=True)  # pyfaidx will index if needed

    # Prepare output
    out_handle, close_out = _open_out_handle(args.out_faa)

    # Translate per mRNA
    n_out = 0
    n_skipped_contig = 0
    n_skipped_pseudo = 0

    try:
        records: List[Tuple[str, str]] = []
        for mrna_id, m in mrnas.items():
            if m.seqid in disallowed_contigs:
                n_skipped_contig += 1
                continue
            if not m.cds_parts:
                logging.warning("mRNA %s has no CDS parts; skipped", mrna_id)
                continue

            # Fetch CDS nt sequence (handles phase and strand)
            cds_nt = fetch_cds_sequence(fa, m.cds_parts, m.strand)

            # Sanity checks
            non_acgt = count_non_acgt(cds_nt)
            if non_acgt > 0:
                logging.warning("CDS for %s contains %d non-ACGT bases", mrna_id, non_acgt)

            if len(cds_nt) % 3 != 0:
                logging.warning("CDS length not divisible by 3 for %s (len=%d) — truncating remainder",
                                mrna_id, len(cds_nt))
                cds_nt = cds_nt[: len(cds_nt) - (len(cds_nt) % 3)]

            # Translate
            try:
                prot = translate_nt(cds_nt, args.transl_table)
            except ValueError as e:
                logging.error(str(e))
                return 2

            # Pseudogene logic
            internal_stop = "*" in prot[:-1] if prot else False
            if internal_stop and args.exclude_pseudogenes:
                n_skipped_pseudo += 1
                continue

            # FASTA header
            target = m.attrs.get("Target", "")
            header_parts = [
                mrna_id,
                f"{m.seqid}:{m.start}-{m.end}({m.strand})",
                f"table={args.transl_table}",
            ]
            if target:
                header_parts.append(f"target={target}")
            if internal_stop:
                header_parts.append("pseudogene=1")
            header = " ".join(header_parts)

            records.append((header, prot))
            n_out += 1

        # Write output FASTA (streaming keeps memory small; here we batch in one go)
        write_fasta(records, out_handle)
    finally:
        if close_out:
            out_handle.close()
        fa.close()

    logging.info("Wrote %d protein(s)", n_out)
    if n_skipped_contig:
        logging.info("Skipped %d protein(s) due to contig annotation limit", n_skipped_contig)
    if n_skipped_pseudo:
        logging.info("Excluded %d pseudogene translation(s) with internal stops", n_skipped_pseudo)

    return 0


if __name__ == "__main__":
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
