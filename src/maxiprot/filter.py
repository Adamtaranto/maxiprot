#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pick the best miniprot alignment per locus from a miniprot GFF3.

Default outputs:
- GFF3 → **stdout** (unless --out-gff3 PATH is provided)
- TSV  → **stderr** (unless --out-tsv PATH is provided)

This tool reads a miniprot GFF3 that includes ``##PAF`` header lines, computes
alignment metrics, applies user-configurable scoring and gating, clusters
alignments into loci on the same sequence/strand, and selects one "best" or
"longest" candidate per locus.

Outputs
-------
1) A GFF3 (to stdout by default) with a valid gene → mRNA → CDS hierarchy for
   the chosen candidate in each locus. Multi-line CDS features **share the same ID**,
   as required by the GFF3 specification.
2) A TSV summary (to stderr by default) with one row per chosen locus.

Notes
-----
- Requires miniprot GFF with ``##PAF`` header lines (contains tags like
  ``AS:i``, ``ms:i``, ``np:i``, ``fs:i``, ``st:i``, ``cg:Z``, ``cs:Z``).
- Coverage gate is applied to *reference (query) protein* coverage.
- Identity is computed from ``cs:Z:`` as identical-AA count divided by aligned AA.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import logging
import re
import sys
from typing import Dict, List, Optional, Sequence, TextIO, Tuple

import pandas as pd

from maxiprot._version import __version__

# ---------------------------------------------------------------------------
# Regexes reused across functions (compiled once)
# ---------------------------------------------------------------------------

_RE_CG = re.compile(r'(\d+)([MIDNUVFG])')
_RE_CS_MATCHES = re.compile(r':(\d+)')  # e.g., ':128' blocks -> identical AA counts
_RE_CS_SUBS = re.compile(
    r'\*'
)  # '*' ops -> mismatches/substitutions (not used in score)


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ParsedCg:
    """
    Parsed components of miniprot PAF ``cg:Z`` string.

    Parameters
    ----------
    M, I, D, G : int
        Amino acid operation lengths for Match, Insertion, Deletion, and 'G'.
        We treat M/I/D/G as contributing to aligned AA length.
    intron_nt : int
        Total intronic nucleotides (N/U/V ops) for reporting; not used in scoring.
    """

    M: int = 0
    I: int = 0
    D: int = 0
    G: int = 0
    intron_nt: int = 0


# ---------------------------------------------------------------------------
# Utilities & parsing
# ---------------------------------------------------------------------------


def parse_cg(cg: Optional[str]) -> ParsedCg:
    """
    Parse a miniprot ``cg:Z`` string into operation counts."""
    if not isinstance(cg, str):
        return ParsedCg()
    M = I = D = G = intron_nt = 0
    for n_str, op in _RE_CG.findall(cg):
        n = int(n_str)
        if op == 'M':
            M += n
        elif op == 'I':
            I += n
        elif op == 'D':
            D += n
        elif op == 'G':
            G += n
        elif op in ('N', 'U', 'V'):
            intron_nt += n
    return ParsedCg(M=M, I=I, D=D, G=G, intron_nt=intron_nt)


def parse_cs(cs: Optional[str]) -> Tuple[int, int]:
    """
    Parse a miniprot ``cs:Z`` string to count identity and substitutions.

    Returns
    -------
    (int, int)
        (identical_aa_count, substitution_count)
    """
    if not isinstance(cs, str):
        return 0, 0
    aa_ident = sum(map(int, _RE_CS_MATCHES.findall(cs)))
    aa_subs = len(_RE_CS_SUBS.findall(cs))
    return aa_ident, aa_subs


def read_lines(source: str) -> List[str]:
    """
    Read all lines from a file path or stdin."""
    if source == '-' or source is None:
        return sys.stdin.read().splitlines()
    with open(source, 'r', encoding='utf-8') as fh:
        return fh.read().splitlines()


def read_paf_from_gff_lines(lines: Sequence[str]) -> pd.DataFrame:
    """
    Extract PAF-like fields from GFF3 ``##PAF\t...`` header lines."""
    rows: List[Dict[str, object]] = []

    def to_int(k: str) -> int:
        try:
            return int(tagd.get(k, 0))
        except Exception:
            return 0

    for ln in lines:
        if not ln:
            continue
        if ln.startswith('##PAF'):
            parts = ln.split('\t')[1:]
        else:
            continue
        if len(parts) < 12:
            continue
        core = parts[:12]
        tags = parts[12:]
        qname = core[0]
        qlen = int(core[1])
        qstart = int(core[2])
        qend = int(core[3])
        strand = core[4]
        tname = core[5]
        tlen = int(core[6])
        tstart = int(core[7])
        tend = int(core[8])
        nmatch_nt = int(core[9])
        aln_nt_no_introns = int(core[10])
        mapq = int(core[11])
        tagd: Dict[str, str] = {}
        for kv in tags:
            if kv.count(':') >= 2:
                k, _typ, val = kv.split(':', 2)
                tagd[k] = val

        rows.append(
            {
                'qname': qname,
                'qlen': qlen,
                'qstart': qstart,
                'qend': qend,
                'strand': strand,
                'tname': tname,
                'tlen': tlen,
                'tstart': tstart,
                'tend': tend,
                'nmatch_nt': nmatch_nt,
                'aln_nt_no_introns': aln_nt_no_introns,
                'mapq': mapq,
                'AS': to_int('AS'),
                'ms': to_int('ms'),
                'np': to_int('np'),
                'fs': to_int('fs'),
                'st': to_int('st'),
                'cg': tagd.get('cg'),
                'cs': tagd.get('cs'),
            }
        )
    if not rows:
        raise SystemExit(
            'ERROR: No PAF header lines (##PAF) found. '
            'Run miniprot with GFF output that includes PAF headers.'
        )
    return pd.DataFrame(rows)


def read_gff_features(lines: Sequence[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse mRNA and CDS features from miniprot GFF3 body lines.
    """
    mrows: List[Dict[str, object]] = []
    crows: List[Dict[str, object]] = []
    for ln in lines:
        if not ln or ln.startswith('#'):
            continue
        parts = ln.split('\t')
        if len(parts) < 9:
            continue
        seqid, source, ftype, start, end, score, strand, phase, attrs = parts[:9]
        start_i, end_i = int(start), int(end)

        # Attribute parsing (tolerant)
        ad: Dict[str, str] = {}
        for kv in attrs.strip().split(';'):
            if '=' in kv:
                k, v = kv.split('=', 1)
                ad[k] = v

        if ftype == 'mRNA':
            tgt = ad.get('Target', '')
            tgt_qname = tgt.split()[0] if tgt else None
            mrows.append(
                {
                    'seqid': seqid,
                    'source': source,
                    'start': start_i,
                    'end': end_i,
                    'strand': strand,
                    'score': score,
                    'phase': phase,
                    'attrs': ad,
                    'mrna_id': ad.get('ID'),
                    'target_qname': tgt_qname,
                }
            )
        elif ftype == 'CDS':
            crows.append(
                {
                    'seqid': seqid,
                    'source': source,
                    'start': start_i,
                    'end': end_i,
                    'strand': strand,
                    'score': score,
                    'phase': phase,
                    'attrs': ad,
                    'parent': ad.get('Parent'),
                    'raw_attr': attrs,
                }
            )

    mdf = pd.DataFrame(mrows)
    cdf = pd.DataFrame(crows)
    if not cdf.empty:
        cdf['len_nt'] = cdf['end'] - cdf['start'] + 1
    return mdf, cdf


# ---------------------------------------------------------------------------
# Metrics & scoring
# ---------------------------------------------------------------------------


def add_alignment_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Append derived alignment metrics to the PAF dataframe.
    """
    cg_parsed = df['cg'].apply(parse_cg)
    df[['M', 'I', 'D', 'G', 'intron_nt']] = pd.DataFrame(
        [(p.M, p.I, p.D, p.G, p.intron_nt) for p in cg_parsed], index=df.index
    )
    df['aa_aligned'] = (df['M'] + df['I'] + df['D'] + df['G']).clip(lower=1)
    ident_subs = df['cs'].apply(parse_cs)
    df[['aa_ident', 'aa_subs']] = pd.DataFrame(list(ident_subs), index=df.index)
    df['cov_aa'] = (df['qend'] - df['qstart']) / df['qlen'].clip(lower=1)
    df['pid_aa'] = (df['aa_ident'] / df['aa_aligned']).clip(0, 1)
    df['positives'] = (df['np'] / df['aa_aligned']).clip(0, 1)
    df['ms_per_qlen'] = df['ms'] / df['qlen'].clip(lower=1)
    df['AS_per_qlen'] = df['AS'] / df['qlen'].clip(lower=1)
    df['len_frac'] = df['aa_aligned'] / df['qlen'].clip(lower=1)
    df['cds_aa_len'] = df['aa_aligned']
    return df


def _linear(row: pd.Series, args: argparse.Namespace, length_metric: str) -> float:
    return (
        args.w_pid * float(row['pid_aa'])
        + args.w_cov * float(row['cov_aa'])
        + args.w_pos * float(row['positives'])
        + args.w_ms * float(row['ms_per_qlen'])
        + args.w_AS * float(row['AS_per_qlen'])
        + args.w_len * float(row[length_metric])
    )


def _geom(
    row: pd.Series, args: argparse.Namespace, length_metric: str, eps: float = 1e-9
) -> float:
    factors: List[float] = []
    for w, val in (
        (args.w_pid, float(row['pid_aa'])),
        (args.w_cov, float(row['cov_aa'])),
        (args.w_pos, float(row['positives'])),
        (args.w_ms, float(row['ms_per_qlen'])),
        (args.w_AS, float(row['AS_per_qlen'])),
        (args.w_len, float(row[length_metric])),
    ):
        if w != 0:
            factors.append(max(val, eps) ** w)
    if not factors:
        return 0.0
    prod = 1.0
    for f in factors:
        prod *= f
    return float(prod)


def score_alignment(row: pd.Series, mode: str, args: argparse.Namespace) -> float:
    """
    Compute a score for a candidate row.
    """
    length_metric = 'len_frac' if args.length_metric == 'frac' else 'cds_aa_len'

    if mode == 'ms_cov_pos':
        return (
            float(row['ms_per_qlen'])
            * float(row['cov_aa'])
            * (0.5 + 0.5 * float(row['positives']))
        )
    if mode == 'AS':
        return float(row['AS_per_qlen'])
    if mode == 'ms':
        return float(row['ms_per_qlen'])
    if mode == 'pid_cov':
        return float(row['pid_aa']) * float(row['cov_aa'])
    if mode == 'pid_cov_len':
        return (
            (float(row['pid_aa']) ** max(args.w_pid, 0.0))
            * (float(row['cov_aa']) ** max(args.w_cov, 0.0))
            * (float(row[length_metric]) ** max(args.w_len, 0.0))
        )
    if mode == 'length':
        return float(row[length_metric])
    if mode == 'linear':
        return _linear(row, args, length_metric)
    if mode == 'geom':
        return _geom(row, args, length_metric)

    return (
        float(row['ms_per_qlen'])
        * float(row['cov_aa'])
        * (0.5 + 0.5 * float(row['positives']))
    )


def apply_gates(
    df: pd.DataFrame, cov_min: float = 0.60, pid_min: float = 0.30
) -> pd.DataFrame:
    """Apply gating thresholds for coverage and identity."""
    out = df.copy()
    out['passes'] = (out['cov_aa'] >= cov_min) & (out['pid_aa'] >= pid_min)
    return out


def cluster_into_loci(df: pd.DataFrame, pad_nt: int = 5000) -> pd.DataFrame:
    """Cluster hits into loci along each target sequence and strand."""
    out = df.copy()
    out['locus'] = None
    for (tname, strand), sub in df.sort_values(
        ['tname', 'strand', 'tstart', 'tend']
    ).groupby(['tname', 'strand'], sort=False):
        current_end = -(10**18)
        cid = 0
        for i, r in sub.iterrows():
            if r.tstart <= current_end + pad_nt:
                current_end = max(current_end, r.tend)
            else:
                cid += 1
                current_end = r.tend
            out.at[i, 'locus'] = f'{tname}:{strand}:{cid}'
    return out


# ---------------------------------------------------------------------------
# Mapping winners back to GFF features
# ---------------------------------------------------------------------------


def jaccard(a0: int, a1: int, b0: int, b1: int) -> float:
    """
    Compute Jaccard index for two closed intervals [a0,a1], [b0,b1].
    """
    inter = max(0, min(a1, b1) - max(a0, b0) + 1)
    uni = (a1 - a0 + 1) + (b1 - b0 + 1) - inter
    return inter / uni if uni > 0 else 0.0


def attach_mrna_and_cds_length(
    winners: pd.DataFrame, mdf: pd.DataFrame, cdf: pd.DataFrame
) -> pd.DataFrame:
    """
    Find the matching mRNA in the GFF and compute total CDS span for winners.
    """
    mrna_ids: List[Optional[str]] = []
    cds_sums: List[Optional[int]] = []

    for _, r in winners.iterrows():
        cand = mdf[
            (mdf['seqid'] == r['tname'])
            & (mdf['strand'] == r['strand'])
            & (mdf['target_qname'] == r['qname'])
        ]
        if cand.empty:
            mrna_ids.append(None)
            cds_sums.append(None)
            continue
        cand = cand.copy()
        cand['ovl'] = cand.apply(
            lambda x, tstart=int(r['tstart']), tend=int(r['tend']): jaccard(
                tstart, tend, int(x['start']), int(x['end'])
            ),
            axis=1,
        )
        cand = cand.sort_values('ovl', ascending=False)
        mrna_id = (
            str(cand.iloc[0]['mrna_id'])
            if not cand.empty and cand.iloc[0]['ovl'] > 0
            else None
        )
        if mrna_id and not cdf.empty:
            cds_sum = int(cdf.loc[cdf['parent'] == mrna_id, 'len_nt'].sum())
        else:
            cds_sum = None
        mrna_ids.append(mrna_id)
        cds_sums.append(cds_sum)

    out = winners.copy()
    out['mrna_id'] = mrna_ids
    out['gff_cds_nt_len'] = cds_sums
    return out


# ---------------------------------------------------------------------------
# GFF3 writing
# ---------------------------------------------------------------------------


def gff3_escape(s: str) -> str:
    """
    Escape attribute values for GFF3 (percent-first order).
    """
    return (
        s.replace('%', '%25')
        .replace(';', '%3B')
        .replace('=', '%3D')
        .replace(',', '%2C')
        .replace('&', '%26')
    )


def _open_out_handle(path: Optional[str], default: TextIO) -> Tuple[TextIO, bool]:
    """
    Return a writable handle and whether we own/should close it.
    """
    if path is None or path == '-' or path == '':
        return default, False
    fh = open(path, 'w', encoding='utf-8')
    return fh, True


def write_best_gff3(
    out_path: Optional[str],
    winners: pd.DataFrame,
    mdf: pd.DataFrame,
    cdf: pd.DataFrame,
    id_prefix: str = 'PBM',
) -> None:
    """
    Write the best-per-locus annotations as a valid GFF3 hierarchy.

    If ``out_path`` is None/'-' -> write to stdout.
    """
    out, close_me = _open_out_handle(out_path, sys.stdout)
    try:
        out.write('##gff-version 3\n')
        for _, r in winners.iterrows():
            if pd.isna(r.get('mrna_id')):
                # Synthesize a minimal hierarchy from PAF coords.
                gene_id = f'{id_prefix}:gene:{gff3_escape(str(r["locus"]))}'
                mrna_id = f'{id_prefix}:mrna:{gff3_escape(str(r["locus"]))}'
                seqid = str(r['tname'])
                strand = str(r['strand'])
                start = int(min(r['tstart'], r['tend']))
                end = int(max(r['tstart'], r['tend']))
                out.write(
                    f'{seqid}\tminiprot\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_id};Name={gene_id}\n'
                )
                attrs = f'ID={mrna_id};Parent={gene_id};Target={gff3_escape(str(r["qname"]))} {int(r["qstart"]) + 1} {int(r["qend"])}'
                out.write(
                    f'{seqid}\tminiprot\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n'
                )
                cds_id = f'{id_prefix}:cds:{gff3_escape(str(r["locus"]))}'
                out.write(
                    f'{seqid}\tminiprot\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID={cds_id};Parent={mrna_id}\n'
                )
                continue

            mrna_id = str(r['mrna_id'])
            m = mdf.loc[mdf['mrna_id'] == mrna_id].iloc[0]
            seqid = str(m['seqid'])
            strand = str(m['strand'])
            start = int(min(m['start'], m['end']))
            end = int(max(m['start'], m['end']))

            gene_id = f'{id_prefix}:gene:{gff3_escape(str(r["locus"]))}'
            out.write(
                f'{seqid}\tminiprot\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_id};Name={gene_id}\n'
            )

            # mRNA: preserve Target and Identity if present
            attrs = [f'ID={mrna_id}', f'Parent={gene_id}']
            tgt = m['attrs'].get('Target')
            if tgt:
                attrs.append(f'Target={gff3_escape(str(tgt))}')
            ident = m['attrs'].get('Identity')
            if ident:
                attrs.append(f'Identity={ident}')
            out.write(
                f'{seqid}\tminiprot\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t{";".join(attrs)}\n'
            )

            # Multi-line CDS must share the same ID
            cds_id = f'{id_prefix}:cds:{gff3_escape(mrna_id)}'
            cds_parts = cdf.loc[cdf['parent'] == mrna_id].sort_values(['start', 'end'])
            for _, c in cds_parts.iterrows():
                cstart, cend = int(c['start']), int(c['end'])
                phase = str(c['phase']) if str(c['phase']) in ('0', '1', '2') else '0'
                out.write(
                    f'{seqid}\tminiprot\tCDS\t{cstart}\t{cend}\t.\t{strand}\t{phase}\tID={cds_id};Parent={mrna_id}\n'
                )
    finally:
        if close_me:
            out.close()


# ---------------------------------------------------------------------------
# CLI / main
# ---------------------------------------------------------------------------


def configure_logging(level: str = 'INFO') -> None:
    """
    Configure root logger with a standard format (logs to stderr).
    """
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format='%(asctime)s [%(levelname)s] %(message)s',
        stream=sys.stderr,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    """
    Build the command-line interface parser.
    """
    ap = argparse.ArgumentParser(
        description=(
            'Pick best miniprot alignments per locus from a GFF3 (with ##PAF headers). '
            f'Version {__version__}'
        )
    )
    ap.add_argument(
        'gff',
        nargs='?',
        default='-',
        help="miniprot GFF3 path or '-' for stdin (default: '-')",
    )

    # Scoring
    ap.add_argument(
        '--score-mode',
        choices=[
            'ms_cov_pos',
            'AS',
            'ms',
            'pid_cov',
            'pid_cov_len',
            'length',
            'linear',
            'geom',
        ],
        default='ms_cov_pos',
        help='Scoring function to rank candidates (default: ms_cov_pos)',
    )
    ap.add_argument(
        '--length-metric',
        choices=['frac', 'aa'],
        default='frac',
        help="Use 'frac' (aa_len/qlen) or 'aa' (absolute aa length) when a score uses length (default: frac)",
    )
    ap.add_argument(
        '--w-pid',
        type=float,
        default=1.0,
        help='Weight/exponent for identity where relevant (default: 1.0)',
    )
    ap.add_argument(
        '--w-cov',
        type=float,
        default=1.0,
        help='Weight/exponent for coverage where relevant (default: 1.0)',
    )
    ap.add_argument(
        '--w-len',
        type=float,
        default=1.0,
        help='Weight/exponent for length where relevant (default: 1.0)',
    )
    ap.add_argument(
        '--w-pos',
        type=float,
        default=0.0,
        help='Weight for positives in linear/geom modes (default: 0.0)',
    )
    ap.add_argument(
        '--w-ms',
        type=float,
        default=0.0,
        help='Weight for ms_per_qlen in linear/geom modes (default: 0.0)',
    )
    ap.add_argument(
        '--w-AS',
        type=float,
        default=0.0,
        help='Weight for AS_per_qlen in linear/geom modes (default: 0.0)',
    )

    # Gates
    ap.add_argument(
        '--min-cov',
        type=float,
        default=0.60,
        help='Minimum AA coverage of query to pass (default: 0.60)',
    )
    ap.add_argument(
        '--min-pid',
        type=float,
        default=0.80,
        help='Minimum AA identity to pass (default: 0.30)',
    )
    ap.add_argument(
        '--strict',
        action='store_true',
        help='Drop loci with no passing candidates (default: pick best anyway)',
    )

    # Locus clustering
    ap.add_argument(
        '--locus-pad',
        dest='locus_pad',
        type=int,
        default=5000,
        help='Max gap (nt) to cluster non-overlapping hits to the same locus on the same strand (default: 5000)',
    )
    ap.add_argument(
        '--locus-gap', dest='locus_pad', type=int, help='Alias for --locus-pad'
    )

    # Selection policy
    ap.add_argument(
        '--selection-mode',
        choices=['best', 'prefer_intact', 'longest', 'longest_prefer_intact'],
        default='best',
        help=(
            'How to choose the winner within a locus:\n'
            '  best: highest score (gated first)\n'
            '  prefer_intact: best among intact if available; else best overall\n'
            '  longest: longest AA length; tie-break by score\n'
            '  longest_prefer_intact: longest among intact; else longest overall'
        ),
    )

    # Output / misc
    ap.add_argument(
        '--id-prefix',
        default='PBM',
        help='Prefix for synthesized Gene/CDS IDs (default: PBM)',
    )
    ap.add_argument(
        '--out-gff3',
        default='-',
        help="Output GFF3 path (default: '-' = stdout)",
    )
    ap.add_argument(
        '--out-tsv',
        default=None,
        help="Output TSV path (default: write TSV to stderr). Use '-' to also force stderr.",
    )
    ap.add_argument(
        '--log-level',
        default='INFO',
        help='Logging level (DEBUG, INFO, WARNING, ERROR; default: INFO)',
    )
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    CLI entrypoint.

    Returns
    -------
    int
        Exit status code.
    """
    ap = build_arg_parser()
    args = ap.parse_args(argv)

    configure_logging(args.log_level)
    logging.info('maxiprot version %s', __version__)
    logging.info('Reading input: %s', args.gff)

    lines = read_lines(args.gff)

    logging.info('Parsing PAF headers (##PAF lines)')
    paf = read_paf_from_gff_lines(lines)
    logging.info('Found %d candidate alignments in PAF', len(paf))

    logging.info('Parsing GFF features (mRNA/CDS)')
    mdf, cdf = read_gff_features(lines)
    logging.info('GFF mRNAs: %d; CDS parts: %d', len(mdf), len(cdf))

    logging.info('Computing alignment metrics and gates')
    paf = add_alignment_metrics(paf)
    paf['score_raw'] = paf.apply(
        lambda r: score_alignment(r, args.score_mode, args), axis=1
    )
    paf['passes'] = (paf['cov_aa'] >= args.min_cov) & (paf['pid_aa'] >= args.min_pid)

    logging.info(
        'Clustering into loci with pad=%d nt (same seq & strand)', args.locus_pad
    )
    paf = cluster_into_loci(paf, pad_nt=args.locus_pad)

    # Selection per locus
    winners_rows: List[pd.Series] = []
    logging.info(
        'Selection mode: %s | score-mode: %s | gates: cov>=%.2f pid>=%.2f',
        args.selection_mode,
        args.score_mode,
        args.min_cov,
        args.min_pid,
    )
    if args.score_mode in ('pid_cov_len', 'linear', 'geom'):
        logging.info(
            'Weights: w_pid=%.3f w_cov=%.3f w_len=%.3f w_pos=%.3f w_ms=%.3f w_AS=%.3f | length-metric=%s',
            args.w_pid,
            args.w_cov,
            args.w_len,
            args.w_pos,
            args.w_ms,
            args.w_AS,
            args.length_metric,
        )

    for locus, sub in paf.groupby('locus', sort=False):
        num = len(sub)
        logging.info('[Locus %s] %d candidates', locus, num)

        # Sort: passers first, then by score, then by length (as tie-break)
        sub = sub.copy().sort_values(
            ['passes', 'score_raw', 'cds_aa_len'], ascending=[False, False, False]
        )

        # Per-candidate log
        for _, r in sub.iterrows():
            status = (
                'pseudogene' if (int(r['fs']) > 0 or int(r['st']) > 0) else 'intact'
            )
            logging.info(
                '[Locus %s] cand q=%s score=%.5f cov=%.3f pid=%.3f lenAA=%d pos=%.3f ms=%d AS=%d fs=%d st=%d pass=%s status=%s',
                locus,
                r['qname'],
                float(r['score_raw']),
                float(r['cov_aa']),
                float(r['pid_aa']),
                int(r['cds_aa_len']),
                float(r['positives']),
                int(r['ms']),
                int(r['AS']),
                int(r['fs']),
                int(r['st']),
                bool(r['passes']),
                status,
            )

        # Optional strict drop of non-pass loci
        if args.strict and not bool(sub['passes'].any()):
            logging.warning(
                '[Locus %s] no candidates pass gates -> skipped due to --strict', locus
            )
            continue

        # Selection policy
        def choose(df: pd.DataFrame) -> pd.Series:
            if args.selection_mode in ('longest', 'longest_prefer_intact'):
                return df.sort_values(
                    ['passes', 'cds_aa_len', 'score_raw'],
                    ascending=[False, False, False],
                ).iloc[0]
            return df.iloc[0]

        working = sub
        if args.selection_mode in ('prefer_intact', 'longest_prefer_intact'):
            intact = sub[(sub['fs'] == 0) & (sub['st'] == 0)]
            if not intact.empty:
                logging.info(
                    '[Locus %s] intact candidates available: %d -> prefer intact',
                    locus,
                    len(intact),
                )
                working = intact
            else:
                logging.info(
                    '[Locus %s] no intact candidates -> fallback to all', locus
                )

        pick = choose(working)
        logging.info(
            '[Locus %s] SELECTED q=%s (score=%.5f, cov=%.3f, pid=%.3f, lenAA=%d, fs=%d, st=%d)',
            locus,
            pick['qname'],
            float(pick['score_raw']),
            float(pick['cov_aa']),
            float(pick['pid_aa']),
            int(pick['cds_aa_len']),
            int(pick['fs']),
            int(pick['st']),
        )
        winners_rows.append(pick)

    if not winners_rows:
        logging.warning('No winners selected. Exiting.')
        return 0

    winners = pd.DataFrame(winners_rows)
    winners = winners.assign(
        status=winners.apply(
            lambda r: 'pseudogene'
            if (int(r['fs']) > 0 or int(r['st']) > 0)
            else 'intact',
            axis=1,
        )
    )

    # Map back to GFF features for chosen winners
    winners = attach_mrna_and_cds_length(winners, mdf, cdf)

    # Write TSV (file if given; else stderr)
    tsv_cols = [
        'locus',
        'tname',
        'tstart',
        'tend',
        'strand',
        'qname',
        'qlen',
        'cov_aa',
        'pid_aa',
        'positives',
        'ms',
        'AS',
        'score_raw',
        'cds_aa_len',
        'fs',
        'st',
        'status',
        'passes',
        'mapq',
        'mrna_id',
        'gff_cds_nt_len',
    ]

    if args.out_tsv and args.out_tsv != '-':
        winners[tsv_cols].to_csv(args.out_tsv, sep='\t', index=False)
        logging.info('Wrote TSV summary: %s', args.out_tsv)
    else:
        # Pure TSV to stderr (logs may interleave; consider --log-level ERROR or redirect)
        winners[tsv_cols].to_csv(sys.stderr, sep='\t', index=False)
        logging.info('Wrote TSV summary to stderr')

    # Write GFF3 (stdout by default)
    write_best_gff3(args.out_gff3, winners, mdf, cdf, id_prefix=args.id_prefix)
    if args.out_gff3 and args.out_gff3 != '-':
        logging.info('Wrote best-per-locus GFF3: %s', args.out_gff3)
    else:
        logging.info('Wrote best-per-locus GFF3 to stdout')

    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        # allow piping into head/tail without noisy tracebacks
        try:
            sys.stderr.close()
        except Exception:
            pass
        try:
            sys.stdout.close()
        except Exception:
            pass
        raise
