#!/usr/bin/env python3
"""
Pick the best miniprot alignment per locus from a miniprot GFF3.

Reads a miniprot GFF3 that includes ``##PAF`` header lines, computes
alignment metrics (identity from ``cs:Z``, query coverage, positives from
``np:i``), applies user-configurable scoring and gating, clusters alignments
into loci on the same sequence/strand, and selects candidates per locus.

Outputs
-------
1) A GFF3 (stdout by default) with a valid gene -> mRNA -> CDS hierarchy for
   each emitted candidate. Multi-line CDS features **share the same ID**, as
   required by the GFF3 specification.
2) An optional TSV summary (only with ``--out-tsv``; ``-`` = stdout) with one
   row per emitted candidate.

Notes
-----
- Requires miniprot GFF with ``##PAF`` header lines (tags like ``AS:i``,
  ``ms:i``, ``np:i``, ``fs:i``, ``st:i``, ``cs:Z``).
- Coverage is computed on the *reference protein* (query): ``(qe - qs)/qlen``.
"""

from __future__ import annotations

import argparse
import logging
from typing import Optional, Sequence

import pandas as pd

from maxiprot.emit import write_best_gff3
from maxiprot.gffio import GffFeature, PafRecord, iter_gff3
from maxiprot.ioutils import guard_broken_pipe, open_output, read_lines
from maxiprot.logs import LOG_LEVELS, init_logging
from maxiprot.scoring import (
    SCORE_MODES,
    Weights,
    compute_alignment_metrics,
    compute_scores,
)

logger = logging.getLogger(__name__)

PAF_CORE_COLS = (
    'qname',
    'qlen',
    'qs',
    'qe',
    'strand',
    'tname',
    'tlen',
    'ts',
    'te',
    'nmatch',
    'alen',
    'mapq',
)

TSV_COLUMNS = (
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
    'mrna_id',
    'gff_cds_nt_len',
    'passes',
    'mapq',
)


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------


def load_input(
    path: Optional[str],
) -> tuple[pd.DataFrame, dict[str, GffFeature], dict[str, list[GffFeature]]]:
    """
    Stream a miniprot GFF3 and collect PAF records plus mRNA/CDS features.

    Parameters
    ----------
    path : str or None
        Input path, or ``'-'``/None for stdin.

    Returns
    -------
    tuple
        ``(paf_df, mrna_by_id, cds_by_parent)`` where ``paf_df`` has the PAF
        core columns plus one column per tag seen in the input.
    """
    paf_rows: list[dict[str, object]] = []
    mrna_by_id: dict[str, GffFeature] = {}
    cds_by_parent: dict[str, list[GffFeature]] = {}

    for item in iter_gff3(read_lines(path)):
        if isinstance(item, PafRecord):
            row: dict[str, object] = {col: getattr(item, col) for col in PAF_CORE_COLS}
            row.update(item.tags)
            paf_rows.append(row)
        elif item.type == 'mRNA':
            fid = item.attrs.get('ID')
            if fid:
                if fid in mrna_by_id:
                    logger.warning('Duplicate mRNA ID in input: %s', fid)
                mrna_by_id[fid] = item
        elif item.type == 'CDS':
            parent = item.attrs.get('Parent')
            if parent:
                cds_by_parent.setdefault(parent, []).append(item)

    paf = pd.DataFrame.from_records(paf_rows)
    return paf, mrna_by_id, cds_by_parent


# ---------------------------------------------------------------------------
# Gating
# ---------------------------------------------------------------------------


def apply_gates(df: pd.DataFrame, min_cov: float, min_pid: float) -> pd.Series:
    """
    Apply pass/fail gates for coverage and percent identity.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with columns ``cov_aa`` and ``pid_aa``.
    min_cov : float
        Minimum fraction of reference coverage required (0-1).
    min_pid : float
        Minimum AA identity required (0-1).

    Returns
    -------
    pandas.Series
        Boolean series indicating pass/fail per row.
    """
    return (df['cov_aa'] >= min_cov) & (df['pid_aa'] >= min_pid)


def gate_intact_v_pseudo(df: pd.DataFrame) -> pd.Series:
    """
    Label candidates as intact (no frameshifts / in-frame stops) or not.

    Missing ``fs``/``st`` tags are treated as 0 (intact).

    Parameters
    ----------
    df : pandas.DataFrame
        PAF dataframe; ``fs`` and ``st`` columns are optional.

    Returns
    -------
    pandas.Series
        Boolean series, True for intact.
    """
    fs = _tag_series(df, 'fs')
    st = _tag_series(df, 'st')
    return (fs == 0) & (st == 0)


def _tag_series(df: pd.DataFrame, col: str) -> pd.Series:
    """Return an integer tag column, defaulting missing values to 0."""
    if col in df.columns:
        return pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
    return pd.Series(0, index=df.index, dtype=int)


# ---------------------------------------------------------------------------
# Locus clustering
# ---------------------------------------------------------------------------


def cluster_into_loci(paf: pd.DataFrame, locus_pad: int = 5000) -> pd.DataFrame:
    """
    Group alignments into loci on the same (tname, strand) within distance.

    Clustering is single-linkage chaining: a new locus starts when the next
    hit begins more than ``locus_pad`` nt after the running end of the current
    locus. Overlapping hits always share a locus. Tandem copies spaced closer
    than ``locus_pad`` therefore merge into one locus.

    Parameters
    ----------
    paf : pandas.DataFrame
        PAF-like dataframe (must have 'tname', 'strand', 'ts', 'te').
    locus_pad : int, default 5000
        Maximum gap (nt) joining consecutive hits into one locus.

    Returns
    -------
    pandas.DataFrame
        Input sorted by (tname, strand, ts, te) with an added ``locus`` label
        of the form ``tname:strand:index`` (index is 1-based per
        (tname, strand) group, in coordinate order).
    """
    df = paf.copy()
    if df.empty:
        df['locus'] = pd.Series(dtype=str)
        return df

    df = df.sort_values(['tname', 'strand', 'ts', 'te'], kind='mergesort')
    group_keys = [df['tname'], df['strand']]
    running_end = df.groupby(['tname', 'strand'], sort=False)['te'].cummax()
    prev_end = running_end.groupby(group_keys, sort=False).shift()
    new_locus = prev_end.isna() | (df['ts'] > prev_end + locus_pad)
    index_in_group = new_locus.groupby(group_keys, sort=False).cumsum().astype(int)
    df['locus'] = (
        df['tname'].astype(str)
        + ':'
        + df['strand'].astype(str)
        + ':'
        + index_in_group.astype(str)
    )
    return df


# ---------------------------------------------------------------------------
# Per-locus selection
# ---------------------------------------------------------------------------


def select_per_locus(
    df: pd.DataFrame,
    selection_mode: str = 'best',
    emit_mode: str = 'best',
    max_per_locus: Optional[int] = None,
    strict: bool = False,
) -> pd.DataFrame:
    """
    Select candidates to emit for each locus.

    Gates always exclude failing candidates from selection. A locus with no
    passing candidate is dropped under ``strict``; otherwise its single
    best-scoring candidate is emitted (flagged ``passes=False`` in the TSV).

    Parameters
    ----------
    df : pandas.DataFrame
        Scored, gated, clustered candidates (columns ``locus``, ``passes``,
        ``intact``, ``score_raw``, ``cds_aa_len``).
    selection_mode : str
        'best', 'prefer_intact', 'longest' or 'longest_prefer_intact'.
    emit_mode : str
        'best' (one winner per locus) or 'all_passing'.
    max_per_locus : int or None
        Cap on emitted candidates per locus for 'all_passing'.
    strict : bool
        Drop loci where no candidate passes gates.

    Returns
    -------
    pandas.DataFrame
        Winners with an added 1-based ``emit_rank`` column per locus.
    """
    frames: list[pd.DataFrame] = []
    for locus, sub in df.groupby('locus', sort=False):
        logger.debug('[Locus %s] %d candidate(s)', locus, len(sub))
        passers = sub[sub['passes']]
        if passers.empty:
            if strict:
                logger.info(
                    '[Locus %s] no candidates pass gates -> dropped (--strict)',
                    locus,
                )
                continue
            logger.info(
                '[Locus %s] no candidates pass gates -> emitting best-scoring '
                'candidate as fallback',
                locus,
            )
            pool = sub
            fallback = True
        else:
            pool = passers
            fallback = False

        if emit_mode == 'all_passing' and not fallback:
            emit = pool.sort_values(
                ['score_raw', 'cds_aa_len'], ascending=[False, False]
            )
            if max_per_locus is not None and max_per_locus > 0:
                emit = emit.head(max_per_locus)
        else:
            working = pool
            if selection_mode in ('prefer_intact', 'longest_prefer_intact'):
                intact = working[working['intact']]
                if not intact.empty:
                    working = intact
            if selection_mode in ('longest', 'longest_prefer_intact'):
                keys = ['cds_aa_len', 'score_raw']
            else:
                keys = ['score_raw', 'cds_aa_len']
            emit = working.sort_values(keys, ascending=[False, False]).head(1)

        for _, pick in emit.iterrows():
            logger.info(
                '[Locus %s] SELECTED q=%s (score=%.5f cov=%.3f pid=%.3f lenAA=%d passes=%s)',
                locus,
                pick['qname'],
                float(pick['score_raw']),
                float(pick['cov_aa']),
                float(pick['pid_aa']),
                int(pick['cds_aa_len']),
                'Y' if pick['passes'] else 'N',
            )
        frames.append(emit)

    if not frames:
        return df.head(0).assign(emit_rank=pd.Series(dtype=int))
    winners = pd.concat(frames)
    winners['emit_rank'] = winners.groupby('locus').cumcount() + 1
    return winners


# ---------------------------------------------------------------------------
# Linking winners back to input mRNA features
# ---------------------------------------------------------------------------


def _parse_target_name(s: Optional[str]) -> Optional[str]:
    """Return the first token of a GFF3 Target attribute."""
    if not s:
        return None
    s = str(s).strip()
    if not s:
        return None
    return s.split()[0]


def _interval_overlap(a1: int, a2: int, b1: int, b2: int) -> int:
    """1-based inclusive intervals; return overlap length (>=0)."""
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def attach_mrna_ids_to_winners(
    winners: pd.DataFrame,
    mrna_by_id: dict[str, GffFeature],
    cds_by_parent: dict[str, list[GffFeature]],
    overlap_thresh: float = 0.5,
) -> pd.DataFrame:
    """
    Attach ``mrna_id`` to winners by Target name and/or coordinate overlap.

    Strategy:
      1) Target name == PAF qname on the same seqid/strand, requiring a
         positive coordinate overlap (largest overlap wins).
      2) Otherwise the mRNA with the largest overlap of its CDS-union span
         with the PAF interval, accepted when
         ``overlap / min(len(paf), len(mrna_span)) >= overlap_thresh``.

    Parameters
    ----------
    winners : pandas.DataFrame
        Winner rows (``tname``, ``strand``, ``ts``, ``te``, ``qname``).
    mrna_by_id : dict of str to GffFeature
        Input mRNA features keyed by ID.
    cds_by_parent : dict of str to list of GffFeature
        Input CDS features grouped by Parent.
    overlap_thresh : float
        Minimum overlap fraction for the coordinate fallback.

    Returns
    -------
    pandas.DataFrame
        Winners with an ``mrna_id`` column (None where unmatched).
    """
    winners = winners.copy()

    # (seqid, strand) -> list of (mrna_id, span_start, span_end, target_name)
    by_key: dict[tuple[str, str], list[tuple[str, int, int, Optional[str]]]] = {}
    for fid, feat in mrna_by_id.items():
        parts = cds_by_parent.get(fid, [])
        if parts:
            span_start = min(p.start for p in parts)
            span_end = max(p.end for p in parts)
        else:
            span_start, span_end = feat.start, feat.end
        target_name = _parse_target_name(feat.attrs.get('Target'))
        by_key.setdefault((feat.seqid, feat.strand), []).append(
            (fid, span_start, span_end, target_name)
        )

    mrna_ids: list[Optional[str]] = []
    for _, w in winners.iterrows():
        candidates = by_key.get((str(w['tname']), str(w['strand'])), [])
        paf_start = int(w['ts']) + 1
        paf_end = int(w['te'])
        paf_len = max(0, paf_end - paf_start + 1)
        qname = str(w['qname'])

        chosen: Optional[str] = None
        best_ov = 0
        for fid, s, e, tn in candidates:
            if tn is not None and tn == qname:
                ov = _interval_overlap(paf_start, paf_end, s, e)
                if ov > best_ov:
                    best_ov = ov
                    chosen = fid

        if chosen is None:
            best = (0.0, 0)
            best_fid: Optional[str] = None
            for fid, s, e, _tn in candidates:
                ov = _interval_overlap(paf_start, paf_end, s, e)
                if ov <= 0:
                    continue
                mrna_len = e - s + 1
                frac = ov / max(1, min(paf_len, mrna_len))
                if (frac, ov) > best:
                    best = (frac, ov)
                    best_fid = fid
            if best_fid is not None and best[0] >= overlap_thresh:
                chosen = best_fid

        mrna_ids.append(chosen)

    winners['mrna_id'] = mrna_ids
    n_matched = sum(x is not None for x in mrna_ids)
    logger.info(
        'Linked %d/%d winners to existing mRNA features.', n_matched, len(winners)
    )
    return winners


# ---------------------------------------------------------------------------
# TSV output
# ---------------------------------------------------------------------------


def build_tsv(
    winners: pd.DataFrame,
    cds_by_parent: dict[str, list[GffFeature]],
) -> pd.DataFrame:
    """
    Build the TSV summary table (one row per emitted candidate).

    Parameters
    ----------
    winners : pandas.DataFrame
        Winner rows with metric columns and ``mrna_id``.
    cds_by_parent : dict of str to list of GffFeature
        Input CDS features grouped by Parent (for ``gff_cds_nt_len``).

    Returns
    -------
    pandas.DataFrame
        Table with :data:`TSV_COLUMNS`, floats rounded to 4 decimal places
        and missing values rendered as ``NA``.
    """
    df = winners.copy()
    df['tstart'] = df['ts'].astype(int) + 1
    df['tend'] = df['te'].astype(int)
    df['status'] = df['intact'].map({True: 'intact', False: 'pseudogene'})
    df['fs'] = _tag_series(df, 'fs')
    df['st'] = _tag_series(df, 'st')
    for col in ('ms', 'AS'):
        if col not in df.columns:
            df[col] = pd.NA
    df['gff_cds_nt_len'] = df['mrna_id'].map(
        lambda mid: (
            sum(p.end - p.start + 1 for p in cds_by_parent.get(mid, []))
            if pd.notna(mid) and mid in cds_by_parent
            else pd.NA
        )
    )
    for col in ('cov_aa', 'pid_aa', 'positives', 'score_raw'):
        df[col] = pd.to_numeric(df[col], errors='coerce').round(4)
    return df.reindex(columns=list(TSV_COLUMNS))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def add_arguments(ap: argparse.ArgumentParser) -> None:
    """
    Add ``maxiprot filter`` arguments to a parser.

    Parameters
    ----------
    ap : argparse.ArgumentParser
        Parser (top-level subparser or standalone) to populate.
    """
    ap.add_argument(
        'gff',
        nargs='?',
        default='-',
        help="Miniprot GFF3 with ##PAF lines (path or '-' for stdin).",
    )

    # Gating
    ap.add_argument(
        '--min-cov',
        type=float,
        default=0.60,
        help='Minimum reference coverage (fraction, 0-1).',
    )
    ap.add_argument(
        '--min-pid',
        type=float,
        default=0.30,
        help='Minimum AA identity (fraction, 0-1).',
    )
    ap.add_argument(
        '--strict',
        action='store_true',
        help='Drop a locus entirely if no candidate passes gates.',
    )

    # Scoring
    ap.add_argument(
        '--score-mode',
        choices=sorted(SCORE_MODES),
        default='ms_cov_pos',
        help='How candidates are ranked within a locus.',
    )
    ap.add_argument(
        '--length-metric',
        choices=['frac', 'aa'],
        default='frac',
        help=(
            "Length ingredient: 'frac' = aligned AA / reference length; "
            "'aa' = absolute aligned AA length."
        ),
    )
    ap.add_argument('--w-pid', type=float, default=1.0, help='Weight for identity.')
    ap.add_argument('--w-cov', type=float, default=1.0, help='Weight for coverage.')
    ap.add_argument(
        '--w-len', type=float, default=1.0, help='Weight for the length metric.'
    )
    ap.add_argument(
        '--w-pos', type=float, default=0.0, help='Weight for positives (linear/geom).'
    )
    ap.add_argument(
        '--w-ms', type=float, default=0.0, help='Weight for ms/qlen (linear/geom).'
    )
    ap.add_argument(
        '--w-AS', type=float, default=0.0, help='Weight for AS/qlen (linear/geom).'
    )

    # Selection policy
    ap.add_argument(
        '--selection-mode',
        choices=['best', 'prefer_intact', 'longest', 'longest_prefer_intact'],
        default='best',
        help=(
            'How to choose a single winner per locus: '
            "'best' by score; 'prefer_intact' picks best among intact if any; "
            "'longest' picks maximum AA length; "
            "'longest_prefer_intact' prefers intact among the longest."
        ),
    )

    # Emission policy
    ap.add_argument(
        '--emit-mode',
        choices=['best', 'all_passing'],
        default='best',
        help=(
            "'best' emits a single winner per locus; 'all_passing' emits all "
            'candidates that pass gates in each locus.'
        ),
    )
    ap.add_argument(
        '--max-per-locus',
        type=int,
        default=None,
        help='With --emit-mode=all_passing, cap emitted candidates per locus.',
    )

    # Locus clustering
    ap.add_argument(
        '--locus-pad',
        '--locus-gap',
        dest='locus_pad',
        type=int,
        default=5000,
        help='Max gap (nt) joining hits on the same target/strand into one locus.',
    )

    # Output / misc
    ap.add_argument(
        '--out-gff3',
        metavar='FILE',
        default='-',
        help="Write GFF3 to FILE, or '-' for stdout.",
    )
    ap.add_argument(
        '--out-tsv',
        metavar='FILE',
        default=None,
        help="Write TSV summary to FILE ('-' = stdout). Not written by default.",
    )
    ap.add_argument(
        '--id-prefix', default='PBM', help='Prefix for synthesized GFF3 IDs.'
    )
    ap.add_argument(
        '--log-level', default='INFO', choices=LOG_LEVELS, help='Logging level.'
    )


def build_arg_parser() -> argparse.ArgumentParser:
    """
    Build the standalone ``maxiprot filter`` parser.

    Returns
    -------
    argparse.ArgumentParser
        A configured argument parser instance.
    """
    ap = argparse.ArgumentParser(
        prog='maxiprot filter',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Pick the best (or all passing) miniprot alignments per locus.',
    )
    add_arguments(ap)
    return ap


def run(args: argparse.Namespace) -> int:
    """
    Run the filter pipeline for parsed arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments (see :func:`add_arguments`).

    Returns
    -------
    int
        Exit status: 0 success, 2 usage/input error.
    """
    init_logging(loglevel=args.log_level)

    for name, value in (('--min-cov', args.min_cov), ('--min-pid', args.min_pid)):
        if not 0.0 <= value <= 1.0:
            logger.error('%s must be a fraction in [0, 1]; got %s', name, value)
            return 2
    gff3_to_stdout = args.out_gff3 in (None, '-', '')
    if gff3_to_stdout and args.out_tsv == '-':
        logger.error(
            'GFF3 and TSV cannot both go to stdout; give --out-gff3 or '
            '--out-tsv a file path.'
        )
        return 2

    if args.gff in (None, '-', ''):
        logger.info('Reading GFF3 from stdin...')
    else:
        logger.info('Reading GFF3: %s', args.gff)
    try:
        paf, mrna_by_id, cds_by_parent = load_input(args.gff)
    except OSError as e:
        logger.error('Cannot read input: %s', e)
        return 2
    if paf.empty:
        logger.error('No PAF (##PAF) records were found; cannot proceed.')
        return 2

    paf = compute_alignment_metrics(paf)
    weights = Weights(
        pid=args.w_pid,
        cov=args.w_cov,
        len=args.w_len,
        pos=args.w_pos,
        ms=args.w_ms,
        AS=args.w_AS,
    )
    try:
        paf['score_raw'] = compute_scores(
            paf, args.score_mode, weights, args.length_metric
        )
    except ValueError as e:
        logger.error(str(e))
        return 2

    paf['passes'] = apply_gates(paf, args.min_cov, args.min_pid)
    paf['intact'] = gate_intact_v_pseudo(paf)
    paf = cluster_into_loci(paf, locus_pad=args.locus_pad)

    winners = select_per_locus(
        paf,
        selection_mode=args.selection_mode,
        emit_mode=args.emit_mode,
        max_per_locus=args.max_per_locus,
        strict=args.strict,
    )
    if winners.empty:
        logger.warning('No winners selected.')

    winners = attach_mrna_ids_to_winners(winners, mrna_by_id, cds_by_parent)

    tlens = (
        dict(zip(paf['tname'].astype(str), paf['tlen'].astype(int)))
        if not paf.empty
        else {}
    )
    write_best_gff3(
        args.out_gff3,
        winners,
        mrna_by_id,
        cds_by_parent,
        id_prefix=args.id_prefix,
        tlens=tlens,
    )

    if args.out_tsv is not None:
        tsv_df = build_tsv(winners, cds_by_parent)
        out, close_me = open_output(args.out_tsv)
        try:
            tsv_df.to_csv(out, sep='\t', index=False, na_rep='NA')
        finally:
            if close_me:
                out.close()
        if args.out_tsv != '-':
            logger.info('Wrote TSV summary to file: %s', args.out_tsv)

    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Run the ``maxiprot filter`` command-line interface.

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
