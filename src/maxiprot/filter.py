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
- Coverage gate is applied to *reference* coverage (fraction of reference AA
  covered by the mapped CDS), not target coverage.

"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from maxiprot.logs import init_logging

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def gff3_escape(x: str) -> str:
    """
    Escape a string for placement in a GFF3 attribute field.

    Parameters
    ----------
    x : str
        Raw attribute string.

    Returns
    -------
    str
        Escaped string according to GFF3 rules.
    """
    if x is None:
        return ''
    return (
        str(x)
        .replace('%', '%25')
        .replace(';', '%3B')
        .replace('=', '%3D')
        .replace(',', '%2C')
        .replace('\t', ' ')
    )


def _open_out_handle(path: Optional[str]):
    """
    Open an output handle for text writing.

    Parameters
    ----------
    path : str or None
        Path or '-', or None for stdout.

    Returns
    -------
    tuple
        (handle, close_required_bool)
    """
    if path in (None, '-', ''):
        return sys.stdout, False
    return open(path, 'w', encoding='utf-8'), True


# ---------------------------------------------------------------------------
# Parsing miniprot GFF3 (+ PAF header lines)
# ---------------------------------------------------------------------------


PAF_HEADER_RE = re.compile(r'^##PAF\t(.+)$')


def parse_gff3_with_paf(
    gff_text: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Parse a miniprot GFF3 text that contains ``##PAF`` header lines.

    Parameters
    ----------
    gff_text : str
        Raw GFF3 text.

    Returns
    -------
    tuple of (paf_df, mrna_df, cds_df)
        DataFrames holding parsed PAF headers and GFF mRNA/CDS rows.

    Notes
    -----
    - The PAF dataframe includes columns: ['qname','qlen','qs','qe','strand',
      'tname','tlen','ts','te','nmatch','alen','mapq', PAF tags expanded].
    - The mRNA and CDS tables are parsed from feature rows 'mRNA' and 'CDS'.
    """
    paf_records: List[Dict[str, object]] = []
    gff_rows: List[List[str]] = []

    for line in gff_text.splitlines():
        if not line or line.startswith('#') and not line.startswith('##PAF'):
            # retain normal comments only when writing output; here parse features
            if line.startswith('##PAF'):
                # handled below
                pass
            continue

        m = PAF_HEADER_RE.match(line)
        if m:
            payload = m.group(1).split('\t')
            # PAF core first 12 fields
            qname = payload[0]
            qlen = int(payload[1])
            qs = int(payload[2])
            qe = int(payload[3])
            strand = payload[4]
            tname = payload[5]
            tlen = int(payload[6])
            ts = int(payload[7])
            te = int(payload[8])
            nmatch = int(payload[9])  # matches or aligned matches
            alen = int(payload[10])  # aligned block length on query
            mapq = int(payload[11])

            tags = payload[12:]
            tag_dict: Dict[str, object] = {}
            for tag in tags:
                # e.g. AS:i:1234 or cg:Z:...
                try:
                    k, typ, val = tag.split(':', 2)
                except ValueError:
                    # robust to malformed
                    continue
                if typ == 'i':
                    try:
                        tag_dict[k] = int(val)
                    except ValueError:
                        tag_dict[k] = val
                else:
                    tag_dict[k] = val

            rec: Dict[str, object] = {
                'qname': qname,
                'qlen': qlen,
                'qs': qs,
                'qe': qe,
                'strand': strand,
                'tname': tname,
                'tlen': tlen,
                'ts': ts,
                'te': te,
                'nmatch': nmatch,
                'alen': alen,
                'mapq': mapq,
                **tag_dict,
            }
            paf_records.append(rec)
            continue

        # Otherwise try to parse a GFF3 feature row (tab-delimited 9 fields)
        if not line.startswith('#'):
            parts = line.rstrip('\n').split('\t')
            if len(parts) == 9:
                gff_rows.append(parts)

    paf = pd.DataFrame.from_records(paf_records)

    # Parse GFF features (mRNA, CDS) into dataframes
    cols = [
        'seqid',
        'source',
        'type',
        'start',
        'end',
        'score',
        'strand',
        'phase',
        'attributes',
    ]
    gff = (
        pd.DataFrame(gff_rows, columns=cols) if gff_rows else pd.DataFrame(columns=cols)
    )
    if not gff.empty:
        gff['start'] = gff['start'].astype(int)
        gff['end'] = gff['end'].astype(int)

    def parse_attr(attr: str) -> Dict[str, str]:
        d = {}
        for kv in attr.split(';'):
            if not kv:
                continue
            if '=' in kv:
                k, v = kv.split('=', 1)
                d[k] = v
        return d

    if not gff.empty:
        attrs = gff['attributes'].apply(parse_attr)
        gff = pd.concat(
            [gff.drop(columns=['attributes']), attrs.apply(pd.Series)], axis=1
        )

    mdf = (
        gff[gff['type'] == 'mRNA'].copy()
        if not gff.empty
        else pd.DataFrame(columns=gff.columns)
    )
    cdf = (
        gff[gff['type'] == 'CDS'].copy()
        if not gff.empty
        else pd.DataFrame(columns=gff.columns)
    )

    return paf, mdf, cdf


# ---------------------------------------------------------------------------
# Scoring & gating
# ---------------------------------------------------------------------------


def compute_alignment_metrics(paf: pd.DataFrame) -> pd.DataFrame:
    """
    Add derived metrics to a PAF dataframe.

    Parameters
    ----------
    paf : pandas.DataFrame
        Raw PAF dataframe parsed from '##PAF' lines.

    Returns
    -------
    pandas.DataFrame
        Copy of ``paf`` with added metric columns.

    Notes
    -----
    Adds columns:
    - ``cov_aa`` : float
        Coverage on query (reference protein), fraction of AA covered.
    - ``pid_aa`` : float
        Percent identity estimate (matches / aligned length).
    - ``cds_aa_len`` : int
        Inferred amino-acid length of mapped CDS (aligned length on query).
    - ``score_raw`` : float
        Default raw score derived from AS:i if present, otherwise ms:i, else
        number of matches ``nmatch``.
    """
    df = paf.copy()
    # coverage = aligned_length_on_query / qlen
    df['cov_aa'] = np.where(df['qlen'] > 0, df['alen'] / df['qlen'], 0.0)
    # pid: use ms:i matches if present else nmatch
    matches = np.where(
        df.get('ms', pd.Series([np.nan] * len(df))).notna(), df['ms'], df['nmatch']
    )
    df['pid_aa'] = np.where(df['alen'] > 0, matches / df['alen'], 0.0)
    # cds length in AA approx equals aligned AA length (alen)
    df['cds_aa_len'] = (df['alen']).astype(int)
    # score: prefer AS (alignment score), then ms, then nmatch
    df['score_raw'] = np.where(
        df.get('AS', pd.Series([np.nan] * len(df))).notna(),
        df['AS'],
        np.where(
            df.get('ms', pd.Series([np.nan] * len(df))).notna(), df['ms'], df['nmatch']
        ),
    )
    return df


def apply_gates(
    df: pd.DataFrame,
    min_cov: float,
    min_pid: float,
) -> pd.Series:
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

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain columns ``fs`` (frameshifts, int) and ``st`` (in-frame stops, int).

    Returns
    -------
    pandas.Series
        Boolean series True for intact, False for pseudo-like.
    """
    fs = df.get('fs', pd.Series([0] * len(df))).fillna(0).astype(int)
    st = df.get('st', pd.Series([0] * len(df))).fillna(0).astype(int)
    return (fs == 0) & (st == 0)


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
    mdf: pd.DataFrame,
    cdf: pd.DataFrame,
    overlap_thresh: float = 0.5,
) -> pd.DataFrame:
    """
    Attach mrna_id to winners by Target name and/or coordinate overlap.

    Strategy:
      1) Exact Target name == PAF qname on same seqid/strand (prefer).
      2) Otherwise, pick mRNA with max overlap with the PAF [tstart+1, tend] on same seqid/strand,
         accepting if overlap / min(len(paf), len(mrna_span)) >= overlap_thresh.

    Returns
    -------
    pandas.DataFrame
        Winners with 'mrna_id' column filled where matched.
    """
    winners = winners.copy()

    # Normalize PAF intervals to 1-based inclusive for matching.
    winners['paf_start'] = winners['ts'].astype(int) + 1
    winners['paf_end'] = winners['te'].astype(int)
    winners['paf_len'] = (winners['paf_end'] - winners['paf_start'] + 1).clip(lower=0)

    # Prepare mRNA table with a reliable ID + span (prefer union of CDS parts).
    mdf2 = mdf.copy()
    id_col = 'id' if 'id' in mdf2.columns else 'ID'  # your parser uses 'id'
    tgt_col = (
        'target'
        if 'target' in mdf2.columns
        else ('Target' if 'Target' in mdf2.columns else None)
    )

    # Compute union span of CDS per mRNA (safer than mRNA bounds when present)
    if not cdf.empty and 'parent' in cdf.columns:
        cds_span = (
            cdf.groupby('parent', as_index=False)
            .agg(start_min=('start', 'min'), end_max=('end', 'max'))
            .rename(columns={'parent': id_col})
        )
        mdf2 = mdf2.merge(cds_span, on=id_col, how='left')
        # choose CDS union if available, else mRNA start/end
        mdf2['mrna_span_start'] = np.where(
            mdf2['start_min'].notna(), mdf2['start_min'], mdf2['start']
        )
        mdf2['mrna_span_end'] = np.where(
            mdf2['end_max'].notna(), mdf2['end_max'], mdf2['end']
        )
        mdf2.drop(columns=['start_min', 'end_max'], inplace=True)
    else:
        mdf2['mrna_span_start'] = mdf2['start']
        mdf2['mrna_span_end'] = mdf2['end']

    # Parse Target name (first token)
    if tgt_col:
        mdf2['target_name'] = mdf2[tgt_col].apply(_parse_target_name)
    else:
        mdf2['target_name'] = None

    # Build fast lookups by (seqid,strand)
    by_chr_strand = {}
    for _, row in mdf2.iterrows():
        key = (row['seqid'], row['strand'])
        by_chr_strand.setdefault(key, []).append(row)

    # For each winner, try to assign mrna_id
    mrna_ids: list[Optional[str]] = []
    for _, w in winners.iterrows():
        key = (w['tname'], w['strand'])
        candidates = by_chr_strand.get(key, [])
        if not candidates:
            mrna_ids.append(None)
            continue

        # 1) target name match on same seqid/strand
        qname = str(w.get('qname', '')).strip() or None
        if qname:
            name_matches = [r for r in candidates if r.get('target_name') == qname]
        else:
            name_matches = []

        chosen: Optional[pd.Series] = None
        if name_matches:
            # If multiple, choose by largest overlap of span with PAF interval
            best = (-1, None)
            for r in name_matches:
                ov = _interval_overlap(
                    int(w['paf_start']),
                    int(w['paf_end']),
                    int(r['mrna_span_start']),
                    int(r['mrna_span_end']),
                )
                if ov > best[0]:
                    best = (ov, r)
            chosen = best[1]

        # 2) Fallback: coordinate overlap max (if above threshold)
        if chosen is None:
            best = (0.0, -1, None)  # (score, ov_len, row)
            for r in candidates:
                ov = _interval_overlap(
                    int(w['paf_start']),
                    int(w['paf_end']),
                    int(r['mrna_span_start']),
                    int(r['mrna_span_end']),
                )
                if ov <= 0:
                    continue
                paf_len = int(w['paf_len'])
                mrna_len = int(r['mrna_span_end']) - int(r['mrna_span_start']) + 1
                denom = max(1, min(paf_len, mrna_len))
                frac = ov / denom
                score = (frac, ov)
                if score > best[:2]:
                    best = (frac, ov, r)
            if best[2] is not None and best[0] >= overlap_thresh:
                chosen = best[2]

        mrna_ids.append(str(chosen[id_col]) if chosen is not None else None)

    winners['mrna_id'] = mrna_ids
    n_matched = sum(x is not None for x in mrna_ids)
    logging.info(
        'Linked %d/%d winners to existing mRNA features.', n_matched, len(winners)
    )
    return winners


# ---------------------------------------------------------------------------
# Locus clustering
# ---------------------------------------------------------------------------


def jaccard_1d(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    """
    Compute Jaccard index between two 1D intervals (closed).

    Parameters
    ----------
    a_start, a_end, b_start, b_end : int
        Interval endpoints (1-based closed coordinates).

    Returns
    -------
    float
        Jaccard = intersection / union.
    """
    a1, a2 = min(a_start, a_end), max(a_start, a_end)
    b1, b2 = min(b_start, b_end), max(b_start, b_end)
    inter = max(0, min(a2, b2) - max(a1, b1) + 1)
    union = (a2 - a1 + 1) + (b2 - b1 + 1) - inter
    return inter / union if union > 0 else 0.0


def cluster_into_loci(
    paf: pd.DataFrame,
    max_gap: int = 5000,
) -> pd.DataFrame:
    """
    Group alignments into loci on the same (tname, strand) within distance.

    Parameters
    ----------
    paf : pandas.DataFrame
        PAF-like dataframe (must have 'tname', 'strand', 'ts', 'te').
    max_gap : int, default 5000
        Non-overlapping candidates on the same strand within ``max_gap`` bases
        are grouped as the same locus.

    Returns
    -------
    pandas.DataFrame
        Copy of input paf with added ``locus`` label (int per (tname,strand) group).
    """
    df = paf.copy()
    if df.empty:
        df['locus'] = pd.Series(dtype=int)
        return df

    loci: List[int] = []
    next_id = 1

    # process per (tname, strand)
    for (_tname, _strand), sub in df.sort_values(
        ['tname', 'strand', 'ts', 'te']
    ).groupby(['tname', 'strand'], sort=False):
        current_locus = next_id
        last_end = None
        for _, row in sub.iterrows():
            start = min(int(row['ts']), int(row['te']))
            end = max(int(row['ts']), int(row['te']))
            if last_end is None:
                loci.append(current_locus)
                last_end = end
                continue
            # New locus if non-overlap and start - last_end > max_gap
            if start > last_end and (start - last_end) > max_gap:
                next_id += 1
                current_locus = next_id
            loci.append(current_locus)
            last_end = max(last_end, end)
        next_id += 1

    df = df.sort_values(['tname', 'strand', 'ts', 'te']).copy()
    df['locus'] = loci
    return df


# ---------------------------------------------------------------------------
# GFF writing
# ---------------------------------------------------------------------------


def write_best_gff3(
    out_path: Optional[str],
    winners: pd.DataFrame,
    mdf: pd.DataFrame,
    cdf: pd.DataFrame,
    id_prefix: str = 'PBM',
) -> None:
    """
    Write emitted candidates as valid GFF3.

    If ``out_path`` is None/'-' → write to stdout.

    Parameters
    ----------
    out_path : str or None
        Output destination (path or '-' / None for stdout).
    winners : pandas.DataFrame
        Rows corresponding to the emitted candidate(s); may be multiple rows
        per locus if ``--emit-mode all_passing`` was used.
    mdf : pandas.DataFrame
        DataFrame of parsed mRNA features from input GFF.
    cdf : pandas.DataFrame
        DataFrame of parsed CDS features from input GFF.
    id_prefix : str, default "PBM"
        Prefix used when synthesizing IDs.
    """

    if out_path and out_path != '-':
        logging.info('Writing filtered records to GFF3: %s', out_path)
    else:
        logging.info('Writing filtered records to stdout')

    out, close_me = _open_out_handle(out_path)
    try:
        out.write('##gff-version 3\n')

        # Pre-index mRNA / CDS by ID for fast retrieval
        mrna_by_id = (
            {row['ID']: row for _, row in mdf.iterrows()} if not mdf.empty else {}
        )
        cds_by_parent: Dict[str, List[pd.Series]] = {}
        if not cdf.empty:
            for _, row in cdf.iterrows():
                cds_by_parent.setdefault(row.get('Parent', ''), []).append(row)

        for _, r in winners.iterrows():
            # If this winner corresponds to an existing mRNA in the original file,
            # reuse it (and all its CDS parts), but ensure we synthesize a stable
            # gene parent.
            mrna_id = r.get('mrna_id', None)
            if pd.notna(mrna_id) and mrna_id in mrna_by_id:
                mrna = mrna_by_id[mrna_id]
                gene_id = f'{id_prefix}:gene:{gff3_escape(str(r["locus"]))}:{gff3_escape(str(r["qname"]))}'
                # gene feature (synthesized bounds from mRNA)
                gseq = mrna['seqid']
                gstart = int(mrna['start'])
                gend = int(mrna['end'])
                gstrand = mrna['strand']
                gscore = '.'
                gattr = f'ID={gene_id};Name={gff3_escape(r["qname"])};locus={gff3_escape(str(r["locus"]))}'
                out.write(
                    '\t'.join(
                        [
                            gseq,
                            'maxiprot',
                            'gene',
                            str(gstart),
                            str(gend),
                            gscore,
                            gstrand,
                            '.',
                            gattr,
                        ]
                    )
                    + '\n'
                )

                # write original mRNA, reparents to our gene
                m_attr_items = [f'ID={mrna["ID"]}', f'Parent={gene_id}']
                # keep Rank/Identity/etc if present on mrna row
                for keep in (
                    'Rank',
                    'Identity',
                    'Positive',
                    'Frameshift',
                    'StopCodon',
                    'Target',
                ):
                    if keep in mrna and pd.notna(mrna[keep]):
                        m_attr_items.append(f'{keep}={mrna[keep]}')
                m_attr = ';'.join(m_attr_items)
                out.write(
                    '\t'.join(
                        [
                            mrna['seqid'],
                            mrna['source'],
                            'mRNA',
                            str(int(mrna['start'])),
                            str(int(mrna['end'])),
                            '.' if pd.isna(mrna['score']) else str(mrna['score']),
                            mrna['strand'],
                            '.',
                            m_attr,
                        ]
                    )
                    + '\n'
                )

                # write CDS rows under this mRNA
                for cds_row in sorted(
                    cds_by_parent.get(mrna_id, []), key=lambda x: int(x['start'])
                ):
                    cds_attr_items = [f'ID={cds_row["ID"]}', f'Parent={mrna["ID"]}']
                    for keep in (
                        'Rank',
                        'Identity',
                        'Frameshift',
                        'StopCodon',
                        'Target',
                    ):
                        if keep in cds_row and pd.notna(cds_row[keep]):
                            cds_attr_items.append(f'{keep}={cds_row[keep]}')
                    cds_attr = ';'.join(cds_attr_items)
                    out.write(
                        '\t'.join(
                            [
                                cds_row['seqid'],
                                cds_row['source'],
                                'CDS',
                                str(int(cds_row['start'])),
                                str(int(cds_row['end'])),
                                '.'
                                if pd.isna(cds_row['score'])
                                else str(cds_row['score']),
                                cds_row['strand'],
                                str(cds_row['phase']),
                                cds_attr,
                            ]
                        )
                        + '\n'
                    )
                continue

            # Otherwise synthesize a minimal hierarchy from PAF
            logging.warning(
                f'No mrna_id found for {r["tname"]}, reconstructing range from PAF. CDS may be incorrect.'
            )
            seqid = r['tname']
            strand = r['strand']
            # mRNA bounds from ts/te
            mstart = min(int(r['ts']), int(r['te']))
            mend = max(int(r['ts']), int(r['te']))

            gene_id = f'{id_prefix}:gene:{gff3_escape(str(r["locus"]))}:{gff3_escape(str(r["qname"]))}'
            mrna_id = f'{id_prefix}:mrna:{gff3_escape(str(r["locus"]))}:{gff3_escape(str(r["qname"]))}'
            cds_id = f'{id_prefix}:cds:{gff3_escape(str(r["locus"]))}:{gff3_escape(str(r["qname"]))}'

            gattr = f'ID={gene_id};Name={gff3_escape(r["qname"])};locus={gff3_escape(str(r["locus"]))}'
            out.write(
                '\t'.join(
                    [
                        seqid,
                        'maxiprot',
                        'gene',
                        str(mstart),
                        str(mend),
                        '.',
                        strand,
                        '.',
                        gattr,
                    ]
                )
                + '\n'
            )

            mattr = f'ID={mrna_id};Parent={gene_id};Target={gff3_escape(r["qname"])} {int(r["qs"]) + 1} {int(r["qe"])}'
            out.write(
                '\t'.join(
                    [
                        seqid,
                        'maxiprot',
                        'mRNA',
                        str(mstart),
                        str(mend),
                        '.',
                        strand,
                        '.',
                        mattr,
                    ]
                )
                + '\n'
            )

            # A single synthetic CDS segment spanning mRNA bounds (phase set to 0)
            cattr = f'ID={cds_id};Parent={mrna_id}'
            out.write(
                '\t'.join(
                    [
                        seqid,
                        'maxiprot',
                        'CDS',
                        str(mstart),
                        str(mend),
                        '.',
                        strand,
                        '0',
                        cattr,
                    ]
                )
                + '\n'
            )

    finally:
        if close_me:
            out.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """
    Build the command-line interface parser for the ``maxiprot filter`` subcommand.

    The parser defines input handling, scoring and gating parameters, locus
    clustering, selection and emission policies, and output control flags.

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

    # Input
    ap.add_argument(
        '--gff',
        metavar='FILE',
        default='-',
        help='Miniprot GFF3 with ##PAF lines (path or "-" for stdin).',
    )

    # Gating
    ap.add_argument(
        '--min-cov',
        type=float,
        default=0.60,
        help='Minimum reference coverage (fraction).',
    )
    ap.add_argument(
        '--min-pid', type=float, default=0.30, help='Minimum AA identity (fraction).'
    )

    # Scoring
    ap.add_argument(
        '--score-mode',
        choices=['as_ms', 'ms_only', 'weighted', 'length_biased'],
        default='as_ms',
        help=(
            "Scoring mode. 'as_ms': prefer AS then ms; 'ms_only': use ms:i; "
            "'weighted': combine identity, coverage, and length; "
            "'length_biased': combine length with identity."
        ),
    )

    # Selection policy
    ap.add_argument(
        '--selection-mode',
        choices=['best', 'prefer_intact', 'longest', 'longest_prefer_intact'],
        default='best',
        help=(
            'How to choose a single winner per locus when emitting only one: '
            "'best' by score; 'prefer_intact' picks best among intact if any; "
            "'longest' picks maximum AA length; 'longest_prefer_intact' prefers intact among the longest."
        ),
    )

    # Emission policy
    ap.add_argument(
        '--emit-mode',
        choices=['best', 'all_passing'],
        default='best',
        help=(
            "Per-locus emission policy. 'best' (default) emits a single winner per locus. "
            "'all_passing' emits all candidates that pass gates in each locus "
            '(or none with --strict; otherwise falls back to the best single candidate).'
        ),
    )
    ap.add_argument(
        '--max-per-locus',
        type=int,
        default=None,
        help='When --emit-mode=all_passing, cap emitted candidates per locus (highest-scoring first).',
    )

    # Locus clustering
    ap.add_argument(
        '--max-gap',
        type=int,
        default=5000,
        help='Max gap to join non-overlapping on same strand into a locus.',
    )

    # Output / misc
    ap.add_argument(
        '--out-gff3',
        metavar='FILE',
        default='-',
        help='Write GFF3 to FILE or "-" for stdout.',
    )
    ap.add_argument(
        '--out-tsv',
        metavar='FILE',
        default=None,
        help='Write TSV summary to FILE; else write to stderr.',
    )
    ap.add_argument(
        '--id-prefix', default='PBM', help='Prefix for synthesized GFF3 IDs.'
    )
    ap.add_argument(
        '--strict',
        action='store_true',
        help='Drop a locus if no candidate passes gates.',
    )
    ap.add_argument('--log-level', default='INFO', help='Logging level.')

    return ap


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
    # Parse args
    args = build_arg_parser().parse_args(argv)

    # Set up logging
    init_logging(loglevel=args.log_level)

    # Read input
    if args.gff in (None, '-', ''):
        logging.info('Reading GFF3 from stdin...')
        gff_text = sys.stdin.read()
    else:
        logging.info('Reading GFF3: %s', args.gff)
        with open(args.gff, 'r', encoding='utf-8') as fh:
            gff_text = fh.read()

    # Parse
    paf, mdf, cdf = parse_gff3_with_paf(gff_text)
    if paf.empty:
        logging.error('No PAF (##PAF) records were found; cannot proceed.')
        return 2

    # Derive metrics, apply gates
    paf = compute_alignment_metrics(paf)

    # Scoring modes
    mode = args.score_mode
    if mode == 'as_ms':
        paf['score_raw'] = np.where(
            paf.get('AS', pd.Series([np.nan] * len(paf))).notna(),
            paf['AS'],
            np.where(
                paf.get('ms', pd.Series([np.nan] * len(paf))).notna(),
                paf['ms'],
                paf['nmatch'],
            ),
        )
    elif mode == 'ms_only':
        paf['score_raw'] = paf['ms'].fillna(0)
    elif mode == 'weighted':
        # simple weighted combination (identity 0.5, coverage 0.3, length 0.2)
        paf['score_raw'] = (
            0.5 * paf['pid_aa']
            + 0.3 * paf['cov_aa']
            + 0.2 * (paf['cds_aa_len'] / (paf['cds_aa_len'].max() or 1))
        )
    elif mode == 'length_biased':
        paf['score_raw'] = (
            0.6 * paf['cds_aa_len'] / (paf['cds_aa_len'].max() or 1)
        ) + 0.4 * paf['pid_aa']
    else:
        pass  # already set

    paf['passes'] = apply_gates(paf, args.min_cov, args.min_pid)
    paf['intact'] = gate_intact_v_pseudo(paf)

    # Locus clustering
    paf = cluster_into_loci(paf, max_gap=args.max_gap)

    # Per-locus selection/emission
    winners_rows: List[pd.Series] = []

    for locus, sub in paf.groupby('locus', sort=False):
        logging.info(
            '[Locus %s] %d candidate(s) on %s:%s',
            locus,
            len(sub),
            sub['tname'].iloc[0],
            sub['strand'].iloc[0],
        )
        # Detailed logging per candidate
        for _, r in sub.sort_values(
            ['score_raw', 'cds_aa_len'], ascending=[False, False]
        ).iterrows():
            logging.info(
                '  q=%s | score=%.5f cov=%.3f pid=%.3f lenAA=%d fs=%d st=%d passes=%s intact=%s',
                r['qname'],
                float(r['score_raw']),
                float(r['cov_aa']),
                float(r['pid_aa']),
                int(r['cds_aa_len']),
                int(r.get('fs', 0)),
                int(r.get('st', 0)),
                'Y' if r['passes'] else 'N',
                'Y' if r['intact'] else 'N',
            )

        if args.strict and not any(sub['passes']):
            logging.info(
                '[Locus %s] no candidates pass gates -> skipped due to --strict', locus
            )
            continue

        # Emission policy
        if args.emit_mode == 'all_passing':
            emit_df = sub.loc[sub['passes']].copy()
            if emit_df.empty:
                if args.strict:
                    logging.info(
                        '[Locus %s] no candidates pass gates -> skipped due to --strict',
                        locus,
                    )
                    continue
                # Non-strict: still emit the single best candidate
                emit_df = sub.iloc[[0]].copy()
            # Highest-scoring first, then longest as tie-breaker
            emit_df = emit_df.sort_values(
                ['score_raw', 'cds_aa_len'], ascending=[False, False]
            )
            if args.max_per_locus is not None and args.max_per_locus > 0:
                emit_df = emit_df.head(args.max_per_locus)
            for _, pick in emit_df.iterrows():
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
            continue

        # Selection policy (single best)
        def choose(df: pd.DataFrame) -> pd.Series:
            if args.selection_mode in ('best', 'prefer_intact'):
                # rank by score then length
                return df.sort_values(
                    ['score_raw', 'cds_aa_len'], ascending=[False, False]
                ).iloc[0]
            elif args.selection_mode in ('longest', 'longest_prefer_intact'):
                return df.sort_values(
                    ['cds_aa_len', 'score_raw'], ascending=[False, False]
                ).iloc[0]
            else:
                return df.sort_values(
                    ['score_raw', 'cds_aa_len'], ascending=[False, False]
                ).iloc[0]

        working = sub.copy()
        if args.selection_mode in ('prefer_intact', 'longest_prefer_intact'):
            intact = working[working['intact']]
            if not intact.empty:
                working = intact

        pick = None
        # In strict mode, pick must pass gates; otherwise, pick from working regardless
        if args.strict:
            passers = working[working['passes']]
            if passers.empty:
                logging.info(
                    '[Locus %s] all candidates fail gates -> skipped due to --strict',
                    locus,
                )
                continue
            pick = choose(passers)
        else:
            # prefer a passer if any
            passers = working[working['passes']]
            pick = choose(passers if not passers.empty else working)

        logging.info(
            '[Locus %s] SELECTED q=%s (score=%.5f, cov=%.3f, pid=%.3f, lenAA=%d, fs=%d, st=%d)',
            locus,
            pick['qname'],
            float(pick['score_raw']),
            float(pick['cov_aa']),
            float(pick['pid_aa']),
            int(pick['cds_aa_len']),
            int(pick.get('fs', 0)),
            int(pick.get('st', 0)),
        )
        winners_rows.append(pick)

    if not winners_rows:
        logging.warning('No winners selected. Exiting.')
        return 0

    winners = pd.DataFrame(winners_rows).reset_index(drop=True)

    # Try to map emitted rows back to original mrna/cds IDs if possible (heuristic)
    # Here the input GFF mRNA rows usually include Target=... pointing back to qname
    if not mdf.empty:
        mdf_idx = mdf.reset_index(drop=True).copy()
        mdf_idx['mrna_id'] = mdf_idx['ID']
        # Join winners to mRNA by (Target contains qname) OR by coords overlap on same seqid/strand
        # Best-effort: users with robust miniprot GFF should get this mapping correct.
        winners = attach_mrna_ids_to_winners(winners, mdf, cdf)

    # Write GFF3 (stdout by default)
    write_best_gff3(args.out_gff3, winners, mdf, cdf, id_prefix=args.id_prefix)

    # Write TSV summary (stderr by default unless out-tsv is set)
    tsv_df = winners[
        [
            'tname',
            'strand',
            'locus',
            'qname',
            'cov_aa',
            'pid_aa',
            'cds_aa_len',
            'score_raw',
            'fs',
            'st',
        ]
    ].copy()
    tsv_df.rename(
        columns={
            'tname': 'seqid',
            'cds_aa_len': 'len_aa',
            'score_raw': 'score',
        },
        inplace=True,
    )
    tsv_text = tsv_df.to_csv(sep='\t', index=False)

    if args.out_tsv and args.out_tsv != '-':
        with open(args.out_tsv, 'w', encoding='utf-8') as fh:
            fh.write(tsv_text)
        logging.info('Wrote TSV summary to file: %s', args.out_tsv)
    else:
        # Pure TSV to stderr (logs may interleave; consider --log-level ERROR or redirect)
        logging.info('Writing TSV summary to stderr:')
        sys.stderr.write(tsv_text)

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
