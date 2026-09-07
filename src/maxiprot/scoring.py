"""Alignment metrics and scoring modes for ``maxiprot filter``.

Implements the documented scoring model: per-candidate metrics derived from
miniprot ``##PAF`` records (identity from ``cs:Z``, query coverage, positives
from ``np:i``, per-length ``ms``/``AS``) and the eight ranking modes
``ms_cov_pos``, ``AS``, ``ms``, ``pid_cov``, ``pid_cov_len``, ``length``,
``linear`` and ``geom``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import logging
import re

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

#: Minimum factor value used in the geometric score so that zero/negative
#: ingredients do not produce NaN/inf. Documented behavior.
GEOM_FLOOR = 1e-9


@dataclass(frozen=True)
class CsStats:
    """Counts derived from a miniprot ``cs:Z`` alignment string."""

    identical: int
    substituted: int

    @property
    def aligned_aa(self) -> int:
        """Number of aligned AA columns (matches + substitutions)."""
        return self.identical + self.substituted


# miniprot cs:Z tokens. The grammar extends minimap2's cs with introns and
# frameshifts; ops that do not represent aligned AA columns contribute 0.
_CS_TOKEN = re.compile(
    r'(:\d+)'  # :N     -> N identical residues
    r'|(\*[a-zA-Z]+)'  # *xxY   -> one substituted position
    r'|(\+[a-zA-Z]+)'  # +seq   -> insertion in query (unaligned)
    r'|(-[a-zA-Z]+)'  # -seq   -> deletion from query (unaligned)
    r'|(~[a-zA-Z]{2}\d+[a-zA-Z]{2})'  # ~gt123ag -> intron
    r'|(\$\d*)'  # frameshift-style ops
    r'|(\d+)'  # bare digits (miniprot frameshift lengths)
)

_warned_unknown_cs_ops: set[str] = set()


def parse_cs(cs: str) -> CsStats:
    """
    Parse a miniprot ``cs:Z`` string into identity/substitution counts.

    Unknown operators are skipped; each distinct unknown operator triggers a
    single warning per process.

    Parameters
    ----------
    cs : str
        The ``cs:Z`` tag value (without the ``cs:Z:`` prefix).

    Returns
    -------
    CsStats
        Identity and substitution counts; ``aligned_aa`` is their sum.
    """
    identical = 0
    substituted = 0
    pos = 0
    n = len(cs)
    while pos < n:
        m = _CS_TOKEN.match(cs, pos)
        if m is None:
            op = cs[pos]
            if op not in _warned_unknown_cs_ops:
                _warned_unknown_cs_ops.add(op)
                logger.warning(
                    "Unknown cs:Z operator '%s' encountered; skipping "
                    '(further occurrences are silent)',
                    op,
                )
            pos += 1
            continue
        token = m.group(0)
        if token.startswith(':'):
            identical += int(token[1:])
        elif token.startswith('*'):
            substituted += 1
        pos = m.end()
    return CsStats(identical=identical, substituted=substituted)


@dataclass(frozen=True)
class Weights:
    """Scoring weights (exponents for multiplicative modes, coefficients for ``linear``)."""

    pid: float = 1.0
    cov: float = 1.0
    len: float = 1.0
    pos: float = 0.0
    ms: float = 0.0
    AS: float = 0.0


def _series_or_nan(df: pd.DataFrame, col: str) -> pd.Series:
    """Return ``df[col]`` as float, or an all-NaN Series if absent."""
    if col in df.columns:
        return pd.to_numeric(df[col], errors='coerce')
    return pd.Series(np.nan, index=df.index, dtype=float)


def compute_alignment_metrics(paf: pd.DataFrame) -> pd.DataFrame:
    """
    Add derived metric columns to a PAF dataframe.

    Parameters
    ----------
    paf : pandas.DataFrame
        Parsed ``##PAF`` records (one row per alignment). Requires the core
        columns ``qname, qlen, qs, qe, tname, ts, te, nmatch, alen, mapq``;
        tag columns (``AS``, ``ms``, ``np``, ``fs``, ``st``, ``cs``) are
        optional.

    Returns
    -------
    pandas.DataFrame
        Copy of ``paf`` with added columns ``aligned_aa``, ``pid_aa``,
        ``cov_aa``, ``positives``, ``ms_per_qlen``, ``AS_per_qlen``,
        ``cds_aa_len``, ``len_frac``.

    Notes
    -----
    - ``pid_aa`` = identical AA / aligned AA columns from ``cs:Z``;
      fallback ``nmatch / alen`` when ``cs:Z`` is absent.
    - ``cov_aa`` = ``(qe - qs) / qlen`` (query = reference protein).
    - ``positives`` = ``np:i / aligned_aa`` (0.0 when ``np:i`` absent).
    - ``ms_per_qlen`` / ``AS_per_qlen`` fall back along ``AS -> ms -> nmatch``
      when a tag is missing (warned once per input).
    """
    df = paf.copy()
    qlen = pd.to_numeric(df['qlen'], errors='coerce')
    qspan = pd.to_numeric(df['qe'], errors='coerce') - pd.to_numeric(
        df['qs'], errors='coerce'
    )
    nmatch = pd.to_numeric(df['nmatch'], errors='coerce')
    alen = pd.to_numeric(df['alen'], errors='coerce')

    # cs:Z-derived identity
    if 'cs' in df.columns:
        cs_col = df['cs']
        has_cs = cs_col.notna()
    else:
        cs_col = pd.Series(None, index=df.index, dtype=object)
        has_cs = pd.Series(False, index=df.index)

    identical = pd.Series(np.nan, index=df.index, dtype=float)
    aligned = pd.Series(np.nan, index=df.index, dtype=float)
    if has_cs.any():
        stats = cs_col[has_cs].astype(str).map(parse_cs)
        identical.loc[has_cs] = stats.map(lambda s: s.identical).astype(float)
        aligned.loc[has_cs] = stats.map(lambda s: s.aligned_aa).astype(float)
    if not has_cs.all():
        logger.warning(
            '%d/%d PAF records lack a cs:Z tag; using span-based fallbacks '
            'for identity and aligned length.',
            int((~has_cs).sum()),
            len(df),
        )
    aligned = aligned.where(has_cs, qspan)
    df['aligned_aa'] = aligned

    pid_cs = (identical / aligned.replace(0, np.nan)).clip(0.0, 1.0)
    pid_fallback = (nmatch / alen.replace(0, np.nan)).clip(0.0, 1.0)
    df['pid_aa'] = pid_cs.where(has_cs, pid_fallback).fillna(0.0)

    df['cov_aa'] = (qspan / qlen.replace(0, np.nan)).clip(0.0, 1.0).fillna(0.0)

    np_tag = _series_or_nan(df, 'np')
    if np_tag.isna().all() and len(df):
        logger.warning("No np:i tags found; 'positives' is 0.0 for all candidates.")
    df['positives'] = (np_tag / aligned.replace(0, np.nan)).clip(0.0, 1.0).fillna(0.0)

    ms_tag = _series_or_nan(df, 'ms')
    as_tag = _series_or_nan(df, 'AS')
    if ms_tag.isna().any() and len(df):
        logger.warning(
            '%d/%d PAF records lack an ms:i tag; falling back to nmatch.',
            int(ms_tag.isna().sum()),
            len(df),
        )
    if as_tag.isna().any() and len(df):
        logger.warning(
            '%d/%d PAF records lack an AS:i tag; falling back to ms/nmatch.',
            int(as_tag.isna().sum()),
            len(df),
        )
    ms_eff = ms_tag.where(ms_tag.notna(), nmatch)
    as_eff = as_tag.where(as_tag.notna(), ms_eff)
    safe_qlen = qlen.replace(0, np.nan)
    df['ms_per_qlen'] = (ms_eff / safe_qlen).fillna(0.0)
    df['AS_per_qlen'] = (as_eff / safe_qlen).fillna(0.0)

    df['cds_aa_len'] = qspan.fillna(0).astype(int)
    df['len_frac'] = (aligned / safe_qlen).clip(0.0, 1.0).fillna(0.0)
    return df


def _length_metric(df: pd.DataFrame, length_metric: str) -> pd.Series:
    """Return the chosen length ingredient (``frac`` or ``aa``)."""
    if length_metric == 'aa':
        return df['cds_aa_len'].astype(float)
    return df['len_frac']


def _score_ms_cov_pos(df: pd.DataFrame, w: Weights, L: pd.Series) -> pd.Series:
    return df['ms_per_qlen'] * df['cov_aa'] * (0.5 + 0.5 * df['positives'])


def _score_as(df: pd.DataFrame, w: Weights, L: pd.Series) -> pd.Series:
    return df['AS_per_qlen']


def _score_ms(df: pd.DataFrame, w: Weights, L: pd.Series) -> pd.Series:
    return df['ms_per_qlen']


def _score_pid_cov(df: pd.DataFrame, w: Weights, L: pd.Series) -> pd.Series:
    return df['pid_aa'] * df['cov_aa']


def _score_pid_cov_len(df: pd.DataFrame, w: Weights, L: pd.Series) -> pd.Series:
    return (
        df['pid_aa'].clip(lower=GEOM_FLOOR) ** w.pid
        * df['cov_aa'].clip(lower=GEOM_FLOOR) ** w.cov
        * L.clip(lower=GEOM_FLOOR) ** w.len
    )


def _score_length(df: pd.DataFrame, w: Weights, L: pd.Series) -> pd.Series:
    return L


def _score_linear(df: pd.DataFrame, w: Weights, L: pd.Series) -> pd.Series:
    return (
        w.pid * df['pid_aa']
        + w.cov * df['cov_aa']
        + w.pos * df['positives']
        + w.ms * df['ms_per_qlen']
        + w.AS * df['AS_per_qlen']
        + w.len * L
    )


def _score_geom(df: pd.DataFrame, w: Weights, L: pd.Series) -> pd.Series:
    factors = [
        (w.pid, df['pid_aa']),
        (w.cov, df['cov_aa']),
        (w.pos, df['positives']),
        (w.ms, df['ms_per_qlen']),
        (w.AS, df['AS_per_qlen']),
        (w.len, L),
    ]
    active = [(wt, s) for wt, s in factors if wt > 0]
    if not active:
        raise ValueError(
            "score-mode 'geom' requires at least one positive weight "
            '(--w-pid/--w-cov/--w-len/--w-pos/--w-ms/--w-AS)'
        )
    score = pd.Series(1.0, index=df.index)
    for wt, s in active:
        score = score * (s.clip(lower=GEOM_FLOOR) ** wt)
    return score


ScoreFn = Callable[[pd.DataFrame, Weights, pd.Series], pd.Series]

SCORE_MODES: dict[str, ScoreFn] = {
    'ms_cov_pos': _score_ms_cov_pos,
    'AS': _score_as,
    'ms': _score_ms,
    'pid_cov': _score_pid_cov,
    'pid_cov_len': _score_pid_cov_len,
    'length': _score_length,
    'linear': _score_linear,
    'geom': _score_geom,
}


def compute_scores(
    df: pd.DataFrame,
    mode: str,
    weights: Weights,
    length_metric: str = 'frac',
) -> pd.Series:
    """
    Compute the per-candidate ranking score for a scoring mode.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of :func:`compute_alignment_metrics`.
    mode : str
        One of :data:`SCORE_MODES`.
    weights : Weights
        Scoring weights.
    length_metric : str
        ``'frac'`` (aligned fraction of the reference) or ``'aa'``
        (absolute aligned AA length).

    Returns
    -------
    pandas.Series
        ``score_raw`` values aligned to ``df.index``.

    Raises
    ------
    KeyError
        If ``mode`` is not a known scoring mode.
    ValueError
        If ``geom`` is requested with all weights zero.
    """
    fn = SCORE_MODES[mode]
    L = _length_metric(df, length_metric)
    return fn(df, weights, L)
