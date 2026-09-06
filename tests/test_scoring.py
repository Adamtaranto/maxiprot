from __future__ import annotations

import math

import pandas as pd
import pytest

from maxiprot.filter import cluster_into_loci
from maxiprot.scoring import (
    SCORE_MODES,
    Weights,
    compute_alignment_metrics,
    compute_scores,
    parse_cs,
)


# -----------------------------
# cs:Z parser
# -----------------------------


def test_parse_cs_matches_only():
    s = parse_cs(':100')
    assert s.identical == 100
    assert s.substituted == 0
    assert s.aligned_aa == 100


def test_parse_cs_mixed_ops():
    # 50 matches, 2 substitutions, an insertion, a deletion, an intron
    s = parse_cs(':50*ab*cd+eff-gg~gt250ag:25')
    assert s.identical == 75
    assert s.substituted == 2
    assert s.aligned_aa == 77


def test_parse_cs_unknown_op_skipped():
    s = parse_cs(':10^:5')
    assert s.identical == 15


# -----------------------------
# Metrics
# -----------------------------


def _paf_df(**overrides) -> pd.DataFrame:
    row = {
        'qname': 'Q',
        'qlen': 100,
        'qs': 0,
        'qe': 100,
        'strand': '+',
        'tname': 'chr1',
        'tlen': 10_000,
        'ts': 0,
        'te': 300,
        'nmatch': 80,
        'alen': 100,
        'mapq': 60,
        'AS': 300,
        'ms': 290,
        'np': 90,
        'fs': 0,
        'st': 0,
        'cs': ':80' + '*aa' * 20,
    }
    row.update(overrides)
    return pd.DataFrame([row])


def test_metrics_from_cs():
    df = compute_alignment_metrics(_paf_df())
    r = df.iloc[0]
    assert r['aligned_aa'] == 100
    assert r['pid_aa'] == pytest.approx(0.8)
    assert r['cov_aa'] == pytest.approx(1.0)
    assert r['positives'] == pytest.approx(0.9)
    assert r['ms_per_qlen'] == pytest.approx(2.9)
    assert r['AS_per_qlen'] == pytest.approx(3.0)
    assert r['cds_aa_len'] == 100
    assert r['len_frac'] == pytest.approx(1.0)


def test_metrics_bounds_on_gappy_alignment():
    """alen > qe-qs (gaps) must not push cov or pid above 1."""
    df = compute_alignment_metrics(_paf_df(alen=250, cs=None, qs=10, qe=90, nmatch=300))
    r = df.iloc[0]
    assert 0.0 <= r['pid_aa'] <= 1.0
    assert 0.0 <= r['cov_aa'] <= 1.0
    assert r['cov_aa'] == pytest.approx(0.8)


def test_metrics_fallbacks_without_tags():
    """No cs/np/ms/AS tags: fallbacks apply, nothing raises."""
    df = _paf_df()
    df = df.drop(columns=['cs', 'np', 'ms', 'AS'])
    out = compute_alignment_metrics(df)
    r = out.iloc[0]
    assert r['pid_aa'] == pytest.approx(0.8)  # nmatch/alen
    assert r['positives'] == 0.0
    assert r['ms_per_qlen'] == pytest.approx(0.8)  # nmatch fallback
    assert r['AS_per_qlen'] == pytest.approx(0.8)


def test_metrics_zero_qlen_is_safe():
    df = compute_alignment_metrics(_paf_df(qlen=0))
    r = df.iloc[0]
    assert r['cov_aa'] == 0.0
    assert r['ms_per_qlen'] == 0.0


# -----------------------------
# Score modes (hand-computed)
# -----------------------------


@pytest.fixture
def metrics_df() -> pd.DataFrame:
    return compute_alignment_metrics(_paf_df())


def test_score_ms_cov_pos(metrics_df):
    got = compute_scores(metrics_df, 'ms_cov_pos', Weights())
    # ms/qlen * cov * (0.5 + 0.5*positives) = 2.9 * 1.0 * 0.95
    assert got.iloc[0] == pytest.approx(2.9 * 0.95)


def test_score_AS_and_ms(metrics_df):
    assert compute_scores(metrics_df, 'AS', Weights()).iloc[0] == pytest.approx(3.0)
    assert compute_scores(metrics_df, 'ms', Weights()).iloc[0] == pytest.approx(2.9)


def test_score_pid_cov(metrics_df):
    assert compute_scores(metrics_df, 'pid_cov', Weights()).iloc[0] == pytest.approx(
        0.8
    )


def test_score_pid_cov_len_weights(metrics_df):
    w = Weights(pid=2.0, cov=1.0, len=1.0)
    got = compute_scores(metrics_df, 'pid_cov_len', w, 'frac')
    assert got.iloc[0] == pytest.approx(0.8**2 * 1.0 * 1.0)


def test_score_length_metrics(metrics_df):
    assert compute_scores(metrics_df, 'length', Weights(), 'frac').iloc[
        0
    ] == pytest.approx(1.0)
    assert compute_scores(metrics_df, 'length', Weights(), 'aa').iloc[
        0
    ] == pytest.approx(100.0)


def test_score_linear(metrics_df):
    w = Weights(pid=1.0, cov=1.0, len=1.0, pos=1.0, ms=1.0, AS=1.0)
    got = compute_scores(metrics_df, 'linear', w, 'frac')
    assert got.iloc[0] == pytest.approx(0.8 + 1.0 + 1.0 + 0.9 + 2.9 + 3.0)


def test_score_geom(metrics_df):
    w = Weights(pid=1.0, cov=1.0, len=1.0, pos=0.0, ms=0.0, AS=0.0)
    got = compute_scores(metrics_df, 'geom', w, 'frac')
    assert got.iloc[0] == pytest.approx(0.8 * 1.0 * 1.0)


def test_score_geom_all_zero_weights_raises(metrics_df):
    w = Weights(pid=0.0, cov=0.0, len=0.0, pos=0.0, ms=0.0, AS=0.0)
    with pytest.raises(ValueError):
        compute_scores(metrics_df, 'geom', w, 'frac')


def test_score_geom_zero_factor_is_clamped(metrics_df):
    df = metrics_df.copy()
    df.loc[df.index[0], 'pid_aa'] = 0.0
    got = compute_scores(df, 'geom', Weights(), 'frac')
    assert math.isfinite(got.iloc[0])
    assert got.iloc[0] >= 0.0


def test_all_modes_registered():
    assert set(SCORE_MODES) == {
        'ms_cov_pos',
        'AS',
        'ms',
        'pid_cov',
        'pid_cov_len',
        'length',
        'linear',
        'geom',
    }


# -----------------------------
# Locus clustering
# -----------------------------


def _hits(rows):
    return pd.DataFrame(rows, columns=['tname', 'strand', 'ts', 'te']).assign(qname='Q')


def test_cluster_empty():
    df = cluster_into_loci(_hits([]))
    assert 'locus' in df.columns
    assert df.empty


def test_cluster_gap_boundary():
    """Gap exactly == locus_pad joins; one more base splits."""
    df = cluster_into_loci(
        _hits([('c', '+', 0, 100), ('c', '+', 1100, 1200)]), locus_pad=1000
    )
    assert df['locus'].nunique() == 1
    df = cluster_into_loci(
        _hits([('c', '+', 0, 100), ('c', '+', 1101, 1200)]), locus_pad=1000
    )
    assert df['locus'].nunique() == 2


def test_cluster_overlap_always_joins():
    df = cluster_into_loci(
        _hits([('c', '+', 0, 5000), ('c', '+', 4000, 9000)]), locus_pad=10
    )
    assert df['locus'].nunique() == 1


def test_cluster_strands_are_independent():
    df = cluster_into_loci(
        _hits([('c', '+', 0, 100), ('c', '-', 0, 100)]), locus_pad=1000
    )
    assert df['locus'].nunique() == 2
    assert set(df['locus']) == {'c:+:1', 'c:-:1'}


def test_cluster_tandem_chain_merges():
    """Single-linkage: tandem copies < locus_pad apart chain into one locus."""
    df = cluster_into_loci(
        _hits(
            [
                ('c', '+', 0, 1000),
                ('c', '+', 4000, 5000),
                ('c', '+', 8000, 9000),
            ]
        ),
        locus_pad=5000,
    )
    assert df['locus'].nunique() == 1


def test_cluster_locus_label_format():
    df = cluster_into_loci(
        _hits([('c1', '+', 0, 100), ('c1', '+', 50_000, 50_100), ('c2', '-', 0, 9)]),
        locus_pad=1000,
    )
    assert set(df['locus']) == {'c1:+:1', 'c1:+:2', 'c2:-:1'}
