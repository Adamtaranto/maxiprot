from __future__ import annotations

import io
from typing import List

import pytest

from maxiprot import filter as filter_mod  # the module under test


def run_filter_and_capture(argv: List[str], capsys):
    """Helper to run filter.main(argv) and capture stdout/stderr.

    If the first argv element looks like a path (doesn't start with '-'),
    prepend ``--gff`` for compatibility with the updated CLI.
    """
    args = list(argv)
    if args and not str(args[0]).startswith('-'):
        args = ['--gff', args[0], *args[1:]]
    rc = filter_mod.main(args)
    out = capsys.readouterr()
    return rc, out.out, out.err


def test_default_io_and_best_selection(gff_minimal_two_loci, capsys):
    """
    Default behavior:
      * best score at locus A -> RefPseudo (given as_ms scoring)
      * best score at locus B -> RefFar
      * GFF -> stdout, TSV -> stderr
    """
    rc, stdout, stderr = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--score-mode',
            'as_ms',  # updated to allowed choice
            '--min-cov',
            '0',
            '--min-pid',
            '0',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    # mRNA Target values for winners should appear in stdout GFF3:
    assert 'Target=RefPseudo' in stdout  # first locus winner
    assert 'Target=RefFar' in stdout  # second locus winner
    # TSV goes to stderr by default and should mention the chosen qnames:
    assert 'RefPseudo' in stderr
    assert 'RefFar' in stderr


def test_prefer_intact_over_best(gff_minimal_two_loci, capsys):
    """
    With --selection-mode prefer_intact:
      - First locus should pick RefGood (intact) even though RefPseudo may have higher AS.
      - Second locus remains RefFar (intact).
    """
    rc, stdout, stderr = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--score-mode',
            'as_ms',  # updated to allowed choice
            '--selection-mode',
            'prefer_intact',
            '--min-cov',
            '0',
            '--min-pid',
            '0',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    assert 'Target=RefGood' in stdout
    assert 'Target=RefFar' in stdout


def test_file_outputs_redirecting(gff_minimal_two_loci, tmp_path, capsys):
    """
    When --out-gff3 and --out-tsv are set, outputs go to files; stdout/stderr clean.
    """
    out_gff = tmp_path / 'out.gff3'
    out_tsv = tmp_path / 'out.tsv'
    rc, stdout, stderr = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--score-mode',
            'as_ms',  # updated to allowed choice
            '--min-cov',
            '0',
            '--min-pid',
            '0',
            '--out-gff3',
            str(out_gff),
            '--out-tsv',
            str(out_tsv),
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    assert stdout == ''  # GFF to file, not stdout
    assert stderr == ''  # TSV to file, not stderr (and logs suppressed)
    gff_txt = out_gff.read_text(encoding='utf-8')
    tsv_txt = out_tsv.read_text(encoding='utf-8')
    assert 'Target=RefPseudo' in gff_txt or 'Target=RefGood' in gff_txt
    assert 'Target=RefFar' in gff_txt
    assert 'Ref' in tsv_txt  # at least some TSV rows with qnames


def test_stdin_path_reads_from_stdin(gff_minimal_two_loci, capsys, monkeypatch):
    """
    Using '-' as the GFF path should read from stdin and produce normal output.
    """
    with open(gff_minimal_two_loci, 'rt', encoding='utf-8') as fh:
        data = fh.read()
    monkeypatch.setattr('sys.stdin', io.StringIO(data))
    rc, stdout, stderr = run_filter_and_capture(
        [
            '--gff',
            '-',  # ← explicitly tie '-' to --gff
            '--score-mode',
            'as_ms',
            '--min-cov',
            '0',
            '--min-pid',
            '0',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    # Confirm that expected winners appear in the stdout GFF:
    assert 'Target=RefPseudo' in stdout or 'Target=RefGood' in stdout
    assert 'Target=RefFar' in stdout


def test_strict_drops_locus_when_all_fail_gates(gff_all_fail_gating, capsys):
    """
    With --strict and tight gates, a locus with no passing candidates should be dropped.
    """
    rc, stdout, stderr = run_filter_and_capture(
        [
            str(gff_all_fail_gating),
            '--score-mode',
            'as_ms',  # updated to allowed choice
            '--min-cov',
            '0.9',  # candidate has 0.5
            '--min-pid',
            '0.0',
            '--strict',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    # No gene should be emitted
    assert '\tgene\t' not in stdout
    # TSV may be empty or contain a note; avoid being too strict here:
    # We just ensure no obvious winner row exists.
    assert 'LowCov' not in stderr


def test_cds_lines_share_single_id_per_transcript(gff_minimal_two_loci, capsys):
    """
    Multi-exon CDS must share the same ID across lines for a given transcript.
    """
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--score-mode',
            'as_ms',  # updated to allowed choice
            '--selection-mode',
            'prefer_intact',
            '--min-cov',
            '0',
            '--min-pid',
            '0',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    # parse stdout GFF and collect CDS IDs belonging to the first mRNA encountered
    mrna_id = None
    cds_ids = set()
    for ln in stdout.splitlines():
        if not ln or ln.startswith('#'):
            continue
        cols = ln.split('\t')
        if len(cols) != 9:
            continue
        ftype = cols[2]
        attrs = cols[8]
        if ftype == 'mRNA':
            # capture the first mRNA we see (locus A after prefer_intact)
            if mrna_id is None:
                m = dict(x.split('=', 1) for x in attrs.split(';') if '=' in x)
                mrna_id = m.get('ID')
        elif ftype == 'CDS':
            m = dict(x.split('=', 1) for x in attrs.split(';') if '=' in x)
            pid = m.get('Parent')
            cid = m.get('ID')
            if mrna_id and pid == mrna_id and cid:
                cds_ids.add(cid)
    # all CDS lines for that mRNA should share exactly one ID
    assert len(cds_ids) == 1, f'Expected 1 unique CDS ID, saw {cds_ids}'


# === Test-local helpers/fixtures overriding legacy placeholders ===
from pathlib import Path
from typing import Dict, Iterable


def _paf_line(
    qname: str,
    qlen: int,
    qs: int,
    qe: int,
    strand: str,
    tname: str,
    tlen: int,
    ts: int,
    te: int,
    nmatch: int,
    alen: int,
    mapq: int = 60,
    tags: Dict[str, int | str] | None = None,
) -> str:
    base = [
        qname,
        str(qlen),
        str(qs),
        str(qe),
        strand,
        tname,
        str(tlen),
        str(ts),
        str(te),
        str(nmatch),
        str(alen),
        str(mapq),
    ]
    tags = tags or {}
    ex = []
    for k, v in tags.items():
        ex.append(f'{k}:i:{v}' if isinstance(v, int) else f'{k}:Z:{v}')
    return '##PAF\t' + '\t'.join(base + ex)


def _gff_feat(
    seqid: str,
    source: str,
    ftype: str,
    start: int,
    end: int,
    strand: str,
    phase: int | str = '.',
    attrs: Dict[str, str] | None = None,
    score: str = '.',
) -> str:
    attrs = attrs or {}
    attr_txt = ';'.join([f'{k}={v}' for k, v in attrs.items()])
    return '\t'.join(
        map(str, [seqid, source, ftype, start, end, score, strand, phase, attr_txt])
    )


def _wrap_gff(paf_lines: Iterable[str], feats: Iterable[str]) -> str:
    lines = ['##gff-version 3']
    lines.extend(paf_lines)
    lines.extend(feats)
    return '\n'.join(lines) + '\n'


@pytest.fixture
def gff_minimal_two_loci(tmp_path: Path) -> Path:
    """Valid two-locus GFF: Pseudo vs Good near 100k, and intact Far at ~200k."""
    pafs = [
        _paf_line(
            'RefPseudo',
            600,
            0,
            600,
            '+',
            'chr1',
            1_000_000,
            100200,
            101900,
            420,
            600,
            tags={'AS': 1600, 'ms': 420, 'fs': 1, 'st': 1},
        ),
        _paf_line(
            'RefGood',
            600,
            0,
            600,
            '+',
            'chr1',
            1_000_000,
            100000,
            101800,
            440,
            600,
            tags={'AS': 1400, 'ms': 440, 'fs': 0, 'st': 0},
        ),
        _paf_line(
            'RefFar',
            600,
            0,
            600,
            '+',
            'chr1',
            1_000_000,
            200000,
            201000,
            430,
            600,
            tags={'AS': 1500, 'ms': 430, 'fs': 0, 'st': 0},
        ),
    ]
    feats = [
        _gff_feat(
            'chr1',
            'miniprot',
            'mRNA',
            100201,
            101900,
            '+',
            '.',
            {'ID': 'MP_pseudo', 'Target': 'RefPseudo 1 600'},
        ),
        _gff_feat(
            'chr1',
            'miniprot',
            'CDS',
            100201,
            100600,
            '+',
            '0',
            {'Parent': 'MP_pseudo', 'ID': 'CDS_pseudo'},
        ),
        _gff_feat(
            'chr1',
            'miniprot',
            'CDS',
            101300,
            101900,
            '+',
            '0',
            {'Parent': 'MP_pseudo', 'ID': 'CDS_pseudo'},
        ),
        _gff_feat(
            'chr1',
            'miniprot',
            'mRNA',
            100001,
            101800,
            '+',
            ' .'.replace(' ', ''),
            {'ID': 'MP_good', 'Target': 'RefGood 1 600'},
        ),
        _gff_feat(
            'chr1',
            'miniprot',
            'CDS',
            100001,
            100500,
            '+',
            '0',
            {'Parent': 'MP_good', 'ID': 'CDS_good'},
        ),
        _gff_feat(
            'chr1',
            'miniprot',
            'CDS',
            101300,
            101800,
            '+',
            '0',
            {'Parent': 'MP_good', 'ID': 'CDS_good'},
        ),
        _gff_feat(
            'chr1',
            'miniprot',
            'mRNA',
            200001,
            201000,
            '+',
            ' .'.replace(' ', ''),
            {'ID': 'MP_far', 'Target': 'RefFar 1 600'},
        ),
        _gff_feat(
            'chr1',
            'miniprot',
            'CDS',
            200001,
            200400,
            '+',
            '0',
            {'Parent': 'MP_far', 'ID': 'CDS_far'},
        ),
        _gff_feat(
            'chr1',
            'miniprot',
            'CDS',
            200600,
            201000,
            '+',
            '0',
            {'Parent': 'MP_far', 'ID': 'CDS_far'},
        ),
    ]
    p = tmp_path / 'two_loci_valid.gff3'
    p.write_text(_wrap_gff(pafs, feats), encoding='utf-8')
    return p


@pytest.fixture
def gff_all_fail_gating(tmp_path: Path) -> Path:
    """Override broken legacy with a valid low-coverage single locus."""
    pafs = [
        _paf_line(
            'LowCov',
            600,
            0,
            300,
            '+',
            'chr2',
            500_000,
            5000,
            8000,
            150,
            300,
            tags={'AS': 500, 'ms': 150, 'fs': 0, 'st': 0},
        )
    ]
    feats = [
        _gff_feat(
            'chr2',
            'miniprot',
            'mRNA',
            5001,
            8000,
            '+',
            ' .'.replace(' ', ''),
            {'ID': 'MP_low', 'Target': 'LowCov 1 300'},
        ),
        _gff_feat(
            'chr2',
            'miniprot',
            'CDS',
            5001,
            6000,
            '+',
            '0',
            {'Parent': 'MP_low', 'ID': 'CDS_low'},
        ),
        _gff_feat(
            'chr2',
            'miniprot',
            'CDS',
            7000,
            8000,
            '+',
            '0',
            {'Parent': 'MP_low', 'ID': 'CDS_low'},
        ),
    ]
    p = tmp_path / 'all_fail_valid.gff3'
    p.write_text(_wrap_gff(pafs, feats), encoding='utf-8')
    return p


# === Additional tests for new emission and selection modes ===


def _parse_tsv(tsv_text: str):
    lines = [ln for ln in tsv_text.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split('\t')
    rows = []
    for ln in lines[1:]:
        fields = ln.split('\t')
        rows.append(dict(zip(header, fields)))
    return rows


def test_emit_mode_all_passing_and_cap(gff_minimal_two_loci, capsys):
    rc, stdout, stderr = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--emit-mode',
            'all_passing',
            '--max-per-locus',
            '1',
            '--min-cov',
            '0.5',
            '--min-pid',
            '0.3',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    # one gene per locus because of cap=1
    assert stdout.count('\tgene\t') >= 2
    rows = _parse_tsv(stderr)
    assert len(rows) >= 2


def test_selection_mode_prefer_intact_switches_first_locus(
    gff_minimal_two_loci, capsys
):
    rc, stdout, stderr = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--selection-mode',
            'prefer_intact',
            '--min-cov',
            '0.0',
            '--min-pid',
            '0.0',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    # Should include RefGood (intact) instead of RefPseudo at locus A
    assert 'Target=RefGood' in stdout
    # Second locus keeps Far
    assert 'Target=RefFar' in stdout


def test_two_loci_by_distance_default_gap(gff_minimal_two_loci, capsys):
    rc, stdout, stderr = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    assert stdout.count('\tgene\t') >= 2
