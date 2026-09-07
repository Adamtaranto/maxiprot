from __future__ import annotations

import io
from typing import List

from maxiprot import filter as filter_mod  # the module under test

from .conftest import assert_valid_gff3, make_cs, paf_line, wrap_gff


def run_filter_and_capture(argv: List[str], capsys):
    """Run filter.main(argv) and capture stdout/stderr."""
    rc = filter_mod.main(argv)
    out = capsys.readouterr()
    return rc, out.out, out.err


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


def test_default_io_and_best_selection(gff_minimal_two_loci, capsys):
    """
    Default behavior:
      * default ms_cov_pos scoring picks RefPseudo at locus A, RefFar at locus B
      * GFF -> stdout; no TSV without --out-tsv
    """
    rc, stdout, stderr = run_filter_and_capture(
        [str(gff_minimal_two_loci), '--log-level', 'ERROR'],
        capsys,
    )
    assert rc == 0
    assert 'Target=RefPseudo' in stdout
    assert 'Target=RefFar' in stdout
    assert 'Target=RefGood' not in stdout
    # TSV is not written by default (stderr only carries logs, silenced here)
    assert stderr == ''


def test_prefer_intact_over_best(gff_minimal_two_loci, capsys):
    """--selection-mode prefer_intact picks intact RefGood at locus A."""
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--selection-mode',
            'prefer_intact',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    assert 'Target=RefGood' in stdout
    assert 'Target=RefFar' in stdout


def test_score_mode_pid_cov_prefers_higher_identity(gff_minimal_two_loci, capsys):
    """pid_cov ranks by identity*coverage: RefGood beats RefPseudo."""
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--score-mode',
            'pid_cov',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    assert 'Target=RefGood' in stdout
    assert 'Target=RefPseudo' not in stdout


def test_file_outputs_redirecting(gff_minimal_two_loci, tmp_path, capsys):
    """With --out-gff3/--out-tsv set, outputs go to files; stdout/stderr clean."""
    out_gff = tmp_path / 'out.gff3'
    out_tsv = tmp_path / 'out.tsv'
    rc, stdout, stderr = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
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
    assert stdout == ''
    assert stderr == ''
    gff_txt = out_gff.read_text(encoding='utf-8')
    tsv_txt = out_tsv.read_text(encoding='utf-8')
    assert 'Target=RefPseudo' in gff_txt
    assert 'Target=RefFar' in gff_txt
    assert 'RefPseudo' in tsv_txt and 'RefFar' in tsv_txt


def test_tsv_stdout_requires_gff3_file(gff_minimal_two_loci, capsys):
    """--out-tsv - with GFF3 on stdout is a usage error."""
    rc, _, _ = run_filter_and_capture(
        [str(gff_minimal_two_loci), '--out-tsv', '-', '--log-level', 'ERROR'],
        capsys,
    )
    assert rc == 2


def test_tsv_schema_and_metric_bounds(gff_minimal_two_loci, tmp_path, capsys):
    """TSV has the documented columns; pid/cov/positives are fractions."""
    out_gff = tmp_path / 'out.gff3'
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--out-gff3',
            str(out_gff),
            '--out-tsv',
            '-',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    rows = _parse_tsv(stdout)
    assert len(rows) == 2
    expected_cols = set(filter_mod.TSV_COLUMNS)
    assert set(rows[0]) == expected_cols
    for row in rows:
        for col in ('cov_aa', 'pid_aa', 'positives'):
            assert 0.0 <= float(row[col]) <= 1.0, (col, row[col])
        assert row['status'] in {'intact', 'pseudogene'}
        assert row['passes'] in {'True', 'False'}
        # locus label format tname:strand:index
        assert row['locus'].startswith(f'{row["tname"]}:{row["strand"]}:')
        # 1-based coordinates
        assert int(row['tstart']) >= 1
        assert int(row['tend']) >= int(row['tstart'])


def test_stdin_path_reads_from_stdin(gff_minimal_two_loci, capsys, monkeypatch):
    """Using '-' as the GFF path reads from stdin."""
    data = gff_minimal_two_loci.read_text(encoding='utf-8')
    monkeypatch.setattr('sys.stdin', io.StringIO(data))
    rc, stdout, _ = run_filter_and_capture(
        ['-', '--log-level', 'ERROR'],
        capsys,
    )
    assert rc == 0
    assert 'Target=RefPseudo' in stdout
    assert 'Target=RefFar' in stdout


def test_strict_drops_locus_when_all_fail_gates(gff_all_fail_gating, capsys):
    """With --strict and tight gates, a failing locus is dropped entirely."""
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_all_fail_gating),
            '--min-cov',
            '0.9',  # candidate has 0.5
            '--strict',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    assert '\tgene\t' not in stdout


def test_nonstrict_fallback_flags_failing_winner(gff_all_fail_gating, tmp_path, capsys):
    """Without --strict, the best-scoring candidate is emitted with passes=False."""
    out_gff = tmp_path / 'out.gff3'
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_all_fail_gating),
            '--min-cov',
            '0.9',
            '--out-gff3',
            str(out_gff),
            '--out-tsv',
            '-',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    rows = _parse_tsv(stdout)
    assert len(rows) == 1
    assert rows[0]['qname'] == 'LowCov'
    assert rows[0]['passes'] == 'False'
    assert '\tgene\t' in out_gff.read_text(encoding='utf-8')


def test_gates_exclude_failing_candidates(gff_minimal_two_loci, tmp_path, capsys):
    """A gate that only RefGood passes must exclude RefPseudo from selection."""
    out_gff = tmp_path / 'out.gff3'
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--min-pid',
            '0.72',  # RefGood 0.7333 passes; RefPseudo 0.70 fails
            '--out-gff3',
            str(out_gff),
            '--out-tsv',
            '-',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    gff_txt = out_gff.read_text(encoding='utf-8')
    assert 'Target=RefGood' in gff_txt
    assert 'Target=RefPseudo' not in gff_txt


def test_output_is_valid_gff3(gff_minimal_two_loci, capsys):
    """Structural GFF3 validity: shared CDS IDs, resolvable Parents, sorting."""
    rc, stdout, _ = run_filter_and_capture(
        [str(gff_minimal_two_loci), '--log-level', 'ERROR'],
        capsys,
    )
    assert rc == 0
    assert_valid_gff3(stdout)
    assert 'nan' not in stdout


def test_cds_lines_share_single_id_per_transcript(gff_minimal_two_loci, capsys):
    """Multi-exon CDS share one synthesized ID even though input CDS has none."""
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--selection-mode',
            'prefer_intact',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    assert_valid_gff3(stdout)
    cds_ids = set()
    for ln in stdout.splitlines():
        cols = ln.split('\t')
        if len(cols) == 9 and cols[2] == 'CDS' and 'Parent=MP_good' in cols[8]:
            attrs = dict(x.split('=', 1) for x in cols[8].split(';') if '=' in x)
            cds_ids.add(attrs['ID'])
    assert len(cds_ids) == 1


def test_emit_mode_all_passing_and_cap(gff_minimal_two_loci, tmp_path, capsys):
    out_gff = tmp_path / 'out.gff3'
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--emit-mode',
            'all_passing',
            '--max-per-locus',
            '1',
            '--out-gff3',
            str(out_gff),
            '--out-tsv',
            '-',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    gff_txt = out_gff.read_text(encoding='utf-8')
    # one gene per locus because of cap=1
    assert gff_txt.count('\tgene\t') == 2
    rows = _parse_tsv(stdout)
    assert len(rows) == 2


def test_emit_mode_all_passing_uncapped(gff_minimal_two_loci, capsys):
    """all_passing without a cap emits both passers at locus A."""
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--emit-mode',
            'all_passing',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    assert stdout.count('\tgene\t') == 3
    assert_valid_gff3(stdout)


def test_two_loci_by_distance_default_gap(gff_minimal_two_loci, capsys):
    rc, stdout, _ = run_filter_and_capture(
        [str(gff_minimal_two_loci), '--log-level', 'ERROR'],
        capsys,
    )
    assert rc == 0
    assert stdout.count('\tgene\t') == 2


def test_locus_pad_merges_distant_loci(gff_minimal_two_loci, capsys):
    """A huge --locus-pad merges both loci into one gene call."""
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--locus-pad',
            '200000',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    assert stdout.count('\tgene\t') == 1


def test_selection_mode_longest(gff_minimal_two_loci, tmp_path, capsys):
    """'longest' ranks by aligned AA length (all 600 here) with score tiebreak."""
    out_gff = tmp_path / 'out.gff3'
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--selection-mode',
            'longest',
            '--out-gff3',
            str(out_gff),
            '--out-tsv',
            '-',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    rows = _parse_tsv(stdout)
    assert len(rows) == 2
    assert all(int(r['cds_aa_len']) == 600 for r in rows)


def test_invalid_gate_fraction_is_usage_error(gff_minimal_two_loci, capsys):
    rc, _, _ = run_filter_and_capture(
        [str(gff_minimal_two_loci), '--min-pid', '30', '--log-level', 'ERROR'],
        capsys,
    )
    assert rc == 2


def test_no_paf_records_returns_2(tmp_path, capsys):
    p = tmp_path / 'nopaf.gff3'
    p.write_text(
        '##gff-version 3\nchr1\tminiprot\tmRNA\t1\t10\t.\t+\t.\tID=x\n',
        encoding='utf-8',
    )
    rc, _, _ = run_filter_and_capture([str(p), '--log-level', 'ERROR'], capsys)
    assert rc == 2


def test_missing_tags_do_not_crash(tmp_path, capsys):
    """PAF lines without ms/AS/np/fs/st/cs must be handled gracefully."""
    pafs = [
        paf_line('Bare', 100, 0, 100, '+', 'chrB', 10_000, 1000, 1300, 80, 100),
    ]
    p = tmp_path / 'bare.gff3'
    p.write_text(wrap_gff(pafs, []), encoding='utf-8')
    out = run_filter_and_capture(
        [str(p), '--min-cov', '0', '--min-pid', '0', '--log-level', 'ERROR'],
        capsys,
    )
    rc, stdout, _ = out
    assert rc == 0
    assert '\tgene\t' in stdout


def test_synthesized_hierarchy_uses_1based_coords(tmp_path, capsys):
    """Winners without a matching mRNA get PAF coords converted to 1-based."""
    pafs = [
        paf_line(
            'Solo',
            100,
            0,
            100,
            '+',
            'chrS',
            10_000,
            999,
            1300,
            90,
            100,
            tags={'AS': 300, 'ms': 290, 'cs': make_cs(90, 10)},
        ),
    ]
    p = tmp_path / 'solo.gff3'
    p.write_text(wrap_gff(pafs, []), encoding='utf-8')
    rc, stdout, _ = run_filter_and_capture(
        [str(p), '--min-cov', '0', '--min-pid', '0', '--log-level', 'ERROR'],
        capsys,
    )
    assert rc == 0
    gene_lines = [
        ln
        for ln in stdout.splitlines()
        if len(ln.split('\t')) == 9 and ln.split('\t')[2] == 'gene'
    ]
    assert len(gene_lines) == 1
    cols = gene_lines[0].split('\t')
    assert cols[3] == '1000'  # ts=999 (0-based) -> 1000 (1-based)
    assert cols[4] == '1300'
    assert_valid_gff3(stdout)


def test_all_passing_id_collisions_are_deduplicated(tmp_path, capsys):
    """Two alignments of the same query in one locus must not share gene IDs."""
    pafs = [
        paf_line(
            'Dup',
            100,
            0,
            100,
            '+',
            'chrD',
            100_000,
            1000,
            1300,
            90,
            100,
            tags={'AS': 300, 'ms': 290, 'cs': make_cs(90, 10)},
        ),
        paf_line(
            'Dup',
            100,
            0,
            100,
            '+',
            'chrD',
            100_000,
            2000,
            2300,
            85,
            100,
            tags={'AS': 280, 'ms': 270, 'cs': make_cs(85, 15)},
        ),
    ]
    p = tmp_path / 'dup.gff3'
    p.write_text(wrap_gff(pafs, []), encoding='utf-8')
    rc, stdout, _ = run_filter_and_capture(
        [
            str(p),
            '--emit-mode',
            'all_passing',
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
    assert stdout.count('\tgene\t') == 2
    assert_valid_gff3(stdout)  # asserts unique gene/mRNA IDs
