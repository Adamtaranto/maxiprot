from __future__ import annotations

import io
from typing import List

from maxiprot import filter as filter_mod  # the module under test


def run_filter_and_capture(argv: List[str], capsys):
    """Helper to run filter.main(argv) and capture stdout/stderr."""
    rc = filter_mod.main(argv)
    out = capsys.readouterr()
    return rc, out.out, out.err


def test_default_io_and_best_selection(gff_minimal_two_loci, capsys):
    """
    Default behavior:
      - GFF3 -> stdout
      - TSV  -> stderr
      - Selection mode 'best' by AS (we force gates off so AS dominates)
        => choose RefPseudo for first locus (higher AS, but pseudogene),
           choose RefFar for second locus.
    """
    rc, stdout, stderr = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--score-mode',
            'AS',
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
      - First locus should pick RefGood (intact) even though RefPseudo has higher AS.
      - Second locus remains RefFar (intact).
    """
    rc, stdout, stderr = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--score-mode',
            'AS',
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
    assert 'Target=RefPseudo' not in stdout
    assert 'Target=RefFar' in stdout


def test_locus_clustering_counts_genes(gff_minimal_two_loci, capsys):
    """
    Two loci on same contig/strand far apart should yield two gene features.
    """
    rc, stdout, _ = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--score-mode',
            'AS',
            '--min-cov',
            '0',
            '--min-pid',
            '0',
            '--log-level',
            'ERROR',
            # default --locus-pad 5000 is plenty to split 100k vs 200k
        ],
        capsys,
    )
    assert rc == 0
    n_genes = sum(1 for ln in stdout.splitlines() if '\tgene\t' in ln)
    assert n_genes == 2


def test_outfile_options_quiet_streams(gff_minimal_two_loci, tmp_path, capsys):
    """
    When --out-gff3 and --out-tsv are provided, nothing should be emitted to stdout/stderr
    (with --log-level ERROR). Files should be created with expected content.
    """
    out_gff = tmp_path / 'best.gff3'
    out_tsv = tmp_path / 'best.tsv'
    rc, stdout, stderr = run_filter_and_capture(
        [
            str(gff_minimal_two_loci),
            '--score-mode',
            'AS',
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
    assert 'RefFar' in gff_txt
    assert 'Ref' in tsv_txt  # at least some TSV rows with qnames


def test_stdin_path_reads_from_stdin(gff_minimal_two_loci, capsys, monkeypatch):
    """
    Using '-' as the GFF path should read from stdin.
    """
    stdin = gff_minimal_two_loci.read_text(encoding='utf-8')
    monkeypatch.setattr('sys.stdin', io.StringIO(stdin))
    rc, stdout, stderr = run_filter_and_capture(
        [
            '-',  # read from stdin
            '--score-mode',
            'AS',
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
    assert 'Target=RefFar' in stdout
    assert 'RefFar' in stderr


def test_strict_drops_locus_when_all_fail_gates(gff_all_fail_gating, capsys):
    """
    With --strict and tight gates, a locus with no passing candidates should be dropped.
    """
    rc, stdout, stderr = run_filter_and_capture(
        [
            str(gff_all_fail_gating),
            '--score-mode',
            'AS',
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
    print(stdout)
    assert rc == 0
    # No gene should be emitted
    assert '\tgene\t' not in stdout
    # TSV may be empty or contain a note; avoid being too strict here:
    # Accept either empty or a header-only.
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
            'AS',
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
    # map: mRNA ID -> set of CDS IDs seen
    mrna_id = None
    cds_ids = set()
    for ln in stdout.splitlines():
        if not ln or ln.startswith('#'):
            continue
        cols = ln.split('\t')
        if len(cols) < 9:
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
