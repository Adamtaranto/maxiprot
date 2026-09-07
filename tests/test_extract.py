# tests/test_extract.py
from __future__ import annotations

import io
from pathlib import Path
from typing import List

from maxiprot import extract as extract_mod

from .conftest import parse_fasta


def run_extract_and_capture(argv: List[str], capsys):
    rc = extract_mod.main(argv)
    out = capsys.readouterr()
    return rc, out.out, out.err


def test_protein_plus_strand_phase_only_first_exon(
    genome_plus, gff_plus_two_exons, capsys
):
    rc, stdout, _ = run_extract_and_capture(
        [
            str(gff_plus_two_exons),
            '-g',
            str(genome_plus),
            '--extract',
            'protein',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    recs = parse_fasta(stdout)
    assert len(recs) == 1
    assert next(iter(recs.values())) == 'MKF'


def test_protein_minus_strand_phase_on_first_transcript_exon(
    genome_minus, gff_minus_two_exons, capsys
):
    rc, stdout, _ = run_extract_and_capture(
        [
            str(gff_minus_two_exons),
            '-g',
            str(genome_minus),
            '--extract',
            'protein',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    recs = parse_fasta(stdout)
    assert next(iter(recs.values())) == 'MK'


def test_cds_mode_returns_concatenated_nt(genome_plus, gff_plus_two_exons, capsys):
    rc, stdout, _ = run_extract_and_capture(
        [
            str(gff_plus_two_exons),
            '-g',
            str(genome_plus),
            '--extract',
            'cds',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    recs = parse_fasta(stdout)
    assert next(iter(recs.values())) == 'ATGAAATTT'


def test_gene_mode_reverse_complements_span(genome_gene_minus, gff_gene_minus, capsys):
    rc, stdout, _ = run_extract_and_capture(
        [
            str(gff_gene_minus),
            '-g',
            str(genome_gene_minus),
            '--extract',
            'gene',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    recs = parse_fasta(stdout)
    assert next(iter(recs.values())) == 'TTTAAACGT'


def test_exclude_pseudogene_under_table1(genome_pseudo, gff_pseudo, capsys):
    rc, stdout, _ = run_extract_and_capture(
        [
            str(gff_pseudo),
            '-g',
            str(genome_pseudo),
            '--extract',
            'protein',
            '--exclude-pseudogenes',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    assert parse_fasta(stdout) == {}  # filtered out


def test_transl_table_2_rescues_TGA(genome_pseudo, gff_pseudo, capsys):
    rc, stdout, _ = run_extract_and_capture(
        [
            str(gff_pseudo),
            '-g',
            str(genome_pseudo),
            '--extract',
            'protein',
            '--exclude-pseudogenes',
            '--transl-table',
            '2',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    recs = parse_fasta(stdout)
    assert next(iter(recs.values())) == 'MWE'


def test_warn_non_acgt_and_translate_to_X(genome_nonacgt, gff_nonacgt, caplog, capsys):
    caplog.set_level('WARNING')
    rc, stdout, _ = run_extract_and_capture(
        [
            str(gff_nonacgt),
            '-g',
            str(genome_nonacgt),
            '--extract',
            'protein',
            '--log-level',
            'WARNING',
        ],
        capsys,
    )
    assert rc == 0
    recs = parse_fasta(stdout)
    assert next(iter(recs.values())) == 'MX'
    assert any('non-ACGT' in rec.message for rec in caplog.records)


def test_max_annos_per_contig_gates_out_all(
    genome_two_on_one_contig, gff_two_mrnas_same_contig, caplog, capsys
):
    caplog.set_level('WARNING')
    rc, stdout, _ = run_extract_and_capture(
        [
            str(gff_two_mrnas_same_contig),
            '-g',
            str(genome_two_on_one_contig),
            '--extract',
            'protein',
            '--max-annos-per-contig',
            '1',
            '--log-level',
            'WARNING',
        ],
        capsys,
    )
    assert rc == 0
    assert parse_fasta(stdout) == {}
    assert any('Skipping all sequences on' in rec.message for rec in caplog.records)


def test_outfile_path_writes_to_file(genome_plus, gff_plus_two_exons, tmp_path, capsys):
    out = tmp_path / 'out.faa'
    rc, stdout, _ = run_extract_and_capture(
        [
            str(gff_plus_two_exons),
            '-g',
            str(genome_plus),
            '--extract',
            'protein',
            '--out-faa',
            str(out),
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0 and stdout == ''
    recs = parse_fasta(out.read_text(encoding='utf-8'))
    assert next(iter(recs.values())) == 'MKF'


def test_stdin_gff_reads_from_stdin(
    genome_plus, gff_plus_two_exons, monkeypatch, capsys
):
    monkeypatch.setattr(
        'sys.stdin', io.StringIO(gff_plus_two_exons.read_text(encoding='utf-8'))
    )
    rc, stdout, _ = run_extract_and_capture(
        ['-', '-g', str(genome_plus), '--extract', 'cds', '--log-level', 'ERROR'],
        capsys,
    )
    assert rc == 0
    recs = parse_fasta(stdout)
    assert next(iter(recs.values())) == 'ATGAAATTT'


def test_fai_created_if_missing(genome_plus, gff_plus_two_exons, capsys):
    fai = Path(str(genome_plus) + '.fai')
    if fai.exists():
        fai.unlink()
    rc, stdout, _ = run_extract_and_capture(
        [
            str(gff_plus_two_exons),
            '-g',
            str(genome_plus),
            '--extract',
            'protein',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0 and fai.exists()


def test_cds_before_mrna_order_keeps_all_parts(
    genome_plus, gff_plus_cds_before_mrna, capsys
):
    """CDS lines preceding their mRNA line must not lose collected parts."""
    rc, stdout, _ = run_extract_and_capture(
        [
            str(gff_plus_cds_before_mrna),
            '-g',
            str(genome_plus),
            '--extract',
            'cds',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 0
    recs = parse_fasta(stdout)
    assert next(iter(recs.values())) == 'ATGAAATTT'


def test_reverse_complement_preserves_iupac():
    from maxiprot.extract import reverse_complement

    assert reverse_complement('ACGTRYN') == 'NRYACGT'


def test_transl_table_4_TGA_is_tryptophan():
    from maxiprot.extract import translate_nt

    assert translate_nt('ATGTGA', 4) == 'MW'
    assert translate_nt('ATGTGA', 1) == 'M*'


def test_transl_table_11_is_supported():
    from maxiprot.extract import translate_nt

    assert translate_nt('ATGAAA', 11) == 'MK'


def test_unsupported_transl_table_rejected_by_argparse(
    genome_plus, gff_plus_two_exons, capsys
):
    import pytest

    with pytest.raises(SystemExit) as exc:
        extract_mod.main(
            [
                str(gff_plus_two_exons),
                '-g',
                str(genome_plus),
                '--transl-table',
                '99',
            ]
        )
    assert exc.value.code == 2


def test_missing_contig_returns_1(genome_plus, tmp_path, capsys, caplog):
    caplog.set_level('ERROR')
    gff = tmp_path / 'badcontig.gff3'
    gff.write_text(
        '##gff-version 3\n'
        'chrMISSING\tmaxiprot\tmRNA\t1\t9\t.\t+\t.\tID=txX\n'
        'chrMISSING\tmaxiprot\tCDS\t1\t9\t.\t+\t0\tParent=txX\n',
        encoding='utf-8',
    )
    rc, stdout, _ = run_extract_and_capture(
        [
            str(gff),
            '-g',
            str(genome_plus),
            '--extract',
            'protein',
            '--log-level',
            'ERROR',
        ],
        capsys,
    )
    assert rc == 1
    assert parse_fasta(stdout) == {}
    assert any('not found in genome FASTA' in rec.message for rec in caplog.records)


def test_phase_mismatch_on_downstream_part_warns(genome_plus, tmp_path, capsys, caplog):
    """A downstream CDS phase disagreeing with cumulative length warns."""
    caplog.set_level('WARNING')
    gff = tmp_path / 'frameshift.gff3'
    # Parts of length 3 and 3: expected phase of second part is 0, not 2.
    gff.write_text(
        '##gff-version 3\n'
        'chrA\tmaxiprot\tmRNA\t101\t203\t.\t+\t.\tID=txF\n'
        'chrA\tmaxiprot\tCDS\t101\t103\t.\t+\t0\tParent=txF\n'
        'chrA\tmaxiprot\tCDS\t201\t203\t.\t+\t2\tParent=txF\n',
        encoding='utf-8',
    )
    rc, _, _ = run_extract_and_capture(
        [
            str(gff),
            '-g',
            str(genome_plus),
            '--extract',
            'cds',
            '--log-level',
            'WARNING',
        ],
        capsys,
    )
    assert rc == 0
    assert any('frameshifted' in rec.message for rec in caplog.records)


def test_fasta_wrap_width():
    import io as _io

    from maxiprot.extract import write_fasta

    buf = _io.StringIO()
    write_fasta([('long', 'A' * 130)], buf)
    lines = buf.getvalue().splitlines()
    assert lines[0] == '>long'
    assert [len(x) for x in lines[1:]] == [60, 60, 10]
