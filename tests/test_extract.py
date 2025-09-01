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
