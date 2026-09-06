"""Golden-file integration tests on genuine miniprot output.

The fixtures under tests/data/ were generated once with:
    miniprot --gff mini_genome.fa mini_refs.faa > mini_miniprot.gff3
and the goldens with:
    maxiprot filter tests/data/mini_miniprot.gff3 --min-cov 0 --min-pid 0 \
        --out-gff3 expected_best.gff3 --out-tsv expected_best.tsv
"""

from __future__ import annotations

from pathlib import Path
import shutil

from maxiprot import extract as extract_mod, filter as filter_mod

from .conftest import assert_valid_gff3, parse_fasta

DATA = Path(__file__).parent / 'data'


def test_filter_matches_golden_outputs(tmp_path, capsys):
    out_gff = tmp_path / 'best.gff3'
    out_tsv = tmp_path / 'best.tsv'
    rc = filter_mod.main(
        [
            str(DATA / 'mini_miniprot.gff3'),
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
        ]
    )
    capsys.readouterr()
    assert rc == 0
    assert out_gff.read_text() == (DATA / 'expected_best.gff3').read_text()
    assert out_tsv.read_text() == (DATA / 'expected_best.tsv').read_text()


def test_filter_output_on_real_miniprot_is_valid_gff3(tmp_path, capsys):
    out_gff = tmp_path / 'best.gff3'
    rc = filter_mod.main(
        [
            str(DATA / 'mini_miniprot.gff3'),
            '--min-cov',
            '0',
            '--min-pid',
            '0',
            '--out-gff3',
            str(out_gff),
            '--log-level',
            'ERROR',
        ]
    )
    capsys.readouterr()
    assert rc == 0
    text = out_gff.read_text()
    assert_valid_gff3(text)
    assert 'ID=nan' not in text


def test_tsv_metrics_are_fractions_on_real_data(tmp_path, capsys):
    out_gff = tmp_path / 'best.gff3'
    out_tsv = tmp_path / 'best.tsv'
    rc = filter_mod.main(
        [
            str(DATA / 'mini_miniprot.gff3'),
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
        ]
    )
    capsys.readouterr()
    assert rc == 0
    lines = out_tsv.read_text().splitlines()
    header = lines[0].split('\t')
    for ln in lines[1:]:
        row = dict(zip(header, ln.split('\t')))
        for col in ('cov_aa', 'pid_aa', 'positives'):
            assert 0.0 <= float(row[col]) <= 1.0, (col, row[col])


def test_extract_proteins_from_filtered_real_data(tmp_path, capsys):
    # Copy the genome so pyfaidx writes its .fai into tmp_path, not tests/data.
    genome = tmp_path / 'mini_genome.fa'
    shutil.copy(DATA / 'mini_genome.fa', genome)
    rc = extract_mod.main(
        [
            str(DATA / 'expected_best.gff3'),
            '-g',
            str(genome),
            '--extract',
            'protein',
            '--log-level',
            'ERROR',
        ]
    )
    out = capsys.readouterr().out
    assert rc == 0
    recs = parse_fasta(out)
    assert len(recs) == 2
    for seq in recs.values():
        assert len(seq) > 100
        assert set(seq) <= set('ACDEFGHIKLMNPQRSTVWYX*')
    # Both candidates are degraded (frameshifts/stops) -> flagged pseudogene
    assert sum('pseudogene=1' in h for h in recs) == 2
