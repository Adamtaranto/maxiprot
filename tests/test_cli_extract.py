# tests/test_cli_extract.py
from __future__ import annotations

from maxiprot import cli as cli_mod

from .conftest import parse_fasta


def test_cli_dispatches_to_extract(genome_plus, gff_plus_two_exons, capsys):
    rc = cli_mod.main(
        [
            'extract',
            str(gff_plus_two_exons),
            '-g',
            str(genome_plus),
            '--extract',
            'protein',
            '--log-level',
            'ERROR',
        ]
    )
    out = capsys.readouterr().out
    recs = parse_fasta(out)
    assert rc == 0
    assert next(iter(recs.values())) == 'MKF'
