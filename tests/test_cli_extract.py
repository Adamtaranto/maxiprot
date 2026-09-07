# tests/test_cli_extract.py
from __future__ import annotations

import pytest

from maxiprot import __version__, cli as cli_mod

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


def test_cli_dispatches_to_filter(gff_minimal_two_loci, capsys):
    rc = cli_mod.main(['filter', str(gff_minimal_two_loci), '--log-level', 'ERROR'])
    out = capsys.readouterr().out
    assert rc == 0
    assert 'Target=RefPseudo' in out


def test_cli_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        cli_mod.main(['--version'])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_cli_subcommand_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        cli_mod.main(['filter', '--version'])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_cli_requires_subcommand(capsys):
    with pytest.raises(SystemExit) as exc:
        cli_mod.main([])
    assert exc.value.code == 2
