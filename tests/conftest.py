from __future__ import annotations

from pathlib import Path
import textwrap
from typing import Dict, List, Tuple

import pytest


@pytest.fixture
def gff_minimal_two_loci(tmp_path):
    """
    Minimal synthetic miniprot-like GFF3 with two loci on chr1.

    Locus A (100k..~102k):
      - RefPseudo (higher AS, but fs=1, st=1)
      - RefGood   (lower AS, intact: fs=0, st=0)
    Locus B (200k..~201k):
      - RefFar (intact)

    These are crafted so:
      * Default 'best' by AS picks RefPseudo for locus A.
      * '--selection-mode prefer_intact' picks RefGood for locus A.
      * Locus B always picks RefFar.
    """
    gff = textwrap.dedent("""\
        ##gff-version 3
        ##PAF\tRefPseudo\t600\t0\t600\t+\tchr1\t1000000\t100200\t101900\t900\t1800\t0\tAS:i:1600\tms:i:1600\tnp:i:500\tfs:i:1\tst:i:1\tda:i:0\tdo:i:0\tcg:Z:100M\tcs:Z::100
        chr1\tminiprot\tmRNA\t100201\t101900\t1600\t+\t.\tID=MP_pseudo;Target=RefPseudo 1 600
        chr1\tminiprot\tCDS\t100201\t100600\t.\t+\t0\tParent=MP_pseudo;ID=CDS_pseudo
        chr1\tminiprot\tCDS\t101300\t101900\t.\t+\t0\tParent=MP_pseudo;ID=CDS_pseudo
        ##PAF\tRefGood\t600\t0\t600\t+\tchr1\t1000000\t100000\t101800\t900\t1800\t0\tAS:i:1400\tms:i:1400\tnp:i:500\tfs:i:0\tst:i:0\tda:i:0\tdo:i:0\tcg:Z:100M\tcs:Z::100
        chr1\tminiprot\tmRNA\t100001\t101800\t1400\t+\t.\tID=MP_good;Target=RefGood 1 600
        chr1\tminiprot\tCDS\t100001\t100500\t.\t+\t0\tParent=MP_good;ID=CDS_good
        chr1\tminiprot\tCDS\t101300\t101800\t.\t+\t0\tParent=MP_good;ID=CDS_good
        ##PAF\tRefFar\t600\t0\t600\t+\tchr1\t1000000\t200000\t201000\t900\t1800\t0\tAS:i:1000\tms:i:1000\tnp:i:400\tfs:i:0\tst:i:0\tda:i:0\tdo:i:0\tcg:Z:100M\tcs:Z::100
        chr1\tminiprot\tmRNA\t200001\t201000\t1000\t+\t.\tID=MP_far;Target=RefFar 1 600
        chr1\tminiprot\tCDS\t200001\t200400\t.\t+\t0\tParent=MP_far;ID=CDS_far
        chr1\tminiprot\tCDS\t200600\t201000\t.\t+\t0\tParent=MP_far;ID=CDS_far
    """)
    p = tmp_path / 'two_loci.gff3'
    p.write_text(gff, encoding='utf-8')
    return p


@pytest.fixture
def gff_all_fail_gating(tmp_path):
    """
    One locus where coverage should be < min-cov gate if set high.

    qlen=600, qstart=0, qend=300 => cov=0.5, so --min-cov 0.9 will drop it.
    """
    gff = """##gff-version 3
##PAF\tLowCov\t600\t0\t300\t+\tchr2\t500000\t5000\t8000\t300\t1000\t0\tAS:i:500\tms:i:500\tnp:i:200\tfs:i:0\tst:i:0\tda:i:0\tdo:i:0\tcg:Z:100M\tcs:Z::100
chr2\tminiprot\tmRNA\t5001\t8000\t500\t+\t.\tID=MP_low;Target=LowCov 1 300
chr2\tminiprot\tCDS\t5001\t6000\t.\t+\t0\tParent=MP_low;ID=CDS_low
chr2\tminiprot\tCDS\t7000\t8000\t.\t+\t0\tParent=MP_low;ID=CDS_low
"""
    p = tmp_path / 'lowcov.gff3'
    p.write_text(gff, encoding='utf-8')
    return p


def parse_gff_attrs(attr_field: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for kv in attr_field.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
    return d


def split_gff_lines(
    gff_text: str,
) -> List[Tuple[str, str, str, int, int, str, str, str, Dict[str, str]]]:
    """
    Return parsed GFF tuples for non-comment lines:
      (seqid, source, ftype, start, end, score, strand, phase, attrs_dict)
    """
    lines = []
    for ln in gff_text.splitlines():
        if not ln or ln.startswith('#'):
            continue
        cols = ln.rstrip('\n').split('\t')
        if len(cols) < 9:
            continue
        seqid, source, ftype, start_s, end_s, score, strand, phase, attrs = cols[:9]
        lines.append(
            (
                seqid,
                source,
                ftype,
                int(start_s),
                int(end_s),
                score,
                strand,
                phase,
                parse_gff_attrs(attrs),
            )
        )
    return lines


def write_fasta(path: Path, contigs: Dict[str, str]) -> Path:
    """Write a minimal FASTA file."""
    with path.open('w', encoding='utf-8') as fh:
        for name, seq in contigs.items():
            fh.write(f'>{name}\n')
            # short sequences; no need to wrap
            fh.write(seq + '\n')
    return path


def parse_fasta(txt: str) -> Dict[str, str]:
    """Return dict of header->sequence from FASTA text."""
    out: Dict[str, str] = {}
    header = None
    chunks = []
    for ln in txt.splitlines():
        if ln.startswith('>'):
            if header is not None:
                out[header] = ''.join(chunks)
            header = ln[1:].strip()
            chunks = []
        elif ln.strip():
            chunks.append(ln.strip())
    if header is not None:
        out[header] = ''.join(chunks)
    return out


@pytest.fixture
def genome_plus(tmp_path) -> Path:
    """
    Genome with a region supporting a two-exon '+' transcript:

    Exon1 (phase=1): chrA:101-107 = GATGAAA  -> trimmed first exon by 1 => ATGAAA
    Exon2 (phase=1 but must be ignored): chrA:201-203 = TTT
    Combined CDS (transcript 5'→3'): ATGAAA + TTT = ATGAAATTT -> 'MKF'
    """
    seq = ['N'] * 1000
    # exon1 GATGAAA at 1-based [101..107]
    seq[100:107] = list('GATGAAA')
    # exon2 TTT at [201..203]
    seq[200:203] = list('TTT')
    chrA = ''.join(seq)
    return write_fasta(tmp_path / 'genome_plus.fa', {'chrA': chrA})


@pytest.fixture
def gff_plus_two_exons(tmp_path) -> Path:
    """GFF3 for the '+' strand two-exon mRNA (with later exon phase=1 to ensure it's ignored)."""
    gff = textwrap.dedent("""\
        ##gff-version 3
        chrA\tmaxiprot\tmRNA\t101\t203\t.\t+\t.\tID=tx1
        chrA\tmaxiprot\tCDS\t101\t107\t.\t+\t1\tParent=tx1;ID=cds_tx1
        chrA\tmaxiprot\tCDS\t201\t203\t.\t+\t1\tParent=tx1;ID=cds_tx1
    """)
    p = tmp_path / 'tx1_plus.gff3'
    p.write_text(gff, encoding='utf-8')
    return p


@pytest.fixture
def genome_minus(tmp_path) -> Path:
    """
    Genome with two exons on '-' strand:

    We want transcript CDS = 'ATG' + 'AAA' -> 'MK'

    For '-' strand, transcript order is descending genomic start.
    First (transcript) exon: chrB:300-303 = CATC (RC -> GATG), phase=1 => trim 'G' -> 'ATG'
    Second exon: chrB:200-202 = TTT (RC -> AAA)

    Combined: ATG + AAA -> 'MK'
    """
    seq = ['N'] * 1000
    # exonA genomic [300..303] = CATC  (RC -> GATG)
    seq[299:303] = list('CATC')
    # exonB genomic [200..202] = TTT   (RC -> AAA)
    seq[199:202] = list('TTT')
    chrB = ''.join(seq)
    return write_fasta(tmp_path / 'genome_minus.fa', {'chrB': chrB})


@pytest.fixture
def gff_minus_two_exons(tmp_path) -> Path:
    """GFF3 for '-' mRNA with two exons; phase only on the first in transcript order."""
    gff = textwrap.dedent("""\
        ##gff-version 3
        chrB\tmaxiprot\tmRNA\t200\t303\t.\t-\t.\tID=tx2
        chrB\tmaxiprot\tCDS\t300\t303\t.\t-\t1\tParent=tx2;ID=cds_tx2
        chrB\tmaxiprot\tCDS\t200\t202\t.\t-\t0\tParent=tx2;ID=cds_tx2
    """)
    p = tmp_path / 'tx2_minus.gff3'
    p.write_text(gff, encoding='utf-8')
    return p


@pytest.fixture
def genome_pseudo(tmp_path) -> Path:
    """
    Genome with a single-exon CDS that produces an internal stop under table 1:
      'ATG' 'TGA' 'GAA' -> M * E
    """
    seq = ['N'] * 500
    seq[100:109] = list('ATGTGAGAA')  # [101..109]
    chrM = ''.join(seq)
    return write_fasta(tmp_path / 'genome_pseudo.fa', {'chrM': chrM})


@pytest.fixture
def gff_pseudo(tmp_path) -> Path:
    gff = textwrap.dedent("""\
        ##gff-version 3
        chrM\tmaxiprot\tmRNA\t101\t109\t.\t+\t.\tID=txM
        chrM\tmaxiprot\tCDS\t101\t109\t.\t+\t0\tParent=txM;ID=cds_txM
    """)
    p = tmp_path / 'pseudo.gff3'
    p.write_text(gff, encoding='utf-8')
    return p


@pytest.fixture
def genome_nonacgt(tmp_path) -> Path:
    """Genome with NNN in CDS to trigger non-ACGT warning; translates to 'X' for that codon."""
    seq = ['N'] * 500
    seq[100:106] = list('ATGNNN')  # [101..106]
    chrN = ''.join(seq)
    return write_fasta(tmp_path / 'genome_nonacgt.fa', {'chrN': chrN})


@pytest.fixture
def gff_nonacgt(tmp_path) -> Path:
    gff = textwrap.dedent("""\
        ##gff-version 3
        chrN\tmaxiprot\tmRNA\t101\t106\t.\t+\t.\tID=txN
        chrN\tmaxiprot\tCDS\t101\t106\t.\t+\t0\tParent=txN;ID=cds_txN
    """)
    p = tmp_path / 'nonacgt.gff3'
    p.write_text(gff, encoding='utf-8')
    return p


@pytest.fixture
def genome_two_on_one_contig(tmp_path) -> Path:
    """Two mRNAs on the same contig to test --max-annos-per-contig."""
    seq = ['N'] * 1000
    # txA: [101..109] ATGAAAAAA (9 nt: ATG AAA AAA -> M K K)
    seq[100:109] = list('ATGAAAAAA')
    # txB: [201..209] ATGAAAAAA (another)
    seq[200:209] = list('ATGAAAAAA')
    chrZ = ''.join(seq)
    return write_fasta(tmp_path / 'genome_two.fa', {'chrZ': chrZ})


@pytest.fixture
def gff_two_mrnas_same_contig(tmp_path) -> Path:
    gff = textwrap.dedent("""\
        ##gff-version 3
        chrZ\tmaxiprot\tmRNA\t101\t109\t.\t+\t.\tID=txA
        chrZ\tmaxiprot\tCDS\t101\t109\t.\t+\t0\tParent=txA;ID=cds_txA
        chrZ\tmaxiprot\tmRNA\t201\t209\t.\t+\t.\tID=txB
        chrZ\tmaxiprot\tCDS\t201\t209\t.\t+\t0\tParent=txB;ID=cds_txB
    """)
    p = tmp_path / 'two_on_one.gff3'
    p.write_text(gff, encoding='utf-8')
    return p


@pytest.fixture
def gff_gene_minus(tmp_path) -> Path:
    """
    A gene on '-' strand with a span to be reverse-complemented in gene mode.
    Include one mRNA to be realistic, but gene extraction reads gene span only.
    """
    gff = textwrap.dedent("""\
        ##gff-version 3
        chrG\tmaxiprot\tgene\t51\t59\t.\t-\t.\tID=geneG
        chrG\tmaxiprot\tmRNA\t51\t59\t.\t-\t.\tID=txG;Parent=geneG
        chrG\tmaxiprot\tCDS\t51\t59\t.\t-\t0\tParent=txG;ID=cds_txG
    """)
    p = tmp_path / 'gene_minus.gff3'
    p.write_text(gff, encoding='utf-8')
    return p


@pytest.fixture
def genome_gene_minus(tmp_path) -> Path:
    """
    chrG[51..59] = ACGTTTAAA  (RC -> TTTAAACGT)
    """
    seq = ['N'] * 200
    seq[50:59] = list('ACGTTTAAA')
    chrG = ''.join(seq)
    return write_fasta(tmp_path / 'genome_gene.fa', {'chrG': chrG})
