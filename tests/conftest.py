from __future__ import annotations

from pathlib import Path
import textwrap
from typing import Dict, Iterable, List, Tuple

import pytest

# -----------------------------
# Builders for miniprot-shaped GFF3 fixtures
# -----------------------------


def make_cs(identical: int, substituted: int = 0) -> str:
    """Build a miniprot-style cs:Z string with the given composition."""
    cs = f':{identical}' if identical else ''
    cs += '*aa' * substituted
    return cs


def paf_line(
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
    """Build a ##PAF header line as miniprot emits it."""
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


def gff_feat(
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
    """Build a single 9-column GFF3 feature line."""
    attrs = attrs or {}
    attr_txt = ';'.join([f'{k}={v}' for k, v in attrs.items()])
    return '\t'.join(
        map(str, [seqid, source, ftype, start, end, score, strand, phase, attr_txt])
    )


def wrap_gff(paf_lines: Iterable[str], feats: Iterable[str]) -> str:
    """Assemble PAF headers and feature lines into a GFF3 document."""
    lines = ['##gff-version 3']
    lines.extend(paf_lines)
    lines.extend(feats)
    return '\n'.join(lines) + '\n'


# -----------------------------
# GFF3 validity checker (shared by emission tests)
# -----------------------------


def assert_valid_gff3(text: str) -> None:
    """Assert structural GFF3 validity of maxiprot filter output.

    Checks: gff-version header; 1-based start<=end; unique gene/mRNA IDs;
    every Parent resolves to an earlier feature; all CDS parts of one mRNA
    share exactly one ID; features are coordinate-sorted by (seqid, gene
    start) at the gene level.
    """
    lines = text.splitlines()
    assert lines and lines[0] == '##gff-version 3'

    seen_ids: Dict[str, str] = {}  # ID -> ftype
    cds_ids_by_parent: Dict[str, set] = {}
    gene_order: List[Tuple[str, int]] = []
    for ln in lines[1:]:
        if not ln or ln.startswith('#'):
            continue
        cols = ln.split('\t')
        assert len(cols) == 9, f'Bad column count: {ln!r}'
        seqid, _source, ftype, start_s, end_s, _score, strand, _phase, attrs_s = cols
        start, end = int(start_s), int(end_s)
        assert 1 <= start <= end, f'Bad coordinates: {ln!r}'
        assert strand in {'+', '-', '.', '?'}
        attrs = dict(x.split('=', 1) for x in attrs_s.split(';') if '=' in x)
        fid = attrs.get('ID')
        parent = attrs.get('Parent')
        if parent is not None:
            assert parent in seen_ids, f'Unresolved Parent {parent!r} in: {ln!r}'
        if ftype in {'gene', 'mRNA'}:
            assert fid, f'{ftype} without ID: {ln!r}'
            assert fid not in seen_ids, f'Duplicate {ftype} ID {fid!r}'
            seen_ids[fid] = ftype
            if ftype == 'gene':
                gene_order.append((seqid, start))
        elif ftype == 'CDS':
            assert fid, f'CDS without ID: {ln!r}'
            assert parent, f'CDS without Parent: {ln!r}'
            cds_ids_by_parent.setdefault(parent, set()).add(fid)
            seen_ids.setdefault(fid, 'CDS')

    for parent, ids in cds_ids_by_parent.items():
        assert len(ids) == 1, (
            f'CDS parts of {parent!r} carry {len(ids)} distinct IDs: {ids}'
        )
    assert gene_order == sorted(gene_order), 'genes are not coordinate-sorted'


# -----------------------------
# Filter fixtures (miniprot-shaped: CDS lines carry NO ID attribute)
# -----------------------------


@pytest.fixture
def gff_minimal_two_loci(tmp_path: Path) -> Path:
    """Two loci on chr1(+).

    Locus A: RefPseudo (higher ms/AS but fs=1,st=1) vs RefGood (intact,
    higher pid). Locus B: RefFar (intact), 100 kb away.

    Default ms_cov_pos scoring picks RefPseudo at locus A;
    --selection-mode prefer_intact picks RefGood; pid_cov picks RefGood.
    """
    pafs = [
        paf_line(
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
            tags={
                'AS': 1600,
                'ms': 1500,
                'np': 520,
                'fs': 1,
                'st': 1,
                'cs': make_cs(420, 180),
            },
        ),
        paf_line(
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
            tags={
                'AS': 1400,
                'ms': 1400,
                'np': 500,
                'fs': 0,
                'st': 0,
                'cs': make_cs(440, 160),
            },
        ),
        paf_line(
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
            tags={
                'AS': 1500,
                'ms': 1450,
                'np': 510,
                'fs': 0,
                'st': 0,
                'cs': make_cs(430, 170),
            },
        ),
    ]
    feats = [
        gff_feat(
            'chr1',
            'miniprot',
            'mRNA',
            100201,
            101900,
            '+',
            '.',
            {'ID': 'MP_pseudo', 'Target': 'RefPseudo 1 600'},
        ),
        gff_feat(
            'chr1',
            'miniprot',
            'CDS',
            100201,
            100600,
            '+',
            '0',
            {'Parent': 'MP_pseudo', 'Target': 'RefPseudo 1 133'},
        ),
        gff_feat(
            'chr1',
            'miniprot',
            'CDS',
            101300,
            101900,
            '+',
            '0',
            {'Parent': 'MP_pseudo', 'Target': 'RefPseudo 134 600'},
        ),
        gff_feat(
            'chr1',
            'miniprot',
            'mRNA',
            100001,
            101800,
            '+',
            '.',
            {'ID': 'MP_good', 'Target': 'RefGood 1 600'},
        ),
        gff_feat(
            'chr1',
            'miniprot',
            'CDS',
            100001,
            100500,
            '+',
            '0',
            {'Parent': 'MP_good', 'Target': 'RefGood 1 167'},
        ),
        gff_feat(
            'chr1',
            'miniprot',
            'CDS',
            101300,
            101800,
            '+',
            '0',
            {'Parent': 'MP_good', 'Target': 'RefGood 168 600'},
        ),
        gff_feat(
            'chr1',
            'miniprot',
            'mRNA',
            200001,
            201000,
            '+',
            '.',
            {'ID': 'MP_far', 'Target': 'RefFar 1 600'},
        ),
        gff_feat(
            'chr1',
            'miniprot',
            'CDS',
            200001,
            200400,
            '+',
            '0',
            {'Parent': 'MP_far', 'Target': 'RefFar 1 133'},
        ),
        gff_feat(
            'chr1',
            'miniprot',
            'CDS',
            200600,
            201000,
            '+',
            '0',
            {'Parent': 'MP_far', 'Target': 'RefFar 134 600'},
        ),
    ]
    p = tmp_path / 'two_loci.gff3'
    p.write_text(wrap_gff(pafs, feats), encoding='utf-8')
    return p


@pytest.fixture
def gff_all_fail_gating(tmp_path: Path) -> Path:
    """One locus whose only candidate has cov=0.5 (fails --min-cov 0.9)."""
    pafs = [
        paf_line(
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
            tags={
                'AS': 500,
                'ms': 450,
                'np': 200,
                'fs': 0,
                'st': 0,
                'cs': make_cs(150, 150),
            },
        )
    ]
    feats = [
        gff_feat(
            'chr2',
            'miniprot',
            'mRNA',
            5001,
            8000,
            '+',
            '.',
            {'ID': 'MP_low', 'Target': 'LowCov 1 300'},
        ),
        gff_feat(
            'chr2',
            'miniprot',
            'CDS',
            5001,
            6000,
            '+',
            '0',
            {'Parent': 'MP_low', 'Target': 'LowCov 1 150'},
        ),
        gff_feat(
            'chr2',
            'miniprot',
            'CDS',
            7000,
            8000,
            '+',
            '0',
            {'Parent': 'MP_low', 'Target': 'LowCov 151 300'},
        ),
    ]
    p = tmp_path / 'all_fail.gff3'
    p.write_text(wrap_gff(pafs, feats), encoding='utf-8')
    return p


# -----------------------------
# utilities for FASTA generation
# -----------------------------


def parse_gff_attrs(attr_txt: str) -> Dict[str, str]:
    """Parse a GFF attribute column into a dict."""
    out: Dict[str, str] = {}
    for kv in attr_txt.split(';'):
        if not kv or '=' not in kv:
            continue
        k, v = kv.split('=', 1)
        out[k] = v
    return out


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
    chunks: List[str] = []
    for ln in txt.splitlines():
        if ln.startswith('>'):
            if header is not None:
                out[header] = ''.join(chunks)
            header = ln[1:].strip()
            chunks = []
        else:
            chunks.append(ln.strip())
    if header is not None:
        out[header] = ''.join(chunks)
    return out


# -----------------------------
# Fixtures for extract submodule
# -----------------------------
@pytest.fixture
def genome_plus(tmp_path) -> Path:
    """
    Genome with three exons on '+' strand:

    We want transcript CDS = 'ATG' + 'AAA' + 'TTT' -> 'MKF'

    Exon1: chrA:101-103 = ATG (phase 0)
    Exon2: chrA:201-203 = AAA (phase 0)
    Exon3: chrA:301-303 = TTT (phase 0)

    Combined: ATG + AAA + TTT => 'MKF'
    """
    seq = ['N'] * 1000
    seq[100:103] = list('ATG')  # [101..103]
    seq[200:203] = list('AAA')  # [201..203]
    seq[300:303] = list('TTT')  # [301..303]
    chrA = ''.join(seq)
    return write_fasta(tmp_path / 'genome_plus.fa', {'chrA': chrA})


@pytest.fixture
def gff_plus_two_exons(tmp_path) -> Path:
    """
    Three CDS parts on '+' strand for a single transcript.
    """
    gff = textwrap.dedent("""\
        ##gff-version 3
        chrA\tmaxiprot\tmRNA\t101\t303\t.\t+\t.\tID=tx1
        chrA\tmaxiprot\tCDS\t101\t103\t.\t+\t0\tParent=tx1;ID=cds_tx1
        chrA\tmaxiprot\tCDS\t201\t203\t.\t+\t0\tParent=tx1;ID=cds_tx1
        chrA\tmaxiprot\tCDS\t301\t303\t.\t+\t0\tParent=tx1;ID=cds_tx1
    """)
    p = tmp_path / 'tx1_plus.gff3'
    p.write_text(gff, encoding='utf-8')
    return p


@pytest.fixture
def gff_plus_cds_before_mrna(tmp_path) -> Path:
    """
    Same transcript as gff_plus_two_exons but with CDS lines BEFORE the mRNA
    line (legal in GFF3; regression test for shell-mRNA merging).
    """
    gff = textwrap.dedent("""\
        ##gff-version 3
        chrA\tmaxiprot\tCDS\t101\t103\t.\t+\t0\tParent=tx1;ID=cds_tx1
        chrA\tmaxiprot\tCDS\t201\t203\t.\t+\t0\tParent=tx1;ID=cds_tx1
        chrA\tmaxiprot\tCDS\t301\t303\t.\t+\t0\tParent=tx1;ID=cds_tx1
        chrA\tmaxiprot\tmRNA\t101\t303\t.\t+\t.\tID=tx1
    """)
    p = tmp_path / 'tx1_cds_first.gff3'
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
    """
    Two CDS parts on '-' strand with phase on first exon (transcript order).
    """
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
    # txB: [201..209] ATGAAAAAA (another on same contig)
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
    Gene on '-' strand, two CDS parts.
    """
    gff = textwrap.dedent("""\
        ##gff-version 3
        chrG\tmaxiprot\tgene\t51\t59\t.\t-\t.\tID=geneM
        chrG\tmaxiprot\tmRNA\t51\t59\t.\t-\t.\tID=txGM;Parent=geneM
        chrG\tmaxiprot\tCDS\t56\t59\t.\t-\t0\tParent=txGM;ID=cds_txGM
        chrG\tmaxiprot\tCDS\t51\t53\t.\t-\t0\tParent=txGM;ID=cds_txGM
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
