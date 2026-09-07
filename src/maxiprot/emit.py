"""GFF3 emission for ``maxiprot filter`` winners.

Builds a spec-compliant, coordinate-sorted ``gene -> mRNA -> CDS`` hierarchy
for each emitted candidate. All CDS parts of one transcript share a single
``ID`` (the GFF3 discontinuous-feature convention); miniprot's own GFF3 CDS
lines carry no ``ID`` attribute, so CDS IDs are always synthesized here.
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

from maxiprot.gffio import GffFeature, format_gff_line
from maxiprot.ioutils import open_output

logger = logging.getLogger(__name__)

SOURCE = 'maxiprot'

#: mRNA attributes carried over from the input when reusing a feature.
_MRNA_KEEP = ('Rank', 'Identity', 'Positive', 'Frameshift', 'StopCodon', 'Target')
#: CDS attributes carried over from the input when reusing features.
_CDS_KEEP = ('Rank', 'Identity', 'Frameshift', 'StopCodon', 'Target')


def _unique_id(base: str, used: set[str]) -> str:
    """Return ``base`` or ``base:N`` so the result is unique within ``used``."""
    candidate = base
    n = 1
    while candidate in used:
        n += 1
        candidate = f'{base}:{n}'
    used.add(candidate)
    return candidate


def _reused_block(
    winner: pd.Series,
    mrna: GffFeature,
    cds_parts: list[GffFeature],
    gene_id: str,
    mrna_out_id: str,
) -> tuple[str, int, list[str]]:
    """Build the gene block for a winner matched to an input mRNA feature."""
    lines = [
        format_gff_line(
            mrna.seqid,
            SOURCE,
            'gene',
            mrna.start,
            mrna.end,
            '.',
            mrna.strand,
            '.',
            {
                'ID': gene_id,
                'Name': winner['qname'],
                'locus': winner['locus'],
            },
        )
    ]
    m_attrs: dict[str, object] = {'ID': mrna_out_id, 'Parent': gene_id}
    for keep in _MRNA_KEEP:
        if keep in mrna.attrs:
            m_attrs[keep] = mrna.attrs[keep]
    lines.append(
        format_gff_line(
            mrna.seqid,
            mrna.source,
            'mRNA',
            mrna.start,
            mrna.end,
            mrna.score or '.',
            mrna.strand,
            '.',
            m_attrs,
        )
    )
    cds_id = f'{mrna_out_id}.cds'
    for cds in sorted(cds_parts, key=lambda f: f.start):
        c_attrs: dict[str, object] = {'ID': cds_id, 'Parent': mrna_out_id}
        for keep in _CDS_KEEP:
            if keep in cds.attrs:
                c_attrs[keep] = cds.attrs[keep]
        lines.append(
            format_gff_line(
                cds.seqid,
                cds.source,
                'CDS',
                cds.start,
                cds.end,
                cds.score or '.',
                cds.strand,
                cds.phase,
                c_attrs,
            )
        )
    return mrna.seqid, mrna.start, lines


def _synthesized_block(
    winner: pd.Series,
    gene_id: str,
    mrna_out_id: str,
) -> tuple[str, int, list[str]]:
    """Build a minimal gene block from PAF coordinates (no input mRNA match)."""
    logger.warning(
        'No matching mRNA feature for %s on %s; synthesizing a single-span '
        'CDS from PAF coordinates (introns are not represented).',
        winner['qname'],
        winner['tname'],
    )
    seqid = str(winner['tname'])
    strand = str(winner['strand'])
    # PAF ts is 0-based, GFF3 is 1-based inclusive.
    start = int(winner['ts']) + 1
    end = int(winner['te'])
    cds_id = f'{mrna_out_id}.cds'
    target = f'{winner["qname"]} {int(winner["qs"]) + 1} {int(winner["qe"])}'
    lines = [
        format_gff_line(
            seqid,
            SOURCE,
            'gene',
            start,
            end,
            '.',
            strand,
            '.',
            {'ID': gene_id, 'Name': winner['qname'], 'locus': winner['locus']},
        ),
        format_gff_line(
            seqid,
            SOURCE,
            'mRNA',
            start,
            end,
            '.',
            strand,
            '.',
            {'ID': mrna_out_id, 'Parent': gene_id, 'Target': target},
        ),
        format_gff_line(
            seqid,
            SOURCE,
            'CDS',
            start,
            end,
            '.',
            strand,
            '0',
            {'ID': cds_id, 'Parent': mrna_out_id},
        ),
    ]
    return seqid, start, lines


def write_best_gff3(
    out_path: Optional[str],
    winners: pd.DataFrame,
    mrna_by_id: dict[str, GffFeature],
    cds_by_parent: dict[str, list[GffFeature]],
    id_prefix: str = 'PBM',
    tlens: Optional[dict[str, int]] = None,
) -> None:
    """
    Write emitted candidates as coordinate-sorted, spec-compliant GFF3.

    Parameters
    ----------
    out_path : str or None
        Output destination (path, or '-'/None for stdout).
    winners : pandas.DataFrame
        Emitted candidates; requires columns ``locus``, ``qname``, ``tname``,
        ``strand``, ``ts``, ``te``, ``qs``, ``qe`` and optionally ``mrna_id``.
    mrna_by_id : dict of str to GffFeature
        Input mRNA features keyed by their ``ID`` attribute.
    cds_by_parent : dict of str to list of GffFeature
        Input CDS features grouped by their ``Parent`` attribute.
    id_prefix : str
        Prefix used when synthesizing IDs.
    tlens : dict of str to int, optional
        Contig lengths (from PAF ``tlen``) for ``##sequence-region`` lines.
    """
    if out_path and out_path != '-':
        logger.info('Writing filtered records to GFF3: %s', out_path)
    else:
        logger.info('Writing filtered records to stdout')

    used_ids: set[str] = set()
    blocks: list[tuple[str, int, list[str]]] = []

    for _, winner in winners.iterrows():
        base = f'{id_prefix}:{winner["locus"]}:{winner["qname"]}'
        gene_id = _unique_id(f'{base}:gene', used_ids)

        mrna_id = winner.get('mrna_id')
        if pd.notna(mrna_id) and mrna_id in mrna_by_id:
            mrna = mrna_by_id[str(mrna_id)]
            mrna_out_id = _unique_id(str(mrna_id), used_ids)
            blocks.append(
                _reused_block(
                    winner,
                    mrna,
                    cds_by_parent.get(str(mrna_id), []),
                    gene_id,
                    mrna_out_id,
                )
            )
        else:
            mrna_out_id = _unique_id(f'{base}:mrna', used_ids)
            blocks.append(_synthesized_block(winner, gene_id, mrna_out_id))

    blocks.sort(key=lambda b: (b[0], b[1]))

    out, close_me = open_output(out_path)
    try:
        out.write('##gff-version 3\n')
        if tlens:
            for seqid in sorted({b[0] for b in blocks}):
                if seqid in tlens:
                    out.write(f'##sequence-region {seqid} 1 {tlens[seqid]}\n')
        for _, _, lines in blocks:
            for line in lines:
                out.write(line + '\n')
    finally:
        if close_me:
            out.close()
