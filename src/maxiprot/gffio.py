"""GFF3 and miniprot ``##PAF`` parsing/formatting shared by maxiprot subcommands.

This module is the single source of truth for reading and writing GFF3 rows
and for parsing the ``##PAF`` header lines miniprot emits with ``--gff``.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
import logging
import urllib.parse

logger = logging.getLogger(__name__)

GFF3_COLS = (
    'seqid',
    'source',
    'type',
    'start',
    'end',
    'score',
    'strand',
    'phase',
    'attributes',
)

#: Characters that must be percent-encoded in GFF3 attribute values.
#: '%' must be first so already-inserted escapes are not double-escaped.
_ESCAPES = (
    ('%', '%25'),
    (';', '%3B'),
    ('=', '%3D'),
    ('&', '%26'),
    (',', '%2C'),
    ('\t', '%09'),
    ('\n', '%0A'),
    ('\r', '%0D'),
)


@dataclass
class GffFeature:
    """A single GFF3 feature row with parsed attributes."""

    seqid: str
    source: str
    type: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    score: str
    strand: str
    phase: str
    attrs: dict[str, str] = field(default_factory=dict)


@dataclass
class PafRecord:
    """One miniprot ``##PAF`` header line (protein-to-genome alignment)."""

    qname: str
    qlen: int  # query (protein) length in AA
    qs: int  # 0-based query start (AA)
    qe: int  # 0-based-exclusive query end (AA)
    strand: str
    tname: str
    tlen: int
    ts: int  # 0-based target start (nt)
    te: int  # 0-based-exclusive target end (nt)
    nmatch: int  # number of matching residues/bases
    alen: int  # alignment block length
    mapq: int
    tags: dict[str, int | float | str] = field(default_factory=dict)


def gff3_escape(x: object) -> str:
    """
    Escape a value for placement in a GFF3 attribute field.

    Parameters
    ----------
    x : object
        Raw attribute value; ``None`` becomes ``''``.

    Returns
    -------
    str
        Percent-encoded string per the GFF3 specification.
    """
    if x is None:
        return ''
    s = str(x)
    for raw, enc in _ESCAPES:
        s = s.replace(raw, enc)
    return s


def gff3_unescape(x: str) -> str:
    """
    Decode percent-encoded characters in a GFF3 attribute value.

    Parameters
    ----------
    x : str
        Possibly percent-encoded attribute value.

    Returns
    -------
    str
        Decoded string.
    """
    if '%' not in x:
        return x
    return urllib.parse.unquote(x)


def parse_attributes(attr_text: str) -> dict[str, str]:
    """
    Parse a GFF3 column-9 attribute string into a dict.

    Values are URL-decoded; keys are kept verbatim.

    Parameters
    ----------
    attr_text : str
        The raw attributes column.

    Returns
    -------
    dict of str to str
        Attribute key/value mapping.
    """
    attrs: dict[str, str] = {}
    for kv in attr_text.split(';'):
        kv = kv.strip()
        if not kv or '=' not in kv:
            continue
        k, v = kv.split('=', 1)
        attrs[k] = gff3_unescape(v)
    return attrs


def format_attributes(attrs: dict[str, object]) -> str:
    """
    Render an attribute dict as a GFF3 column-9 string with escaping.

    Parameters
    ----------
    attrs : dict
        Attribute key/value mapping (values are escaped; keys are not).

    Returns
    -------
    str
        ``key=value;key=value`` string.
    """
    return ';'.join(f'{k}={gff3_escape(v)}' for k, v in attrs.items())


def format_gff_line(
    seqid: str,
    source: str,
    ftype: str,
    start: int,
    end: int,
    score: str,
    strand: str,
    phase: str,
    attrs: dict[str, object] | str,
) -> str:
    """
    Format one GFF3 feature row (no trailing newline).

    Parameters
    ----------
    seqid, source, ftype, score, strand, phase : str
        GFF3 columns.
    start, end : int
        1-based inclusive coordinates.
    attrs : dict or str
        Attributes; a dict is rendered with escaping, a str is used verbatim.

    Returns
    -------
    str
        Tab-separated GFF3 row.
    """
    attr_text = format_attributes(attrs) if isinstance(attrs, dict) else attrs
    return '\t'.join(
        [
            seqid,
            source,
            ftype,
            str(start),
            str(end),
            score,
            strand,
            phase,
            attr_text,
        ]
    )


def parse_paf_line(payload: str) -> PafRecord | None:
    """
    Parse the payload of a ``##PAF`` line (text after the ``##PAF\\t`` prefix).

    Malformed core fields cause the record to be skipped (returns ``None``
    with a warning); malformed optional tags are skipped individually.

    Parameters
    ----------
    payload : str
        Tab-separated PAF fields.

    Returns
    -------
    PafRecord or None
        Parsed record, or ``None`` if the core fields are malformed.
    """
    fields = payload.rstrip('\n').split('\t')
    if len(fields) < 12:
        logger.warning('Skipping malformed ##PAF line (fewer than 12 fields)')
        return None
    try:
        rec = PafRecord(
            qname=fields[0],
            qlen=int(fields[1]),
            qs=int(fields[2]),
            qe=int(fields[3]),
            strand=fields[4],
            tname=fields[5],
            tlen=int(fields[6]),
            ts=int(fields[7]),
            te=int(fields[8]),
            nmatch=int(fields[9]),
            alen=int(fields[10]),
            mapq=int(fields[11]),
        )
    except ValueError:
        logger.warning('Skipping malformed ##PAF line (non-integer core field)')
        return None

    for tag in fields[12:]:
        parts = tag.split(':', 2)
        if len(parts) != 3:
            continue
        key, typ, val = parts
        if typ == 'i':
            try:
                rec.tags[key] = int(val)
            except ValueError:
                continue
        elif typ == 'f':
            try:
                rec.tags[key] = float(val)
            except ValueError:
                continue
        else:
            rec.tags[key] = val
    return rec


def iter_gff3(lines: Iterable[str]) -> Iterator[GffFeature | PafRecord]:
    """
    Stream GFF3 lines, yielding features and miniprot ``##PAF`` records.

    Comment/directive lines other than ``##PAF`` are skipped, as are rows
    without 9 tab-separated columns or with non-integer coordinates.

    Parameters
    ----------
    lines : iterable of str
        GFF3 lines.

    Yields
    ------
    GffFeature or PafRecord
        Parsed rows in input order.
    """
    for line in lines:
        line = line.rstrip('\n')
        if not line:
            continue
        if line.startswith('#'):
            if line.startswith('##PAF\t'):
                rec = parse_paf_line(line[len('##PAF\t') :])
                if rec is not None:
                    yield rec
            continue
        cols = line.split('\t')
        if len(cols) < 9:
            continue
        try:
            start = int(cols[3])
            end = int(cols[4])
        except ValueError:
            logger.warning('Skipping GFF3 row with non-integer coordinates')
            continue
        yield GffFeature(
            seqid=cols[0],
            source=cols[1],
            type=cols[2],
            start=start,
            end=end,
            score=cols[5],
            strand=cols[6],
            phase=cols[7],
            attrs=parse_attributes(cols[8]),
        )
