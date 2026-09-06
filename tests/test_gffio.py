from __future__ import annotations

from maxiprot.gffio import (
    GffFeature,
    PafRecord,
    format_gff_line,
    gff3_escape,
    gff3_unescape,
    iter_gff3,
    parse_attributes,
    parse_paf_line,
)


def test_escape_specials_roundtrip():
    raw = 'a;b=c,d&e%f\tg\nh'
    escaped = gff3_escape(raw)
    for ch in (';', '=', ',', '&', '\t', '\n'):
        assert ch not in escaped
    assert gff3_unescape(escaped) == raw


def test_escape_none_is_empty():
    assert gff3_escape(None) == ''


def test_parse_attributes_urldecodes():
    attrs = parse_attributes('ID=x%3B1;Target=Ref 1 600;Note=a%2Cb')
    assert attrs['ID'] == 'x;1'
    assert attrs['Target'] == 'Ref 1 600'
    assert attrs['Note'] == 'a,b'


def test_format_gff_line_escapes_dict_attrs():
    line = format_gff_line('c', 's', 'gene', 1, 10, '.', '+', '.', {'ID': 'a;b'})
    assert line.endswith('ID=a%3Bb')
    assert len(line.split('\t')) == 9


def test_parse_paf_line_core_and_tags():
    payload = (
        '\t'.join(
            ['Q', '100', '0', '100', '+', 'c', '1000', '10', '310', '80', '100', '60']
        )
        + '\tAS:i:300\tms:i:290\tcs:Z::80\tbadtag\txx:i:notint'
    )
    rec = parse_paf_line(payload)
    assert isinstance(rec, PafRecord)
    assert rec.qname == 'Q'
    assert rec.ts == 10
    assert rec.tags['AS'] == 300
    assert rec.tags['cs'] == ':80'
    assert 'badtag' not in rec.tags
    assert 'xx' not in rec.tags


def test_parse_paf_line_malformed_returns_none():
    assert parse_paf_line('too\tfew\tfields') is None
    assert (
        parse_paf_line(
            '\t'.join(['Q', 'NaN', '0', '1', '+', 'c', '9', '0', '1', '1', '1', '0'])
        )
        is None
    )


def test_iter_gff3_yields_features_and_paf():
    lines = [
        '##gff-version 3\n',
        '##PAF\tQ\t100\t0\t100\t+\tc\t1000\t10\t310\t80\t100\t60\tAS:i:300\n',
        '# a comment\n',
        'c\tsrc\tmRNA\t11\t310\t.\t+\t.\tID=m1;Target=Q 1 100\n',
        'c\tsrc\tCDS\t11\t100\t.\t+\t0\tParent=m1\n',
        'short\tline\n',
    ]
    items = list(iter_gff3(lines))
    assert len(items) == 3
    assert isinstance(items[0], PafRecord)
    assert isinstance(items[1], GffFeature)
    assert items[1].type == 'mRNA'
    assert items[1].attrs['ID'] == 'm1'
    assert items[2].type == 'CDS'
    assert items[2].attrs['Parent'] == 'm1'
