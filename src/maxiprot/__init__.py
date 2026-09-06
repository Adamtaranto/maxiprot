"""maxiprot: pick the best miniprot alignment per locus from a miniprot GFF3."""

try:
    from maxiprot._version import __version__
except ImportError:  # pragma: no cover - source checkout without a build step
    from importlib.metadata import PackageNotFoundError, version

    try:
        __version__ = version('maxiprot')
    except PackageNotFoundError:
        __version__ = '0.0.0'

__all__ = ['__version__']
