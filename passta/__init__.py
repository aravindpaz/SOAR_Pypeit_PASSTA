"""PASSTA: SOAR/Goodman spectroscopy helpers for PypeIt quicklook workflows."""

from passta.versioning import get_version

__all__ = ['LCOGT', '__version__', 'get_version']


def __getattr__(name: str):
    """Lazy-load :class:`~passta.lcogt.LCOGT` (heavy dependency graph)."""
    if name == 'LCOGT':
        from passta.lcogt import LCOGT

        globals()['LCOGT'] = LCOGT
        return LCOGT
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))


__version__ = get_version()
