"""Route PypeIt's ``msgs`` output into the standard :mod:`logging` tree (PASSTA log)."""

from __future__ import annotations

import contextlib
import contextvars
import io
import logging
from collections.abc import Iterator
from typing import Callable

# Active bridge target; ``None`` means PypeIt uses its default ``msgs`` behaviour.
_bridge_logger: contextvars.ContextVar[logging.Logger | None] = contextvars.ContextVar(
    'passta_pypeit_msgs_bridge_logger',
    default=None,
)

_original_messages_reset: Callable[..., None] | None = None


class _MsgsLoggerWriter(io.TextIOBase):
    """File-like sink so ``pypmsgs.Messages`` writes into a :class:`logging.Logger`."""

    encoding = 'utf-8'

    def __init__(self, log: logging.Logger) -> None:
        super().__init__()
        self._lg = log
        self._buf = ''
        self._closed = False

    @property
    def closed(self) -> bool:  # noqa: D401
        return self._closed

    def writable(self) -> bool:
        return True

    def close(self) -> None:
        if self._closed:
            return
        if self._buf.strip():
            self._emit_line(self._buf)
        self._buf = ''
        self._closed = True

    def write(self, s: str) -> int:
        if self._closed:
            raise ValueError('I/O operation on closed file.')
        if not s:
            return 0
        self._buf += s
        while '\n' in self._buf:
            line, self._buf = self._buf.split('\n', 1)
            self._emit_line(line)
        return len(s)

    def flush(self) -> None:
        if self._closed:
            return
        if self._buf.strip():
            self._emit_line(self._buf)
            self._buf = ''

    def _emit_line(self, line: str) -> None:
        line = line.rstrip('\r\n')
        if not line.strip():
            return
        up = line.upper()
        if 'THIS LOG WAS GENERATED' in up or up.strip('-') == '':
            self._lg.debug('%s', line)
        elif '[ERROR]' in up:
            self._lg.error('%s', line)
        elif '[WARNING]' in up or '[WARN]' in up:
            self._lg.warning('%s', line)
        elif '[BUG]' in up or '[TEST]' in up or '[WORK' in up or '[PROGRESS]' in up:
            self._lg.debug('%s', line)
        elif 'YOU ARE USING ' in up and 'VERSION' in up:
            self._lg.debug('%s', line)
        else:
            self._lg.info('%s', line)


def _patched_messages_reset(
    self,
    log: str | io.IOBase | None = None,
    verbosity: int | None = None,
    colors: bool = True,
    log_to_stderr: bool | None = None,
) -> None:
    assert _original_messages_reset is not None
    target = _bridge_logger.get()
    if target is not None:
        writer = _MsgsLoggerWriter(target)
        return _original_messages_reset(
            self,
            log=writer,
            verbosity=verbosity,
            colors=False,
            log_to_stderr=False,
        )
    return _original_messages_reset(
        self,
        log=log,
        verbosity=verbosity,
        colors=colors,
        log_to_stderr=log_to_stderr,
    )


_msgs_reset_patch_installed = False


def _install_msgs_reset_patch() -> None:
    global _original_messages_reset, _msgs_reset_patch_installed
    from pypeit.pypmsgs import Messages

    if _msgs_reset_patch_installed:
        return
    _original_messages_reset = Messages.reset
    Messages.reset = _patched_messages_reset
    _msgs_reset_patch_installed = True


def _restore_default_msgs() -> None:
    import pypeit

    pypeit.msgs.reset(log=None, verbosity=1, log_to_stderr=True, colors=False)


@contextlib.contextmanager
def pypeit_logging_bridge(log: logging.Logger | None = None) -> Iterator[None]:
    """While active, all ``pypeit.msgs.reset(...)`` calls route to ``log``.

    PypeIt scripts and ``PypeIt`` call ``msgs.reset`` with a log file path or
    ``None``; under this context, those calls are redirected to the given
    logger instead (no separate PypeIt log files from PASSTA).

    Parameters
    ----------
    log : `logging.Logger` or None
        Target logger. Default: ``logging.getLogger('passta')``.

    Notes
    -----
    ``pypeit.msgs`` is process-global. Nested ``with pypeit_logging_bridge()``
    blocks are supported. When the outermost block ends and no bridge remains
    active, ``msgs`` is reset to PypeIt defaults.

    Parallel threads each running PypeIt with different bridge targets can race
    on the shared ``msgs`` object; keep ``--jobs 1`` if you see interleaved
    output, or share one logger for all workers.
    """
    _install_msgs_reset_patch()
    lg = log or logging.getLogger('passta')
    token = _bridge_logger.set(lg)
    import pypeit

    # Apply the bridge immediately; otherwise ``msgs`` keeps prior sinks until
    # the next ``reset`` call inside PypeIt.
    pypeit.msgs.reset(log=None, verbosity=1, log_to_stderr=True, colors=False)
    try:
        yield
    finally:
        _bridge_logger.reset(token)
        if _bridge_logger.get() is None:
            _restore_default_msgs()
        else:
            # Outer bridge still active: re-apply its sink after this inner exit.
            pypeit.msgs.reset(log=None, verbosity=1, log_to_stderr=True, colors=False)
