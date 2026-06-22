"""Route PypeIt user-facing output into the standard :mod:`logging` tree (PASSTA log).

Older PypeIt exposes a process-global :obj:`pypeit.msgs` (:class:`~pypeit.pypmsgs.Messages`).
Current releases use :obj:`pypeit.log` (:class:`~pypeit.logger.PypeItLogger`) instead; this
module supports both.
"""

from __future__ import annotations

import contextlib
import contextvars
import io
import logging
import sys
from collections.abc import Iterator
from pathlib import Path
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

_original_pypeitlogger_init: Callable[..., None] | None = None
_pypeitlogger_init_patch_installed = False


def _legacy_pypeit_msgs_available(pypeit_mod: object) -> bool:
    msgs = getattr(pypeit_mod, 'msgs', None)
    return msgs is not None and getattr(msgs, 'reset', None) is not None


def _install_pypeit_log_init_patch() -> None:
    """Patch :meth:`pypeit.logger.PypeItLogger.init` so an active bridge steers the stream.

    PypeIt scripts call :meth:`pypeit.scripts.scriptbase.ScriptBase.init_log`, which
    reconfigures the ``pypeit`` logger and would otherwise replace any stream we set
    at bridge entry.
    """
    global _original_pypeitlogger_init, _pypeitlogger_init_patch_installed
    if _pypeitlogger_init_patch_installed:
        return
    from pypeit.logger import PypeItLogger

    _original_pypeitlogger_init = PypeItLogger.init

    def _patched_init(  # noqa: ANN001
        self,
        level: int = logging.INFO,
        stream: io.TextIOBase | None = None,
        log_file: str | Path | None = None,
        log_file_level: int | None = None,
    ) -> None:
        assert _original_pypeitlogger_init is not None
        target = _bridge_logger.get()
        if target is not None:
            stream = _MsgsLoggerWriter(target)
            log_file = None
            log_file_level = None
        return _original_pypeitlogger_init(
            self,
            level=level,
            stream=stream,
            log_file=log_file,
            log_file_level=log_file_level,
        )

    PypeItLogger.init = _patched_init  # type: ignore[method-assign]
    _pypeitlogger_init_patch_installed = True


def _pypeit_console_level() -> int:
    """Match legacy bridge (verbosity=1): INFO on the PypeIt console stream."""
    from pypeit import log as pypeit_log

    return pypeit_log.convert_verbosity_to_logging_level(1)


def _apply_modern_pypeit_log_sink(log: logging.Logger) -> None:
    from pypeit import log as pypeit_log

    pypeit_log.init(
        level=_pypeit_console_level(),
        stream=_MsgsLoggerWriter(log),
        log_file=None,
    )


def _restore_modern_pypeit_log_default() -> None:
    from pypeit import log as pypeit_log

    pypeit_log.init(level=logging.DEBUG, stream=sys.stderr, log_file=None)


def _ensure_pypeit_user() -> None:
    """Define ``pypeit.pypeit_user`` if missing (needed before ``pypeit.pypmsgs`` imports).

    Recent PypeIt releases dropped ``pypeit_user`` from :file:`pypeit/__init__.py`, but
    ``pypmsgs.py`` still does ``from pypeit import pypeit_user``, which breaks
    ``from pypeit.pypmsgs import Messages`` on some pip/conda installs. Match the
    historical behaviour (see upstream ``pypeit`` package init on older tags).
    """
    import getpass

    import pypeit

    if not hasattr(pypeit, 'pypeit_user'):
        pypeit.pypeit_user = getpass.getuser()


def _install_msgs_reset_patch() -> None:
    global _original_messages_reset, _msgs_reset_patch_installed
    _ensure_pypeit_user()
    from pypeit.pypmsgs import Messages

    if _msgs_reset_patch_installed:
        return
    _original_messages_reset = Messages.reset
    Messages.reset = _patched_messages_reset
    _msgs_reset_patch_installed = True


def _restore_default_msgs() -> None:
    import pypeit

    if not _legacy_pypeit_msgs_available(pypeit):
        _restore_modern_pypeit_log_default()
        return
    pypeit.msgs.reset(log=None, verbosity=1, log_to_stderr=True, colors=False)


@contextlib.contextmanager
def pypeit_logging_bridge(log: logging.Logger | None = None) -> Iterator[None]:
    """While active, PypeIt script output is routed to ``log``.

    Legacy PypeIt uses ``pypeit.msgs``; current releases use ``pypeit.log``
    (:class:`~pypeit.logger.PypeItLogger`). Under this context, console output
    is redirected to the given logger (and separate PypeIt log files are
    suppressed while the bridge is active).

    Parameters
    ----------
    log : `logging.Logger` or None
        Target logger. Default: ``logging.getLogger('passta')``.

    Notes
    -----
    The legacy ``msgs`` object and the modern ``pypeit`` named logger are
    process-global. Nested ``with pypeit_logging_bridge()`` blocks are
    supported. When the outermost block ends, sinks are reset to PypeIt
    defaults.

    Parallel threads each running PypeIt with different bridge targets can race
    on the shared logger; keep ``--jobs 1`` if you see interleaved output, or
    share one logger for all workers.
    """
    _install_pypeit_log_init_patch()
    lg = log or logging.getLogger('passta')
    token = _bridge_logger.set(lg)
    import pypeit

    legacy = _legacy_pypeit_msgs_available(pypeit)
    if legacy:
        _install_msgs_reset_patch()
        pypeit.msgs.reset(log=None, verbosity=1, log_to_stderr=True, colors=False)
    else:
        # Modern PypeIt: ``PypeItLogger.init`` is patched so ``init_log`` keeps this sink.
        _apply_modern_pypeit_log_sink(lg)
    try:
        yield
    finally:
        _bridge_logger.reset(token)
        parent = _bridge_logger.get()
        if parent is None:
            _restore_default_msgs()
        elif legacy:
            pypeit.msgs.reset(log=None, verbosity=1, log_to_stderr=True, colors=False)
        else:
            _apply_modern_pypeit_log_sink(parent)
