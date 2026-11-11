"""Mirror a step's log output to a file, in addition to the console.

``cli.main`` configures console logging once; steps that accept ``--log PATH``
attach a file handler for the duration of the step with :func:`log_to_file`.
"""

from __future__ import annotations

import logging
from contextlib import contextmanager
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

LOG_FORMAT = "%(asctime)s %(name)s %(levelname)s: %(message)s"


@contextmanager
def log_to_file(path: Path) -> Iterator[None]:
    """Write everything logged inside the block to ``path`` as well as the console.

    The parent directory is created. The handler is removed on exit, so
    repeated in-process calls do not stack handlers.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(path)
    handler.setFormatter(logging.Formatter(LOG_FORMAT))
    root = logging.getLogger()
    root.addHandler(handler)
    try:
        yield
    finally:
        root.removeHandler(handler)
        handler.close()
