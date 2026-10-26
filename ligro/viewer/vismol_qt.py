"""Qt-native VisMol integration for LiGRO.

Requires the patched ``graphics_engine`` package. The GTK host is not imported:
the original VismolSession and VismolGLCore run inside a PyQt5 QOpenGLWidget.
"""
from __future__ import annotations

from typing import Optional
from PyQt5 import QtWidgets

try:
    from vismol.core.vismol_session import VismolSession
    from vismol.gui.vismol_qtwidget import VismolQtWidget
except ImportError as exc:
    raise ImportError(
        "LiGRO's VisMol viewer requires the patched graphics_engine package. "
        "From the LiGRO project root run: pip install -e vendor/graphics_engine"
    ) from exc


GraphicsEngineAdapter = VismolQtWidget


def create_vismol_viewer(parent: Optional[QtWidgets.QWidget] = None):
    """Create and return ``(session, widget)`` ready for a PyQt layout."""
    session = VismolSession(toolkit="Qt5")
    widget = session.vm_widget
    if parent is not None:
        widget.setParent(parent)
    return session, widget


__all__ = ["VismolQtWidget", "GraphicsEngineAdapter", "create_vismol_viewer"]
