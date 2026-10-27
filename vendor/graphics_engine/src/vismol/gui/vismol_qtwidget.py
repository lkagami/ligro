"""PyQt5 QOpenGLWidget host for VisMol's existing graphics engine.

The renderer and molecular model are reused from graphics_engine; only the GTK
host is replaced. All GL calls happen while Qt has made this widget's context
current.
"""
from __future__ import annotations

import logging
import traceback

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QOpenGLWidget

from vismol.libgl.vismol_glcore import VismolGLCore

logger = logging.getLogger(__name__)


def _patch_gl_line_width():
    from OpenGL import GL

    if getattr(GL.glLineWidth, "_vismol_patched", False):
        return

    original_gl_line_width = GL.glLineWidth

    def clamped_gl_line_width(width):
        try:
            low, high = GL.glGetFloatv(GL.GL_ALIASED_LINE_WIDTH_RANGE)
            flags = GL.glGetIntegerv(GL.GL_CONTEXT_FLAGS)
            if int(flags) & GL.GL_CONTEXT_FLAG_FORWARD_COMPATIBLE_BIT:
                high = min(high, 1.0)
            width = min(max(float(width), float(low)), float(high))
        except Exception:
            pass
        return original_gl_line_width(width)

    clamped_gl_line_width._vismol_patched = True
    GL.glLineWidth = clamped_gl_line_width


class VismolQtWidget(QOpenGLWidget):
    """Qt-native host for ``VismolGLCore`` and a ``VismolSession``."""

    viewer_error = QtCore.pyqtSignal(str)

    def __init__(self, vismol_session=None, width=640, height=420, parent=None):
        super().__init__(parent)
        fmt = QtGui.QSurfaceFormat()
        fmt.setRenderableType(QtGui.QSurfaceFormat.OpenGL)
        fmt.setVersion(3, 3)
        fmt.setProfile(QtGui.QSurfaceFormat.CoreProfile)
        fmt.setDepthBufferSize(24)
        fmt.setStencilBufferSize(8)
        self.setFormat(fmt)
        self.setMinimumSize(320, 240)
        self.setFocusPolicy(Qt.StrongFocus)
        self.setMouseTracking(True)
        _patch_gl_line_width()
        self.vm_session = vismol_session
        self.vm_glcore = VismolGLCore(self, vismol_session, width, height)
        self.selection_box_frame = None
        self.glmenu_bg = self.glmenu_sele = self.glmenu_obj = self.glmenu_pick = None
        self.filechooser = None
        self._gl_initialized = False
        self._render_failed = False
        self._last_mouse_pos = QtCore.QPoint()

    # Compatibility methods expected by VismolGLCore.
    def queue_draw(self):
        self.update()

    def process_updates(self):
        self.update()
        application = QtWidgets.QApplication.instance()
        if application is not None:
            application.processEvents()

    def set_has_depth_buffer(self, *_args, **_kwargs):
        # Qt's QSurfaceFormat requests a depth buffer before context creation.
        return None

    def set_has_alpha(self, *_args, **_kwargs):
        return None

    def get_window(self):
        # Used only by the legacy GL context-menu path. Qt integration provides
        # its own context menu in the future; returning None is safe here.
        return None

    def show_gl_menu(self, *_args, **_kwargs):
        return None

    def initializeGL(self):
        try:
            self.vm_glcore.initialize()
            self._gl_initialized = True
        except Exception as exc:
            self.viewer_error.emit(f"VisMol OpenGL initialization failed: {exc}")
            raise

    def resizeGL(self, width, height):
        if self._gl_initialized:
            self.vm_glcore.resize_window(max(1, width), max(1, height))

    def paintGL(self):
        if self._gl_initialized and not self._render_failed:
            try:
                self.vm_glcore.render()
            except Exception as exc:
                self._render_failed = True
                error_details = traceback.format_exc()
                logger.error("VisMol rendering failed:\n%s", error_details)
                self.viewer_error.emit(
                    "VisMol rendering failed: {}\n\n{}".format(exc, error_details)
                )

    def load_structure(self, path, **kwargs):
        if self.vm_session is None:
            raise RuntimeError("VismolQtWidget has no VismolSession")
        return self.vm_session.load_molecule(str(path), **kwargs)

    def mousePressEvent(self, event):
        self.setFocus()
        self._last_mouse_pos = event.pos()
        button = self._button_number(event.button())
        if button:
            self.vm_glcore.mouse_pressed(button, event.x(), event.y())
        self.update()
        event.accept()

    def mouseReleaseEvent(self, event):
        button = self._button_number(event.button())
        if button:
            self.vm_glcore.mouse_released(button, event.x(), event.y())
        self.update()
        event.accept()

    def mouseMoveEvent(self, event):
        # Do not pass Qt's button bitmask as ``button_state``: the legacy
        # core interprets that optional argument as a Gdk.ModifierType mask
        # and imports GTK. Mouse modes are already set by mouse_pressed().
        self.vm_glcore.mouse_motion(event.x(), event.y())
        self._last_mouse_pos = event.pos()
        self.update()
        event.accept()

    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        if delta:
            self.vm_glcore.mouse_scroll(1 if delta > 0 else -1)
            self.update()
        event.accept()

    @staticmethod
    def _button_number(button):
        if button == Qt.LeftButton:
            return 1
        if button == Qt.MiddleButton:
            return 2
        if button == Qt.RightButton:
            return 3
        return 0

    def closeEvent(self, event):
        if self._gl_initialized:
            self.makeCurrent()
            try:
                cleanup = getattr(self.vm_glcore, "cleanup", None)
                if callable(cleanup):
                    cleanup()
            finally:
                self.doneCurrent()
        super().closeEvent(event)
