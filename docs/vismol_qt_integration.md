# VisMol Qt integration

LiGRO now vendors a Qt-host adapter for the upstream VisMol `graphics_engine`.
The molecular model, parser, shaders and `VismolGLCore` are reused; the GTK
`VismolGTKWidget` is not embedded or imported by the Qt viewer.

## Install

From the LiGRO project root, in the environment used by the application:

```bash
python -m pip install -e vendor/graphics_engine
```

The graphics engine's setup declares PyQt5 and PyOpenGL dependencies. Install
the remaining runtime requirements listed by the engine/project if they are
not already present.

## Use inside a PyQt5 window

```python
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout
from ligro.viewer import create_vismol_viewer

window = QMainWindow()
central = QWidget(window)
layout = QVBoxLayout(central)
session, viewer = create_vismol_viewer(central)
layout.addWidget(viewer)
window.setCentralWidget(central)

# Load after the widget/session are created:
session.load_molecule("/absolute/path/to/structure.pdb")
window.show()
```

`VismolSession(toolkit="Qt5")` also creates this widget and exposes it as
`session.vm_widget`. The widget's `initializeGL`, `resizeGL`, and `paintGL`
methods drive the existing `VismolGLCore` with Qt's current OpenGL context.
Mouse rotation/pan/zoom and wheel zoom are bridged to the core.

## Notes

* Requesting OpenGL 3.3 Core is intentional because VisMol's shaders use GLSL
  3.30. The host platform must provide a compatible OpenGL implementation.
* EasyHybrid-only builder hooks are optional; when EasyHybrid's
  `gui.windows.builder` package is unavailable, no-op hooks keep ordinary
  viewing independent of builder functionality.
* GTK and Qt must not share a widget or GL context. Do not instantiate
  `VismolGTKWidget` from LiGRO's Qt UI.
* The adapter is syntax-checked, but an actual interactive OpenGL smoke test
  requires a desktop/display server and a working OpenGL driver.
