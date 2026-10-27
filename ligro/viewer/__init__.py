"""Molecular visualization widgets for LiGRO."""

__all__ = ["PyVistaStructureViewer", "create_pyvista_viewer"]


def __getattr__(name):
    if name in __all__:
        from .pyvista_qt import PyVistaStructureViewer, create_pyvista_viewer

        return {
            "PyVistaStructureViewer": PyVistaStructureViewer,
            "create_pyvista_viewer": create_pyvista_viewer,
        }[name]
    raise AttributeError("module {!r} has no attribute {!r}".format(__name__, name))
