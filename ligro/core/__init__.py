"""LiGRO core: simulation logic with no GUI imports.

Nothing under ligro.core may import PyQt5, GTK or any other UI toolkit.
That rule is what makes the pipelines testable with a FakeRunner.
"""
