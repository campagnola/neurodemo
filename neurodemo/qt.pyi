"""
This stub file informs the IDE about names available in qt.py
"""

from typing import Union, Any

try:
    from PyQt5 import QtCore, QtGui, QtWidgets, QtSvg, QtSvgWidgets
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *
    from PyQt5.QtWidgets import *
    from PyQt5.QtSvg import *
    from PyQt5.QtSvgWidgets import *

    QtCore = QtCore
    QtGui = QtGui
    QtWidgets = QtWidgets
except ImportError:
    try:
        from PyQt6 import QtCore, QtGui, QtWidgets, QtSvg, QtSvgWidgets
        from PyQt6.QtCore import *
        from PyQt6.QtGui import *
        from PyQt6.QtWidgets import *
        from PyQt6.QtSvg import *
        from PyQt6.QtSvgWidgets import *

        QtCore = QtCore
        QtGui = QtGui
        QtWidgets = QtWidgets
    except ImportError:
        try:
            from PySide2 import QtCore, QtGui, QtWidgets, QtSvg, QtSvgWidgets
            from PySide2.QtCore import *
            from PySide2.QtGui import *
            from PySide2.QtWidgets import *
            from PySide2.QtSvg import *
            from PySide2.QtSvgWidgets import *

            QtCore = QtCore
            QtGui = QtGui
            QtWidgets = QtWidgets
        except ImportError:
            try:
                from PySide6 import QtCore, QtGui, QtWidgets, QtSvg, QtSvgWidgets
                from PySide6.QtCore import *
                from PySide6.QtGui import *
                from PySide6.QtWidgets import *
                from PySide6.QtSvg import *
                from PySide6.QtSvgWidgets import *

                QtCore = QtCore
                QtGui = QtGui
                QtWidgets = QtWidgets
            except ImportError as e:
                raise ImportError("No suitable qt binding found") from e

