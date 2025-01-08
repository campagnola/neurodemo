import importlib
import pyqtgraph as pg

if pg.Qt.QT_LIB == 'PyQt6':
    Qt = importlib.__import__(pg.Qt.QT_LIB + '.QtSvgWidgets')
    mod_list = [pg.Qt, pg.Qt.QtCore, pg.Qt.QtGui, pg.Qt.QtWidgets, pg.Qt.QtSvg, Qt.QtSvgWidgets]
else:
    # PyQt5 doesn't need the extra import
    mod_list = [pg.Qt, pg.Qt.QtCore, pg.Qt.QtGui, pg.Qt.QtWidgets, pg.Qt.QtSvg]

# make one large namespace containing everything; pyqtgraph handles translation
# between different Qt versions
for mod in mod_list:
    ns = mod.__dict__.copy()
    # don't copy special variables like __name__, __file__, etc.
    for k in list(ns.keys()):
        if k.startswith('__'):
            ns.pop(k)
    globals().update(ns)
