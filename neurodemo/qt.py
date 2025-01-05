import importlib
import pyqtgraph as pg

Qt = importlib.__import__(pg.Qt.QT_LIB + '.QtSvgWidgets')

# make one large namespace containing everything; pyqtgraph handles translation
# between different Qt versions
for mod in [pg.Qt, pg.Qt.QtCore, pg.Qt.QtGui, pg.Qt.QtWidgets, pg.Qt.QtSvg, Qt.QtSvgWidgets]:
    ns = mod.__dict__.copy()
    # don't copy special variables like __name__, __file__, etc.
    for k in list(ns.keys()):
        if k.startswith('__'):
            ns.pop(k)
    globals().update(ns)
