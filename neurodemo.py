import sys
from pyqtgraph.Qt import QtCore, QtWidgets
from neurodemo.main_window import DemoWindow


if __name__ == '__main__':
    if '--dbg' in sys.argv:
        import pyqtgraph as pg
        pg.dbg()

    win = DemoWindow(multiprocessing=False)
    if (sys.flags.interactive != 1) or not hasattr(QtCore, "PYQT_VERSION"):
        QtWidgets.QApplication.instance().exec()
