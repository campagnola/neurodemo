import sys
from pyqtgraph.Qt import QtCore, QtWidgets
from neurodemo.main_window import DemoWindow


if __name__ == '__main__':
    proc = None
    # Enable running simulation in background process:
    # import pyqtgraph.multiprocess as mp
    # proc = mp.QtProcess(debug=False)

    win = DemoWindow(proc)
    if (sys.flags.interactive != 1) or not hasattr(QtCore, "PYQT_VERSION"):
        QtWidgets.QApplication.instance().exec()
