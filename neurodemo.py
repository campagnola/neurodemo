import sys
from pyqtgraph.Qt import QtCore, QtWidgets
from neurodemo.main_window import DemoWindow

if getattr(sys, 'frozen', False):
    # If we are running in a PyInstaller bundle, then the 'frozen' attribute will be true.
    # If not, then the attribute will not exist, and getattr() will default to False

    # Import module to allow us to close the splash screen when main window opens.
    import pyi_splash

if __name__ == '__main__':
    if '--dbg' in sys.argv:
        import pyqtgraph as pg
        pg.dbg()

    win = DemoWindow(multiprocessing=False)

    if getattr(sys, 'frozen', False):
        # Close splash screen if running from a PyInstaller bundle
        pyi_splash.close()

    if (sys.flags.interactive != 1) or not hasattr(QtCore, "PYQT_VERSION"):
        QtWidgets.QApplication.instance().exec()
