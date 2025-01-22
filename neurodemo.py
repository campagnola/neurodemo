import sys
import platform
from neurodemo import qt
from neurodemo.main_window import DemoWindow

have_splash = False

if getattr(sys, 'frozen', False):
    # If we are running in a PyInstaller bundle, then the 'frozen' attribute will be true.
    # If not, then the attribute will not exist, and getattr() will default to False

    if platform.system() != "Darwin":
        # PyInstaller on MacOS does not support splash screen, so skip this on MacOS

        # On Windows, PyInstaller launches a splash screen to entertain user while everything loads.
        # Import the module that allows us to close the splash screen when main window opens.
        import pyi_splash
        have_splash = True

if __name__ == '__main__':
    if '--dbg' in sys.argv:
        import pyqtgraph as pg
        pg.dbg()

    win = DemoWindow(multiprocessing=False)

    if have_splash:
        # Close splash screen if running from a PyInstaller bundle on Windows
        pyi_splash.close()

    if (sys.flags.interactive != 1) or not hasattr(qt, "PYQT_VERSION"):
        qt.QApplication.instance().exec()
