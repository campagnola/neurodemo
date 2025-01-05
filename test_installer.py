from PyQt6 import QtWidgets
import pyqtgraph as pg
import numpy as np
import scipy as sp

import sys
app = pg.mkQApp()

class MainWindow(qt.QMainWindow):

    def __init__(self):
        super().__init__()

        self.setWindowTitle("Hello World")
        l = qt.QLabel("My simple app.")
        l.setMargin(10)
        self.setCentralWidget(l)
        self.show()

if __name__ == '__main__':
    app = qt.QApplication(sys.argv)
    w = MainWindow()
    app.exec()