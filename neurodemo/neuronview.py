import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np


class NeuronView(pg.GraphicsLayoutWidget):
    """Displays a graphical representation of the neuron and its attached
    mechanisms.
    """
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.view = self.addViewBox(0, 0)
        self.soma = QtGui.QGraphicsEllipseItem(QtCore.QRectF(-10, -10, 20, 20))
        self.soma.setPen(pg.mkPen(0.5))
        self.soma2 = QtGui.QGraphicsEllipseItem(QtCore.QRectF(-10.5, -10.5, 21, 21))
        self.soma2.setPen(pg.mkPen(0.5))
        self.view.addItem(self.soma2)
        self.view.addItem(self.soma)
        self.view.setAspectLocked()

    def update(self, vm):
        rgb = np.clip([
            (vm + 30e-3) * 5e3,
            (vm + 65e-3) * 5e3,
            (-65e-3 - vm) * 10e3], 0, 205) + 50
        self.soma.setBrush(pg.mkBrush((rgb[0], rgb[1], rgb[2], 100)))
        
        