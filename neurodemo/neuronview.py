import os, tempfile
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore, QtSvg
import numpy as np


# for storing dynamically-generated svg files
tmpdir = tempfile.mkdtemp()

def svg_file(name):
    return os.path.join(os.path.dirname(__file__), 'images', name+'.svg')


class NeuronView(pg.GraphicsLayoutWidget):
    """Displays a graphical representation of the neuron and its attached
    mechanisms.
    """
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.setRenderHint(QtGui.QPainter.Antialiasing)
        self.view = self.addViewBox(0, 0)
        self.view.invertY()  # because svg +y points downward
        self.soma = QtGui.QGraphicsEllipseItem(QtCore.QRectF(-10, -10, 20, 20))
        self.soma.setPen(pg.mkPen(0.5))
        self.soma2 = QtGui.QGraphicsEllipseItem(QtCore.QRectF(-10.5, -10.5, 21, 21))
        self.soma2.setPen(pg.mkPen(0.5))
        self.view.addItem(self.soma2)
        self.view.addItem(self.soma)
        self.view.setAspectLocked()
        
        #self.grid = pg.GridItem()
        #self.view.addItem(self.grid)
        
        self.svg = QtSvg.QGraphicsSvgItem(svg_file('cell'))
        self.svg.translate(-200, -150)
        self.view.addItem(self.svg)
        
        self.channels = []
        angle = 0
        for key, color in [('soma.INa', 'dd0000'), 
                           ('soma.IK', '0000dd')]:
            channel = Channel(key, color)
            channel.rotate(angle)
            angle += 30
            channel.translate(0, 50)
            self.view.addItem(channel)
            self.channels.append(channel)
        
    def update_state(self, state):
        vm = state['soma.Vm']
        rgb = np.clip([
            (vm + 65e-3) * 5e3,
            (vm + 30e-3) * 5e3,
            (-65e-3 - vm) * 10e3], 0, 205) + 50
        self.soma.setBrush(pg.mkBrush((rgb[0], rgb[1], rgb[2], 100)))
        
        for ch in self.channels:
            ch.update_state(state)

        
class Channel(QtGui.QGraphicsItemGroup):
    
    @staticmethod
    def get_svg(color):
        """Return a new copy of the channel SVG with the color modified.
        """
        svg = open(svg_file('channel')).read()
        svg = svg.replace('fill:#dc0000', 'fill:#'+color)
        fname = os.path.join(tmpdir, 'channel_' + color + '.svg')
        open(fname, 'w').write(svg)
        return fname
    
    def __init__(self, key, color):
        self.key = key
        QtGui.QGraphicsItemGroup.__init__(self)
        svg = self.get_svg(color)
        self.svg = [QtSvg.QGraphicsSvgItem(svg),
                    QtSvg.QGraphicsSvgItem(svg)]
        self.svg[1].scale(-1, 1)
        
        for svg in self.svg:
            svg.setParentItem(self)
            svg.translate(-9, -10.5)
            
        self.bg = QtGui.QGraphicsRectItem(QtCore.QRectF(-5, -10, 10, 20))
        self.bg.setParentItem(self)
        self.bg.setZValue(-1)
        self.bg.setBrush(pg.mkBrush(100, 0, 0))
    
    def update_state(self, state):
        op = state[self.key + '.OP']
        self.svg[0].setPos(op * -15, 0)
        self.svg[1].setPos(op * 15, 0)
            

