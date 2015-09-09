from __future__ import division
import os, tempfile
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore, QtSvg
from .arrow import ArrowItem

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
        self.view.setAspectLocked()
        self.view.setRange(QtCore.QRectF(-70, -70, 140, 140))
        
        self.cell = Cell('soma')
        self.view.addItem(self.cell)
        self.cell.rotate(90)
        #self.grid = pg.GridItem()
        #self.view.addItem(self.grid)
        
        self.pipette = Pipette('soma.PatchClamp')
        self.view.addItem(self.pipette)
        self.pipette.setZValue(1)
        self.pipette.translate(0, 50)
        
        self.channels = []
        angle = 180
        chan = [
            ('soma.INa', 'dd0000', 0.1), 
            ('soma.IK', '0000dd', 0.1),
            ('soma.Ileak', '00dd00', 1.0),
            ('soma.IH', 'aa00aa', 0.05)
        ]
        for key, color, maxop in chan:
            channel = Channel(key, color, maxop)
            channel.rotate(angle)
            angle += 45
            channel.translate(0, 50)
            self.view.addItem(channel)
            self.channels.append(channel)
        
        self.items = [self.cell, self.pipette] + self.channels
        
        self.mask = QtGui.QGraphicsRectItem(QtCore.QRectF(-1000, -1000, 2000, 2000))
        self.view.addItem(self.mask)
        self.mask.setBrush(pg.mkBrush(0, 0, 0, 180))
        self.mask.setZValue(5)
        
        self.circuit = QtGui.QGraphicsItemGroup()
        self.view.addItem(self.circuit)
        self.circuit.setZValue(10)
        for i in self.items:
            i.circuit.setParentItem(self.circuit)
        self.show_circuit(False)
        
    def update_state(self, state):
        self.cell.update_state(state)
        self.pipette.update_state(state)
        for ch in self.channels:
            ch.update_state(state)

    def show_circuit(self, show):
        self.mask.setVisible(show)
        for i in self.items:
            i.show_circuit(show)


def v_color(v):
    """Return a color corresponding to voltage.
    """
    return np.clip([
        (v + 65e-3) * 5e3,
        (v + 30e-3) * 5e3,
        (-65e-3 - v) * 10e3,
        255], 0, 205) + 50


class NeuronItem(QtGui.QGraphicsItemGroup):
    def __init__(self):
        QtGui.QGraphicsItemGroup.__init__(self)
        self.circuit = None
        self.current = None
        
    def rotate(self, angle):
        QtGui.QGraphicsItemGroup.rotate(self, angle)
        self.circuit.rotate(angle)
        
    def translate(self, x, y):
        QtGui.QGraphicsItemGroup.translate(self, x, y)
        self.circuit.translate(x, y)
        


class Cell(NeuronItem):
    def __init__(self, cell_name):
        self.key = cell_name
        NeuronItem.__init__(self)

        self.soma = QtGui.QGraphicsEllipseItem(QtCore.QRectF(-50, -50, 100, 100))
        self.soma.setPen(pg.mkPen(0.5, width=1, cosmetic=False))
        self.soma.setParentItem(self)
        self.soma2 = QtGui.QGraphicsEllipseItem(QtCore.QRectF(-52, -52, 104, 104))
        self.soma2.setPen(pg.mkPen(0.5, width=1, cosmetic=False))
        self.soma2.setParentItem(self)
        
        self.circuit = QtGui.QGraphicsItemGroup()
        
        self.current = Current(self.key, center=False)
        self.current.setParentItem(self.circuit)
        self.current.setPos(0, 50)
        self.current.setZValue(2)
        
        self.cap = Capacitor(l1=50, l2=50)
        self.cap.setParentItem(self.circuit)
        self.cap.setZValue(1)
        

    def update_state(self, state):
        vm = state[self.key + '.Vm']
        self.soma.setBrush(pg.mkBrush(v_color(vm)))
        self.current.update_state(state)

    def show_circuit(self, show):
        self.cap.setVisible(show)

        
class Channel(NeuronItem):
    
    @staticmethod
    def get_svg(color):
        """Return a new copy of the channel SVG with the color modified.
        """
        svg = open(svg_file('channel')).read()
        svg = svg.replace('fill:#dc0000', 'fill:#'+color)
        fname = os.path.join(tmpdir, 'channel_' + color + '.svg')
        open(fname, 'w').write(svg)
        return fname
    
    def __init__(self, channel_name, color, maxop):
        self.key = channel_name + '.OP'
        self.maxop = maxop
        NeuronItem.__init__(self)
        svg = self.get_svg(color)
        self.svg = [QtSvg.QGraphicsSvgItem(svg),
                    QtSvg.QGraphicsSvgItem(svg)]
        self.svg[0].scale(1, -1)
        self.svg[1].scale(-1, -1)
        
        for svg in self.svg:
            svg.setParentItem(self)
            svg.translate(-9, -10.5)
            
        self.bg = QtGui.QGraphicsRectItem(QtCore.QRectF(-5, -10, 10, 20))
        self.bg.setParentItem(self)
        self.bg.setZValue(-1)
        color = pg.mkColor(color)
        self.bg.setBrush(pg.mkBrush(color.red()//2, color.green()//2, color.blue()//2, 255))
        
        self.circuit = QtGui.QGraphicsItemGroup()
        
        self.current = Current(channel_name)
        self.current.setParentItem(self.circuit)
        self.current.setZValue(2)
        
        self.batt = Capacitor(l1=30, l2=5, w1=7, w2=11, gap=3)
        self.batt.setParentItem(self.circuit)
        self.batt.translate(0, -50)
        
        self.res = Resistor(l1=15, l2=50)
        self.res.setParentItem(self.circuit)
        self.res.translate(0, -15)
    
    def update_state(self, state):
        try:
            op = state[self.key]
            self.setVisible(True)
            self.circuit.setVisible(True)
        except KeyError:
            # channel is disabled
            self.setVisible(False)
            self.circuit.setVisible(False)
            return
        nop = np.clip(op / self.maxop, 0, 1)
        self.svg[0].setPos(nop * -4, 0)
        self.svg[1].setPos(nop * 4, 0)
        self.current.update_state(state)

    def show_circuit(self, show):
        self.batt.setVisible(show)
        self.res.setVisible(show)


class Current(QtGui.QGraphicsItemGroup):
    def __init__(self, channel_name, color='y', center=True):
        QtGui.QGraphicsItemGroup.__init__(self)
        self.key = channel_name + '.I'
        self.center = center
        self.arrow = ArrowItem(brush=(255, 255, 0, 200), angle=90, tailLen=20, pxMode=False, 
                               pen={'color': 'k', 'width': 1, 'cosmetic': False})
        self.arrow.setParentItem(self)

    def update_state(self, state):
        I = state[self.key]
        if abs(I) < 100e-15:
            self.setVisible(False)
            return
        else:
            self.setVisible(True)
        idir = 1 if I > 0 else -1
        amp = np.clip((abs(I) / 10e-9), 0, 1) ** 0.2  # normalize to 10 nA range
        self.arrow.setStyle(
            angle=90 * idir,
            tailLen=20,
            tailWidth=5,
            headLen=20,
            headWidth=10,
        )
        self.resetTransform()
        self.scale(amp, amp)
        if self.center:
            self.arrow.setPos(0, -20 * idir)
        else:
            if idir < 0:
                self.arrow.setPos(0, 0)
            else:
                self.arrow.setPos(0, -40)


class Pipette(NeuronItem):
    def __init__(self, key, color='y'):
        NeuronItem.__init__(self)
        self.key = key
        self.svg = QtSvg.QGraphicsSvgItem(svg_file('pipette'))
        self.svg.scale(1, -1)
        self.svg.translate(-50, -255.67)
        self.svg.setParentItem(self)
        
        path = QtGui.QPainterPath()
        path.moveTo(42, 4)
        path.lineTo(58, 4)
        path.lineTo(80, 260)
        path.lineTo(14, 260)
        path.closeSubpath()
        
        self.voltage = QtGui.QGraphicsPathItem(path)
        self.voltage.setPen(pg.mkPen(None))
        self.voltage.translate(-50, -5)
        self.voltage.setParentItem(self)
        self.voltage.setZValue(-1)
        
        self.circuit = QtGui.QGraphicsItemGroup()
        
        self.current = Current(key)
        self.current.setParentItem(self.circuit)
        self.current.setZValue(2)
        
        self.res = Resistor(l1=50, l2=150)
        self.res.setParentItem(self.circuit)
        self.res.translate(0, -50)
        
        self.cap = Capacitor(l1=15, l2=30)
        self.cap.translate(0, 40)
        self.cap.rotate(-90)
        self.cap.setParentItem(self.circuit)
        
    def update_state(self, state):
        ve = state[self.key + '.Ve']
        self.voltage.setBrush(pg.mkBrush(v_color(ve)))
        self.current.update_state(state)

    def show_circuit(self, show):
        self.cap.setVisible(show)
        self.res.setVisible(show)


class Capacitor(QtGui.QGraphicsItemGroup):
    def __init__(self, l1, l2, w1=10, w2=10, gap=6):
        QtGui.QGraphicsItemGroup.__init__(self)
        
        g2 = gap / 2
        path = QtGui.QPainterPath()
        path.moveTo(0, 0)
        path.lineTo(0, l1-g2)
        path.moveTo(-w1/2, l1-g2)
        path.lineTo(w1/2, l1-g2)
        path.moveTo(-w2/2, l1+g2)
        path.lineTo(w2/2, l1+g2)
        path.moveTo(0, l1+g2)
        path.lineTo(0, l1+l2)
        
        self.line = QtGui.QGraphicsPathItem(path)
        self.line.setBrush(pg.mkBrush(None))
        self.line.setPen(pg.mkPen('w', width=1, cosmetic=False))
        self.line.setParentItem(self)


class Resistor(QtGui.QGraphicsItemGroup):
    def __init__(self, l1, l2):
        QtGui.QGraphicsItemGroup.__init__(self)
        
        w = 3
        h = w / 3**0.5
        
        path = QtGui.QPainterPath()
        path.moveTo(0, 0)
        y = l1 - h*6
        path.lineTo(0, y)
        y += h
        for i in range(3):
            path.lineTo(-w, y)
            path.lineTo(w, y+h*2)
            y += h*4
        path.lineTo(0, y-h)
        path.lineTo(0, l1+l2)
        
        self.line = QtGui.QGraphicsPathItem(path)
        self.line.setBrush(pg.mkBrush(None))
        self.line.setPen(pg.mkPen('w', width=1, cosmetic=False))
        self.line.setParentItem(self)
    
    
