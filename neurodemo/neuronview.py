# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""
from __future__ import division, unicode_literals
import os, tempfile
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore, QtSvg
from PyQt6 import QtSvgWidgets
from PyQt6 import QtGui as QtGui6
from pyqtgraph import ArrowItem

# for storing dynamically-generated svg files
tmpdir = tempfile.mkdtemp()

def svg_file(name):
    return os.path.join(os.path.dirname(__file__), 'images', name+'.svg')


class NeuronView(pg.GraphicsLayoutWidget):
    """Displays a graphical representation of the neuron and its attached
    mechanisms.
    """
    def __init__(self, soma, mechanisms):
        pg.GraphicsLayoutWidget.__init__(self)
        self.soma = soma
        self.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        self.view = self.addViewBox(0, 0)
        self.view.setAspectLocked()
        self.view.setRange(QtCore.QRectF(-70, -70, 140, 140))
        
        self.cell = Cell(soma)
        self.view.addItem(self.cell)
        self.cell.rotate(90)
        #self.grid = pg.GridItem()
        #self.view.addItem(self.grid)
        
        self.items = [self.cell]
        self.channels = []
        angle = 180
        for mech in mechanisms:
            if mech.type == 'PatchClamp':
                item = Pipette(mech)
                self.view.addItem(item)
                item.setZValue(1)
                item.translate(0, 50)
            else:
                item = Channel(mech)
                item.rotate(angle)
                angle += 45
                item.translate(0, 50)
                self.view.addItem(item)
                self.channels.append(item)
                
            self.items.append(item)
        
        # translucent mask to obscure cell when circuit is visible
        self.mask = QtGui.QGraphicsRectItem(QtCore.QRectF(-1000, -1000, 2000, 2000))
        self.view.addItem(self.mask)
        self.mask.setBrush(pg.mkBrush(0, 0, 0, 180))
        self.mask.setZValue(5)
        
        # circuit items are added separately so they appear above the mask
        self.circuit = QtGui.QGraphicsItemGroup()
        self.view.addItem(self.circuit)
        self.circuit.setZValue(10)
        for i in self.items:
            i.circuit.setParentItem(self.circuit)
        self.show_circuit(False)
        
    def update_state(self, state):
        for item in self.items:
            item.update_state(state)

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
        QtGui.QGraphicsItemGroup.setRotation(self, angle)
        self.circuit.setRotation(angle)
        
    def translate(self, x, y):
        QtGui.QGraphicsItemGroup.setTransform(self, QtGui6.QTransform().translate(x, y))
        self.circuit.setTransform(QtGui6.QTransform().translate(x, y))
        
    def setVisible(self, v):
        QtGui.QGraphicsItemGroup.setVisible(self, v)
        self.circuit.setVisible(v)


class Cell(NeuronItem):
    def __init__(self, section):
        self.key = section.name
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
        vm = state[self.key + '.V']
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
    
    def __init__(self, channel):
        self.channel = channel
        self.key = channel.name + '.OP'
        self.maxop = channel.max_op
        color = {
            'INa': 'dd0000',
            'IK': '0000dd',
            'Ileak': '00dd00',
            'IH': 'aaaaaa',            
        }.get(channel.type, '999999')
        NeuronItem.__init__(self)
        svg = self.get_svg(color)
        self.svg = [QtSvgWidgets.QGraphicsSvgItem(svg),
                    QtSvgWidgets.QGraphicsSvgItem(svg)]
        self.svg[0].setTransform(QtGui6.QTransform().scale(1, -1))
        self.svg[1].setTransform(QtGui6.QTransform().scale(-1, -1))
        for svg in self.svg:
            svg.setParentItem(self)
            svg.setTransform(QtGui6.QTransform().translate(-9, -10.5))
            # svg.translate(-9, -10.5)

        self.bg = QtGui.QGraphicsRectItem(QtCore.QRectF(-5, -10, 10, 20))
        self.bg.setParentItem(self)
        self.bg.setZValue(-1)
        color = pg.mkColor(color)
        self.bg.setBrush(pg.mkBrush(color.red()//2, color.green()//2, color.blue()//2, 255))
        
        self.circuit = QtGui.QGraphicsItemGroup()
        
        self.current = Current(channel.name)
        self.current.setParentItem(self.circuit)
        self.current.setZValue(2)
        
        self.res = Resistor(l1=50, l2=15)
        self.res.setParentItem(self.circuit)
        self.res.setTransform(QtGui6.QTransform().translate(0, -50))
        #translate(0, -50)
    
        self.batt = Capacitor(l1=10, l2=40, w1=11, w2=7, gap=4)
        self.batt.setParentItem(self.circuit)
        self.batt.setTransform(QtGui6.QTransform().translate(0, 15))
        #self.batt.translate(0, 15)
        
    def update_state(self, state):
        try:
            op = state[self.key]
            self.setVisible(True)
        except KeyError:
            # channel is disabled
            self.setVisible(False)
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
        self.scale(amp)
        if self.center:
            self.arrow.setPos(0, -20 * idir)
        else:
            if idir < 0:
                self.arrow.setPos(0, 0)
            else:
                self.arrow.setPos(0, -40)


class Pipette(NeuronItem):
    def __init__(self, clamp):
        self.clamp = clamp
        NeuronItem.__init__(self)
        self.key = clamp.name
        self.svg = QtSvgWidgets.QGraphicsSvgItem(svg_file('pipette'))
        self.svg.setScale(1.0)
        self.svg.setRotation(-180.)
        self.svg.setTransform(QtGui6.QTransform().translate(-50.0, -255.67))
        self.svg.setParentItem(self)
        
        path = QtGui.QPainterPath()
        path.moveTo(42, 4)
        path.lineTo(58, 4)
        path.lineTo(80, 260)
        path.lineTo(14, 260)
        path.closeSubpath()
        
        self.voltage = QtGui.QGraphicsPathItem(path)
        self.voltage.setPen(pg.mkPen(None))
        self.voltage.setTransform(QtGui6.QTransform().translate(-50, -5))
        self.voltage.setParentItem(self)
        self.voltage.setZValue(-1)
        
        self.circuit = QtGui.QGraphicsItemGroup()
        
        self.current = Current(self.key)
        self.current.setParentItem(self.circuit)
        self.current.setZValue(2)
        
        self.res = Resistor(l1=50, l2=150)
        self.res.setParentItem(self.circuit)
        self.res.setTransform(QtGui6.QTransform().translate(0, -50))
        
        self.cap = Capacitor(l1=15, l2=30)
        self.cap.setTransform(QtGui6.QTransform().translate(0, 40))
        self.cap.setRotation(-90.0)
        self.cap.setParentItem(self.circuit)
        
    def update_state(self, state):
        try:
            ve = state[self.key + '.V']
            vm = state['.'.join(self.key.split('.')[:-1]) + '.V']
            self.setVisible(True)
        except KeyError:
            self.setVisible(False)
            return
        y = self.voltage.boundingRect().top()
        grad = QtGui.QLinearGradient(QtCore.QPointF(0, y), QtCore.QPointF(0, y+25))
        grad.setColorAt(0, pg.mkColor(v_color(vm)))
        grad.setColorAt(1, pg.mkColor(v_color(ve)))
        self.voltage.setBrush(QtGui.QBrush(grad))
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
    
    
