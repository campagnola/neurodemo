"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015

Updated 2022 pbmanis for Python 3.10 and PyQt6


"""
import os, tempfile
from pathlib import Path
from typing import Union, List, Tuple
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore, QtSvg
from PyQt6 import QtSvgWidgets
from PyQt6 import QtGui as QtGui6
from pyqtgraph import ArrowItem

# for storing dynamically-generated svg files
tmpdir = tempfile.mkdtemp()


def svg_file(name):
    return Path(Path(__file__).parent, "images", name + ".svg")


class NeuronView(pg.GraphicsLayoutWidget):
    """Displays a graphical representation of the neuron and its attached
    mechanisms.


    """

    def __init__(self, soma, mechanisms):
        """Create a NeuronView instance

        Parameters
        ----------
        soma : neuronsim object
            The soma object, with mechanisms methods.
        mechanisms : _type_
            _description_
        """
        pg.GraphicsLayoutWidget.__init__(self)
        self.soma = soma
        self.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        self.view = self.addViewBox(0, 0)
        self.view.setAspectLocked()
        self.view.setRange(QtCore.QRectF(-70, -70, 140, 140))

        self.cell = Cell(soma)
        self.view.addItem(self.cell)
        self.cell.setRotation(90)
        # self.grid = pg.GridItem()
        # self.view.addItem(self.grid)

        self.items = [self.cell]
        self.channels = []
        angle = 0.0
        nmech = len(mechanisms)
        anglestep = 360./nmech
        for mech in mechanisms:
            if mech.type == "PatchClamp":
                item = Pipette(mech)
                self.view.addItem(item)
                item.setZValue(1)
                item.setTransform(QtGui6.QTransform().translate(0, 50))
            else:  # add ion channels
                item = Channel(
                    mech, translation=(0, 50), angle=angle
                )  # set all transforms at once
                self.view.addItem(item)  # add to graphic view of the cell
                # item.setTransform(QtGui6.QTransform().rotate(angle))
                angle += anglestep
                if 195 > angle > 165:
                    angle += anglestep  # skip region aroun pipette
                # item.setTransform(QtGui6.QTransform().translate(0, 50))
                self.channels.append(item)  # keep a list

            self.items.append(item)  # all of the items in the view

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
    """Return a color corresponding to voltage."""
    return (
        np.clip(
            [(v + 65e-3) * 5e3, (v + 30e-3) * 5e3, (-65e-3 - v) * 10e3, 255], 0, 205
        )
        + 50
    )


class NeuronItem(QtGui.QGraphicsItemGroup):
    def __init__(self):
        QtGui.QGraphicsItemGroup.__init__(self)
        self.circuit = None
        self.current = None

    def rotate(self, angle):
        QtGui.QGraphicsItemGroup.setRotation(
            self, angle
        )  # QtGui6.QTransform().rotate(angle)) #  angle)
        self.circuit.setRotation(angle)

    def translate(self, x, y):
        QtGui.QGraphicsItemGroup.setTransform(
            self, QtGui6.QTransform().fromTranslate(x, y)
        )
        self.circuit.setTransform(QtGui6.QTransform().fromTranslate(x, y))

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
        vm = state[self.key + ".V"]
        self.soma.setBrush(pg.mkBrush(v_color(vm)))
        self.current.update_state(state)

    def show_circuit(self, show):
        self.cap.setVisible(show)


class Channel(NeuronItem):
    """Draw a channel on the cell membrane
    The channel is an svg. Two of the svg's, mirrored
    are needed to make the channel, then they are rotated
    and placed.
    Each channel has an associated current object (arrow),
    and possibly a conductance (resistor) and battery ("cap")

    """

    @staticmethod
    def get_svg(color):
        """Return a new copy of the channel SVG with the color modified."""
        with open(svg_file("channel2")) as fh:
            svg = fh.read()
        svg = svg.replace("fill:#dc0000", "fill:#" + color)
        fname = Path(tmpdir, f"channel_{color:s}.svg")
        with open(fname, "w") as fh:
            fh.write(svg)
        return str(fname)

    def __init__(self, channel, translation: Tuple = (0.0, 0.0), angle: float = 0.0):
        self.channel = channel
        self.key = channel.name + ".OP"
        self.maxop = channel.max_op
        self.translation = translation
        self.angle = angle
        self.svg_items = [None, None]
        color = {
            "INa": "dd0000",  # HH
            "IK": "0000dd",   # HH
            "Ileak": "00dd00",  # HH
            "IH": "aa00aa",  # Not in HH model, but included here
            "INa1": "880000",  # LG
            "IKf": "0088ff", # LG
            "IKs": "8800ff", # LG
        }.get(channel.type, "999999")
        NeuronItem.__init__(self)
        self.channel_svg = [
            QtSvgWidgets.QGraphicsSvgItem(self.get_svg(color)),
            # QtSvgWidgets.QGraphicsSvgItem(self.get_svg(color)),   # not needed with channel2.svg
        ]
        scale = [[1, -1], [1, -1]]
        transl = list(self.translation)
        transl[0] -= 9.066
        for i, svg in enumerate(self.channel_svg):
            svg.setParentItem(self)
            self.svg_items[i] = transform = self.set_transform(
                svg, scale[i], angle=self.angle, translate=transl
            )

        self.bg = QtGui.QGraphicsRectItem(QtCore.QRectF(-5, -10, 10, 20))
        self.bg.setParentItem(self)
        self.bg.setZValue(-1)
        color = pg.mkColor(color)
        self.bg.setBrush(
            pg.mkBrush(color.red() // 2, color.green() // 2, color.blue() // 2, 255)
        )

        self.circuit = QtGui.QGraphicsItemGroup()

        self.current = Current(channel.name, center=False)
        self.current.setParentItem(self.circuit)
        self.current_translation = self.translation
        self.current_angle = self.angle
        self.current_transform = self.set_transform(
            self.current,
            scale=scale[0],
            angle=self.current_angle,
            translate=self.current_translation,
        )
        self.current.setZValue(2)

        self.res = Resistor(l1=50, l2=15)
        self.res.setParentItem(self.circuit)
        self.set_transform(self.res, [1, 1], angle=0.0, translate=[0, -50])

        self.batt = Battery(l1=10, l2=40, w1=11, w2=7, gap=4)
        self.batt.setParentItem(self.circuit)
        self.set_transform(self.batt, [1, 1], angle=0.0, translate=[0, 15])
        self.scale = [1, 1]
        self.I_angle = 0.0
        self.I_translation = [0, 15]

    # def update_currents(self, scale):
    #     """Update the current arrow associated with this channel

    #     Parameters
    #     ----------
    #     scale : float
    #         scaled current
    #     """
    #     pass
    # print("update currents: in Channels")
    # self.set_transform(self.current,
    #     scale = scale,
    #     angle = self.current_angle,
    #     translate = self.current_translation)
    # new_transform = QtGui6.QTransform()
    # new_transform.scale(1, scale)
    # new_transform.rotate(self.angle)
    # new_transform.translate(self.translation[0], self.translation[1])

    # self.current.setTransform(new_transform)

    def set_transform(
        self,
        item: object,
        scale: Union[List, Tuple] = [1, 1],
        angle: float = 0.0,
        translate: Union[List, Tuple] = (0, 50),
    ):
        new_transform = QtGui6.QTransform()
        new_transform.scale(scale[0], scale[1])
        new_transform.rotate(angle)
        new_transform.translate(translate[0], translate[1])
        item.setTransform(new_transform)

    def update_state(self, state):
        try:
            op = state[self.key]
            self.setVisible(True)
        except KeyError:
            # channel is disabled
            self.setVisible(False)
            return
        nop = np.clip(op / self.maxop, 0, 1)
        # self.channel_svg[0].setPos(nop * -4, 0)
        # self.svg[1].setPos(nop * 4, 0)
        self.current.update_state(state)

    def show_circuit(self, show):
        self.batt.setVisible(show)
        self.res.setVisible(show)


class Current(QtGui.QGraphicsItemGroup):
    def __init__(self, channel_name, color="y", center=False):
        QtGui.QGraphicsItemGroup.__init__(self)
        self.key = channel_name + ".I"
        self.center = center
        self.length = 20
        self.arrow = ArrowItem(
            brush=(255, 255, 0, 200),
            angle=90,
            tailLen=self.length,
            tailWidth=4,
            headLen=15,
            headWidth=8,
            pxMode=False,
            pen={"color": "k", "width": 1, "cosmetic": False},
        )
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
        if amp < 1e-12:
            amp = 0.
        self.arrow.setStyle(
            angle=90 * idir,
            tailLen=self.length,
            tailWidth=4,
            headLen=15,
            headWidth=8,
        )
        self.setScale(amp)
        if self.center:
            self.arrow.setPos(0, -self.length * idir)
        else:
            if idir < 0:
                self.arrow.setPos(0, 0)
            else:
                self.arrow.setPos(0, -2*self.length)


class Pipette(NeuronItem):
    def __init__(self, clamp):
        self.clamp = clamp
        NeuronItem.__init__(self)
        self.key = clamp.name
        self.pipette_svg = QtSvgWidgets.QGraphicsSvgItem(str(svg_file("pipette")))
        pip_transform = QtGui6.QTransform()
        pip_transform.scale(1.0, 1.0)  # put pipette at the top of the cell
        pip_transform.rotate(-180.0)
        pip_transform.translate(-50.0, -255.67)
        self.pipette_svg.setTransform(pip_transform)
        self.pipette_svg.setParentItem(self)

        path = QtGui.QPainterPath()
        path.moveTo(42, 4)
        path.lineTo(58, 4)
        path.lineTo(80, 260)
        path.lineTo(14, 260)
        path.closeSubpath()

        self.voltage = QtGui.QGraphicsPathItem(path)
        self.voltage.setParentItem(self)
        self.voltage.setPen(pg.mkPen(None))
        voltage_transform = QtGui6.QTransform()
        voltage_transform.translate(-50, -5)
        self.voltage.setTransform(voltage_transform)
        self.voltage.setZValue(-1)

        self.circuit = QtGui.QGraphicsItemGroup()

        self.current = Current(self.key)
        self.current.setParentItem(self.circuit)
        self.current.setZValue(2)

        self.res = Resistor(l1=50, l2=150)
        self.res.setParentItem(self.circuit)
        res_transform = QtGui6.QTransform()
        res_transform.translate(0, -50)
        self.res.setTransform(res_transform)

        self.cap = Capacitor(l1=15, l2=30)
        self.cap.setParentItem(self.circuit)
        cap_transform = QtGui6.QTransform()
        cap_transform.translate(0, 40)
        cap_transform.rotate(-90.0)
        self.cap.setTransform(cap_transform)

    def update_state(self, state):
        try:
            ve = state[self.key + ".V"]
            vm = state[".".join(self.key.split(".")[:-1]) + ".V"]
            self.setVisible(True)
        except KeyError:
            self.setVisible(False)
            return
        y = self.voltage.boundingRect().top()
        grad = QtGui.QLinearGradient(QtCore.QPointF(0, y), QtCore.QPointF(0, y + 25))
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
        path.lineTo(0, l1 - g2)
        path.moveTo(-w1 / 2, l1 - g2)
        path.lineTo(w1 / 2, l1 - g2)
        path.moveTo(-w2 / 2, l1 + g2)
        path.lineTo(w2 / 2, l1 + g2)
        path.moveTo(0, l1 + g2)
        path.lineTo(0, l1 + l2)

        self.line = QtGui.QGraphicsPathItem(path)
        self.line.setBrush(pg.mkBrush(None))
        self.line.setPen(pg.mkPen("w", width=1, cosmetic=False))
        self.line.setParentItem(self)


class Battery(QtGui.QGraphicsItemGroup):
    def __init__(self, l1, l2, w1=10, w2=10, gap=4):
        QtGui.QGraphicsItemGroup.__init__(self)

        g2 = gap / 2
        path = QtGui.QPainterPath()
        path.moveTo(0, 0)
        path.lineTo(0, l1 + g2)
        path.moveTo(-w1 / 2, l1 - g2)
        path.lineTo(w1 / 2, l1 - g2)
        path.moveTo(-w2 / 2, l1 + g2)
        path.lineTo(w2 / 2, l1 + g2)
        path.moveTo(0, l1 + g2)
        path.lineTo(0, l1 + l2)

        self.line = QtGui.QGraphicsPathItem(path)
        self.line.setBrush(pg.mkBrush(None))
        self.line.setPen(pg.mkPen("w", width=1, cosmetic=False))
        self.line.setParentItem(self)


class Resistor(QtGui.QGraphicsItemGroup):
    def __init__(self, l1, l2):
        QtGui.QGraphicsItemGroup.__init__(self)

        w = 3
        h = w / 3**0.5

        path = QtGui.QPainterPath()
        path.moveTo(0, 0)
        y = l1 - h * 6
        path.lineTo(0, y)
        y += h
        for i in range(3):
            path.lineTo(-w, y)
            path.lineTo(w, y + h * 2)
            y += h * 4
        path.lineTo(0, y - h)
        path.lineTo(0, l1 + l2)

        self.line = QtGui.QGraphicsPathItem(path)
        self.line.setBrush(pg.mkBrush(None))
        self.line.setPen(pg.mkPen("w", width=1, cosmetic=False))
        self.line.setParentItem(self)
