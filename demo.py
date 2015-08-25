import numpy as np
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp
import pyqtgraph.parametertree as pt
from neurodemo.units import *

pg.setConfigOption('antialias', True)
app = pg.mkQApp()


class DemoWindow(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.resize(1024, 768)
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.splitter, 0, 0)
        
        self.ptree = pt.ParameterTree()
        self.splitter.addWidget(self.ptree)
        
        self.scroll_plot = ScrollingPlot(parent=self, labels={'left': ('Membrane Potential', 'V'), 
                                                              'bottom': ('Time', 's')})
        self.scroll_plot.setYRange(-90*mV, 50*mV)
        self.scroll_plot.setXRange(0*ms, 800*ms)
        self.splitter.addWidget(self.scroll_plot)
        self.splitter.setSizes([224, 800])

        self.show()
        
        # Set up remote process to run simulation
        self.proc = mp.QtProcess()
        ndemo = self.proc._import('neurodemo')
        self.dt = 10*us
        self.sim = ndemo.Sim(temp=6.3, dt=self.dt)
        self.sim._setProxyOptions(deferGetattr=True)
        self.neuron = ndemo.Section(name='soma')
        self.sim.add(self.neuron)
        self.hhna = self.neuron.add(ndemo.HHNa())
        self.leak = self.neuron.add(ndemo.Leak())
        self.hhk = self.neuron.add(ndemo.HHK())
        self.clamp = self.neuron.add(ndemo.MultiClamp(mode='ic'))
        
        cmd = np.zeros(10000)
        cmd[2000:6000] = 200e-12
        self.clamp.set_command(cmd, dt=10*us)
        
        self.scroll_plot.addLine(y=self.neuron.ek)
        self.scroll_plot.addLine(y=self.neuron.ena)


        self.params = pt.Parameter.create(name='params', type='group', children=[
            dict(name='Run', type='bool', value=True),
            dict(name='Speed', type='float', value=1.0, limits=[0, 1], step=0.1),
            dict(name='Hodgkin Huxley', type='group', children=[
                ChannelParameter(self.leak),
            ]),
        ])
        self.ptree.setParameters(self.params)
        self.params.sigTreeStateChanged.connect(self.params_changed)
        
        self.runner = ndemo.SimRunner(self.sim)
        self.runner.new_result.connect(mp.proxy(self.new_result))
        self.start()

    def params_changed(self, root, changes):
        for param, change, val in changes:
            if change != 'value':
                continue
            if param is self.params.child('Run'):
                if val:
                    self.start()
                else:
                    self.stop()
            elif param is self.params.child('Speed'):
                self.runner.set_speed(val)
            
        
    def start(self):
        self.runner.start()
        
    def stop(self):
        self.runner.stop()
        
    def new_result(self, res):
        self.last_result = res
        vm = res['soma', 'Vm']
        self.scroll_plot.append(vm)


class ChannelParameter(pt.parameterTypes.SimpleParameter):
    def __init__(self, channel):
        self.channel = channel
        pt.parameterTypes.SimpleParameter.__init__(self, name=channel.name, type='bool', 
                                                   value=channel.enabled, children=[
            dict(name='ḡ', type='float', value=0.1, suffix='S/cm²', siPrefix=True, step=0.01),
            dict(name='Erev', type='float', value=-50*mV, suffix='V', siPrefix=True, step=5*mV),
        ])
        self.sigTreeStateChanged.connect(self.treeChange)
        
    def treeChange(self, root, changes):
        for param, change, val in changes:
            if change != 'value':
                continue
            if param is self:
                self.channel.enabled = val
            elif param is self.child('ḡ'):
                self.channel.gbar = val
            elif param is self.child('Erev'):
                self.channel.erev = val


class ScrollingPlot(pg.PlotWidget):
    def __init__(self, **kwds):
        pg.PlotWidget.__init__(self, **kwds)
        self.data_curve = self.plot()
        self.data = []
        
    def append(self, data):
        self.data.append(data)
        while len(self.data) > 150:
            self.data.pop(0)
        data = np.concatenate(self.data)
        t = np.arange(len(data)) * 1e-5
        self.data_curve.setData(t, data)
        

        
if __name__ == '__main__':
    import sys
    win = DemoWindow()
    if sys.flags.interactive == 0:
        app.exec_()
