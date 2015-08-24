import numpy as np
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp
import pyqtgraph.parametertree as pt
from neurodemo.units import *

app = pg.mkQApp()


class DemoWindow(pg.QtGui.QWidget):
    def __init__(self):
        pg.QtGui.QWidget.__init__(self)
        self.resize(1024, 768)
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.ptree = pt.ParameterTree()
        self.layout.addWidget(self.ptree, 0, 0)
        
        self.scroll_plot = ScrollingPlot(parent=self, labels={'left': ('Membrane Potential', 'V'), 
                                                              'bottom': ('Time', 's')})
        self.layout.addWidget(self.scroll_plot, 0, 1)

        self.show()
        
        # Set up remote process to run simulation
        self.proc = mp.Process()
        ndemo = self.proc._import('neurodemo')
        
        self.sim = ndemo.Sim(temp=6.3)
        self.sim._setProxyOptions(deferGetattr=True)
        self.neuron = ndemo.Section()
        self.sim.add(self.neuron)
        self.hhna = self.neuron.add(ndemo.HHNa())
        self.leak = self.neuron.add(ndemo.Leak(gbar=0.1*mS/cm**2))
        self.hhk = self.neuron.add(ndemo.HHK())
        self.clamp = self.neuron.add(ndemo.MultiClamp(mode='ic'))
        
        cmd = np.zeros(10000)
        cmd[2000:6000] = 200e-12
        self.clamp.set_command(cmd, dt=10*us)
        
        self.run_once()
        
        self.timer = pg.QtCore.QTimer()
        self.timer.timeout.connect(self.update)
        self.timer.start(30)
        
    def update(self):
        if self.next_result is None:
            return
        if self.next_result.hasResult():
            res = self.next_result.result()
            res._setProxyOptions(deferGetattr=True)
            self.last_result = res
            self.next_result = None
            vm = res[self.neuron, 'Vm']._getValue()
            self.scroll_plot.append(vm)
            self.run_once()

    def run_once(self):
        self.next_result = self.sim.run(dt=1e-5, dur=100e-3, _callSync='async', _returnType='proxy')


class ScrollingPlot(pg.PlotWidget):
    def __init__(self, **kwds):
        pg.PlotWidget.__init__(self, **kwds)
        self.data_curve = self.plot()
        self.data = []
        
    def append(self, data):
        self.data.append(data)
        while len(self.data) > 10:
            self.data.pop(0)
        data = np.concatenate(self.data)
        t = np.arange(len(data)) * 1e-5
        self.data_curve.setData(t, data)
        
        
        
if __name__ == '__main__':
    import sys
    win = DemoWindow()
    if sys.flags.interactive == 0:
        app.exec_()
