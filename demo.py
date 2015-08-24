import numpy as np
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp
import pyqtgraph.parametertree as pt
from neurodemo.units import *

pg.setConfigOption('antialias', True)
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
        self.scroll_plot.setYRange(-90*mV, 50*mV)
        self.layout.addWidget(self.scroll_plot, 0, 1)

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
        
        self.runner = ndemo.SimRunner(self.sim)
        self.runner.new_result.connect(mp.proxy(self.new_result))
        self.runner.start()
        
        
        #self.run_once()
        #self.timer = pg.QtCore.QTimer()
        #self.timer.timeout.connect(self.update)
        #self.timer.start(10)
        
    #def update(self):
        #if self.next_result is None:
            #return
        #if self.next_result.hasResult():
            #nr = self.next_result
            #self.run_once()
            #res = nr.result()
            #res._setProxyOptions(deferGetattr=True)
            #self.last_result = res
            #vm = res[self.neuron, 'Vm']._getValue()
            #self.scroll_plot.append(vm)

    #def run_once(self):
        #self.next_result = self.sim.run(samples=200, _callSync='async', _returnType='proxy')

    def new_result(self, res):
        #res._setProxyOptions(deferGetattr=True)
        self.last_result = res
        vm = res['soma', 'Vm']
        self.scroll_plot.append(vm)
        


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
