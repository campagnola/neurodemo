import numpy as np
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp
import pyqtgraph.parametertree as pt
from units import *

app = pg.mkQApp()


class DemoWindow(pg.QtGui.QWidget):
    def __init__(self):
        pg.QtGui.QWidget.__init__(self)
        self.resize(1024, 768)
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.ptree = pt.ParameterTree()
        self.layout.addWidget(self.ptree, 0, 0)
        
        self.scroll_plot = pg.PlotWidget(parent=self, labels={'left': ('Membrane Potential', 'V'), 
                                                              'bottom': ('Time', 's')})
        self.layout.addWidget(self.scroll_plot, 0, 1)

        self.show()
        
        # Set up remote process to run simulation
        self.proc = mp.Process()
        ndemo = self.proc._import('neurodemo')
        
        self.sim = ndemo.Sim(temp=6.3)
        self.neuron = ndemo.Section()
        self.hhna = self.neuron.add(ndemo.HHNa())
        self.leak = self.neuron.add(ndemo.Leak(gbar=0.1*mS/cm**2))
        self.hhk = self.neuron.add(ndemo.HHK())
        self.clamp = self.neuron.add(ndemo.MultiClamp(mode='ic'))
        
        cmd = np.zeros(10000)
        cmd[2000:6000] = 200e-12
        self.clamp.set_command(cmd)
        
        self.run_once()
        
        self.timer = pg.QtCore.QTimer()
        self.timer.timeout.connect(self.update)
        
    def update(self):
        if self.next_result is None:
            return
        if self.next_result.hasResult():
            res = self.next_result.result()
            self.last_result = res
            self.next_result = None
            self.scroll_plot.plot(res['t'], res[self.neuron, 'Vm'], clear=True)
            self.run_once()

    def run_once(self):
        self.next_result = self.sim.run(dt=1e-5, dur=100e-3, _sync='async')
        
        
if __name__ == '__main__':
    import sys
    win = DemoWindow()
    if sys.flags.interactive == 0:
        app.exec_()
