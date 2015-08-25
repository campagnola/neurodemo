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
        # set up simulation in remote process
        self.dt = 100*us
        self.proc = mp.QtProcess()
        ndemo = self.proc._import('neurodemo')
        self.sim = ndemo.Sim(temp=6.3, dt=self.dt)
        self.sim._setProxyOptions(deferGetattr=True)
        self.neuron = ndemo.Section(name='soma')
        self.sim.add(self.neuron)
        
        self.hhna = self.neuron.add(ndemo.HHNa())
        self.leak = self.neuron.add(ndemo.Leak())
        self.hhk = self.neuron.add(ndemo.HHK())
        self.dexh = self.neuron.add(ndemo.IH())
        self.dexh.enabled = False
        
        self.clamp = self.neuron.add(ndemo.MultiClamp(mode='ic'))
        cmd = np.zeros(int(1/self.dt))
        cmd[int(0.4/self.dt):int(0.8/self.dt)] = 200e-12
        self.clamp.set_command(cmd, dt=self.dt)

        
        # set up GUI
        QtGui.QWidget.__init__(self)
        self.resize(1024, 768)
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.splitter, 0, 0)
        
        self.ptree = pt.ParameterTree()
        self.splitter.addWidget(self.ptree)
        
        self.plot_splitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.splitter.addWidget(self.plot_splitter)
        
        self.vm_plot = ScrollingPlot(dt=self.dt, npts=int(1.0/self.dt),
                                         parent=self, labels={'left': ('Membrane Potential', 'V'), 
                                                              'bottom': ('Time', 's')})
        self.vm_plot.setYRange(-90*mV, 50*mV)
        self.vm_plot.setXRange(-1000*ms, 0*ms)
        self.vm_plot.addLine(y=self.neuron.ek)
        self.vm_plot.addLine(y=self.neuron.ena)
        self.plot_splitter.addWidget(self.vm_plot)
        self.splitter.setSizes([300, 700])

        self.show()
        
        self.channel_params = [
            ChannelParameter(self.leak, 'Ileak'),
            ChannelParameter(self.hhna, 'INa'),
            ChannelParameter(self.hhk, 'IK'),
            ChannelParameter(self.dexh, 'IH'),
        ]
        for ch in self.channel_params:
            ch.plots_changed.connect(self.channel_plots_changed)
        self.channel_plots = {}
        
        self.params = pt.Parameter.create(name='params', type='group', children=[
            #dict(name='Preset', type='list', values=['', 'Basic Membrane', 'Hodgkin & Huxley']),
            dict(name='Run', type='bool', value=True),
            dict(name='Speed', type='float', value=0.3, limits=[0, 1], step=1, dec=True),
            dict(name='Temp', type='float', value=self.sim.temp, suffix='C', step=1.0),
            dict(name='Ion Channels', type='group', children=self.channel_params),
        ])
        self.ptree.setParameters(self.params)
        self.params.sigTreeStateChanged.connect(self.params_changed)
        
        self.runner = ndemo.SimRunner(self.sim)
        self.runner.add_request('soma.Vm', (self.neuron, 'Vm')) 
        self.runner.new_result.connect(mp.proxy(self.new_result, autoProxy=False))
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
            elif param is self.params.child('Temp'):
                self.sim.temp = val
            #elif param is self.params.child('Preset'):
                #self.load_preset(val)
        
    def channel_plots_changed(self, param, name, plot):
        key = param.name() + '.' + name
        if plot:
            # decide on y range, label , and units for new plot
            yranges = {
                'I': (-1*nA, 1*nA),
                'G': (0, 100*nS),
                'OP': (0, 1),
                'm': (0, 1),
                'h': (0, 1),
                'n': (0, 1),
            }
            color = {'I': 'c', 'G': 'y', 'OP': 'g'}.get(name, 0.7)
            units = {'I': 'A', 'G': 'S'}
            label = param.name() + ' ' + name
            if name in units:
                label = (label, units[name])
                
            # create new plot
            plt = ScrollingPlot(dt=self.vm_plot.dt, npts=self.vm_plot.npts,
                                labels={'left': label}, pen=color)
            plt.setXLink(self.vm_plot)
            plt.setYRange(*yranges.get(name, (0, 1)))
            
            # register this plot for later..
            self.channel_plots[key] = plt
            self.runner.add_request(key, (param.channel, name))
            
            # add new plot to splitter and resize all accordingly
            sizes = self.plot_splitter.sizes()
            self.plot_splitter.addWidget(plt)
            size = self.plot_splitter.height() / (len(sizes) + 1.)
            r = len(sizes) / (len(sizes)+1)
            sizes = [s * r for s in sizes] + [size]
            self.plot_splitter.setSizes(sizes)
        else:
            plt = self.channel_plots.pop(key)
            self.runner.remove_request(key)
            #self.plot_splitter.removeWidget(plt)
            plt.setParent(None)
            plt.close()
        
    def start(self):
        self.runner.start(blocksize=100)
        
    def stop(self):
        self.runner.stop()
        
    def new_result(self, result):
        self.last_result = result
        vm = result['soma.Vm'][1:]
        self.vm_plot.append(vm)
        
        for k, plt in self.channel_plots.items():
            if k not in result:
                continue
            plt.append(result[k][1:])

    def load_preset(self, preset):
        if preset == 'Basic Membrane':
            self.params['Temp'] = 6.3
            for p in self.params.child('Hodgkin & Huxley').children():
                p.setValue(False)
            for p in self.params.child('Lewis & Gerstner').children():
                p.setValue(False)
            self.params['Preset'] = ''


class ChannelParameter(pt.parameterTypes.SimpleParameter):
    
    # emitted when a plot should be shown or hidden
    plots_changed = QtCore.Signal(object, object, object)  # self, param, value
    
    def __init__(self, channel, name):
        self.channel = channel
        ch_params = [
            dict(name='ḡ', type='float', value=channel.gbar*cm**2, suffix='S/cm²', siPrefix=True, step=0.1, dec=True),
            dict(name='Erev', type='float', value=channel.erev, suffix='V', siPrefix=True, step=5*mV),
            dict(name='Plot I', type='bool', value=False),
            dict(name='Plot G', type='bool', value=False),
            dict(name='Plot OP', type='bool', value=False),
        ]
        for sv in channel.state_vars:
            ch_params.append(dict(name='Plot ' + sv, type='bool', value=False))
        
        pt.parameterTypes.SimpleParameter.__init__(self, name=name, type='bool', 
                                                   value=channel.enabled, children=ch_params,
                                                   expanded=False)
        self.sigTreeStateChanged.connect(self.treeChange)
        
    def treeChange(self, root, changes):
        for param, change, val in changes:
            if change != 'value':
                continue
            if param is self:
                self.channel.enabled = val
            elif param is self.child('ḡ'):
                self.channel.gbar = val/cm**2
            elif param is self.child('Erev'):
                self.channel.erev = val
            elif param.name().startswith('Plot'):
                self.plots_changed.emit(self, param.name()[5:], val)


class ScrollingPlot(pg.PlotWidget):
    def __init__(self, dt, npts, pen='w', **kwds):
        pg.PlotWidget.__init__(self, **kwds)
        self.data_curve = self.plot(pen=pen)
        self.data = np.array([], dtype=float)
        self.npts = npts
        self.dt = dt
        
    def append(self, data):
        self.data = np.append(self.data, data)
        if len(self.data) > self.npts:
            self.data = self.data[-self.npts:]
        t = np.arange(len(self.data)) * self.dt
        t -= t[-1]
        self.data_curve.setData(t, self.data)
        

        
if __name__ == '__main__':
    import sys
    win = DemoWindow()
    if sys.flags.interactive == 0:
        app.exec_()
