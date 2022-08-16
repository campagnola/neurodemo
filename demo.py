# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""

# make sure we get the right pyqtgraph.
from dataclasses import dataclass
import sys
import platform
# import appnope
import numpy as np
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp
import pyqtgraph.parametertree as pt
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
from pyqtgraph.debug import ThreadTrace

import neurodemo
import neurodemo.units as NU
from neurodemo.channelparam import ChannelParameter
from neurodemo.channelparam import IonConcentrations
from neurodemo.clampparam import ClampParameter
from neurodemo.neuronview import NeuronView

pg.setConfigOption('antialias', True)

# Disable obnoxious app nap on OSX 
# Many thanks to https://github.com/minrk/appnope
app = pg.mkQApp()
if sys.platform == 'darwin':
    # v = [int(x) for x in platform.mac_ver()[0].split('.')]
    # if (v[0] == 10 and v[1] >= 9) or v[0] >= 11:
    #     import appnope
    #     appnope.nope()
    app.setStyle("Fusion")  # necessary to remove double labels on mac os w/pyqtgraph until PR is done

@dataclass
class IonClass:
    name: str=""
    Cout: float = 1.0
    Cin: float = 1.0
    valence: float = 1.0
    Erev: float = 0.0
    enabled: bool=False

class DemoWindow(QtWidgets.QWidget):
    def __init__(self, proc):

        self.dt = 25*NU.us
        self.proc = proc
        # set up simulation in remote process
        # self.ndemo = self.proc._import('neurodemo')
        # or do not use remote process:
        self.ndemo = neurodemo
        self.sim = self.ndemo.Sim(temp=6.3, dt=self.dt)
        # self.sim._setProxyOptions(deferGetattr=True)  # only if using remote process
        self.neuron = self.ndemo.Section(name='soma')
        self.sim.add(self.neuron)
        
        self.hhna = self.neuron.add(self.ndemo.HHNa())
        self.leak = self.neuron.add(self.ndemo.Leak())
        self.hhk = self.neuron.add(self.ndemo.HHK())
        self.dexh = self.neuron.add(self.ndemo.IH())
        self.dexh.enabled = False
        self.lgna = self.neuron.add(self.ndemo.LGNa())
        self.lgkf = self.neuron.add(self.ndemo.LGKfast())
        self.lgks = self.neuron.add(self.ndemo.LGKslow())
        self.lgna.enabled = False
        self.lgkf.enabled = False
        self.lgks.enabled = False

    
        
        self.clamp = self.neuron.add(self.ndemo.PatchClamp(mode='ic'))
        
        mechanisms = [self.clamp, self.hhna, self.leak, self.hhk, self.dexh, self.lgna, self.lgkf, self.lgks]

        # loop to run the simulation indefinitely
        self.running = False
        self.runner = self.ndemo.SimRunner(self.sim)
        self.runner.add_request('t') 
        # if using remote process:
        # self.runner.new_result.connect(mp.proxy(self.new_result, autoProxy=False, callSync='off'))
        self.runner.new_result.connect(self.new_result,) 

        # set up GUI
        QtGui.QWidget.__init__(self)
        self.fullscreen_widget = None
        self.resize(1024, 768)
        self.layout = QtWidgets.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        self.layout.addWidget(self.splitter, 0, 0)
        
        self.ptree = pt.ParameterTree(showHeader=False)
        self.splitter.addWidget(self.ptree)
        
        self.plot_splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        self.splitter.addWidget(self.plot_splitter)
        
        self.neuronview = NeuronView(self.neuron, mechanisms)
        self.plot_splitter.addWidget(self.neuronview)
        
        self.channel_plots = {}
        
        self.channel_params = [
            ChannelParameter(self.leak),
            ChannelParameter(self.hhna),
            ChannelParameter(self.hhk),
            ChannelParameter(self.dexh),
            ChannelParameter(self.lgna),
            ChannelParameter(self.lgkf),
            ChannelParameter(self.lgks),
        ]

        for ch in self.channel_params:
            ch.plots_changed.connect(self.plots_changed)

        self.ion_concentrations = [
            IonConcentrations(IonClass(name='Na', Cout=140.0, Cin=8.0, valence=+1, enabled=False)),
            IonConcentrations(IonClass(name='K', Cout=4., Cin=140., valence=+1, enabled=False)),
        ]
        for ion in self.ion_concentrations:
            ion.updateErev(self.sim.temp)
        
        self.clamp_param = ClampParameter(self.clamp, self)
        self.clamp_param.plots_changed.connect(self.plots_changed)

        self.vm_plot = self.add_plot('soma.V', 'Membrane Potential', 'V')
        
        self.splitter.setSizes([300, 650])

        self.params = pt.Parameter.create(name='Parameters', type='group', children=[
            dict(name='Preset', type='list', values=['HH Action Potential', 'Passive Membrane', 'LG Action Potential']),
            dict(name='Run/Stop', type='action', value=False),
            dict(name='Speed', type='float', value=1.0, limits=[0, 10], step=1, dec=True),
            dict(name='Temp', type='float', value=self.sim.temp, suffix='C', step=1.0),
            dict(name='Capacitance', type='float', value=self.neuron.cap, suffix='F', siPrefix=True, dec=True),
            dict(name='Ions', type='group', children=self.ion_concentrations),            
            dict(name='Cell Schematic', type='bool', value=True, children=[
                dict(name='Show Circuit', type='bool', value=False),
            ]),
            self.clamp_param,
            dict(name='Ion Channels', type='group', children=self.channel_params),
        ])
        self.ptree.setParameters(self.params)
        self.params.sigTreeStateChanged.connect(self.params_changed)
        # make Run/Stop button change color to indicate running state
        p = self.params.child("Run/Stop")
        rsbutton = list(p.items.keys())[0].button
        rsbutton.setCheckable(True)  # toggle
        rsbutton.setStyleSheet("QPushButton { background-color: #225522}"
                               "QPushButton:checked { background-color: #882222}" )

        # self.start()  # if autostart desired

        self.clamp_param['Plot Current'] = True
        self.plot_splitter.setSizes([300, 500, 200])
        
        self.pause_shortcut = QtGui.QShortcut(QtGui.QKeySequence('Space'), self)
        self.pause_shortcut.activated.connect(self.pause)
        self.slow_shortcut = QtGui.QShortcut(QtGui.QKeySequence('-'), self)
        self.slow_shortcut.activated.connect(self.slower)
        self.fast_shortcut = QtGui.QShortcut(QtGui.QKeySequence('='), self)
        self.fast_shortcut.activated.connect(self.faster)
        self.fullscreen_shortcut = QtGui.QShortcut(QtGui.QKeySequence('F11'), self)
        self.fullscreen_shortcut.activated.connect(self.fullscreen)

        #self.fullscreen_shortcut.setContext().Qt.ShortcutContext(QtCore.Qt.ApplicationShortcut)
        self.show()


    def params_changed(self, root, changes):
        for param, change, val in changes:
            path = self.params.childPath(param)
            if path[0] == "Run/Stop":
                if self.running is True:
                    self.stop()
                    self.running = False
                else:
                    self.start()
                    self.running = True
            if change != 'value':
                continue

            elif param is self.params.child('Speed'):
                self.runner.set_speed(val)
            elif param is self.params.child('Temp'):
                self.sim.temp = val
                # also update the ion channel values = specifically Erev
                for ion in self.ion_concentrations:
                    ion.updateErev(self.sim.temp)
            elif param is self.params.child('Capacitance'):
                self.neuron.cap = val
            elif param is self.params.child('Preset'):
                self.load_preset(val)
            elif param is self.params.child('Cell Schematic'):
                self.neuronview.setVisible(val)
            elif param is self.params.child('Cell Schematic', 'Show Circuit'):
                self.neuronview.show_circuit(val)
            elif param is self.params.child('Ions', 'Na'):
                if val:
                    self.use_calculated_erev()
                else:
                    self.use_default_erev()
        
    def plots_changed(self, param, channel, name, plot):
        key = channel.name + '.' + name
        if plot:
            self.add_plot(key, param.name(), name)
        else:
            self.remove_plot(key)

    def add_plot(self, key, pname, name):
        # decide on y range, label, and units for new plot
        yranges = {
            'V': (-100*NU.mV, 50*NU.mV),
            'I': (-1*NU.nA, 1*NU.nA),
            'G': (0, 100*NU.nS),
            'OP': (0, 1),
            'm': (0, 1),
            'h': (0, 1),
            'n': (0, 1),
        }
        color = {'I': 'c', 'G': 'y', 'OP': 'g', 'V': 'w'}.get(name, 0.7)
        units = {'I': 'A', 'G': 'S', 'V': 'V'}
        label = pname + ' ' + name
        if name in units:
            label = (label, units[name])
            
        # create new plot
        plt = ScrollingPlot(dt=self.dt, npts=int(1.0 / self.dt),
                            labels={'left': label}, pen=color)
        if hasattr(self, 'vm_plot'):
            plt.setXLink(self.vm_plot)
        else:
            plt.setXRange(-1, 0)
        plt.setYRange(*yranges.get(name, (0, 1)))
        
        # register this plot for later..
        self.channel_plots[key] = plt
        self.runner.add_request(key)
        
        # add new plot to splitter and resize all accordingly
        sizes = self.plot_splitter.sizes()
        self.plot_splitter.addWidget(plt)
        size = self.plot_splitter.height() / (len(sizes) + 1.)
        r = len(sizes) / (len(sizes)+1)
        sizes = [int(s * r) for s in sizes] + [int(size)]
        self.plot_splitter.setSizes(sizes)
        
        # Ask sequence plotter to update as well
        self.clamp_param.add_plot(key, label)
        return plt
            
    def remove_plot(self, key):
        plt = self.channel_plots.pop(key)
        self.runner.remove_request(key)
        self.clamp_param.remove_plot(key)
        plt.setParent(None)
        plt.close()
        
    def start(self):
        self.runner.start(blocksize=500)
        # set button color
        
    def stop(self):
        self.runner.stop()
        # reset button color
        
    def pause(self):
        self.params['Run'] = not self.params['Run']

    def slower(self):
        self.params['Speed'] *= 0.5
        
    def faster(self):
        self.params['Speed'] = min(self.params['Speed'] * 2.0, 10)

    def fullscreen(self):
        if self.fullscreen_widget is None:
            w = QtGui.QApplication.focusWidget()
            ind = self.plot_splitter.indexOf(w)
            if ind < 0:
                return
            self.fs_widget_index = ind
            w.setParent(None)
            w.showFullScreen()
            self.fullscreen_widget = w
        else:
            self.fullscreen_widget.showNormal()
            self.plot_splitter.insertWidget(self.fs_widget_index, self.fullscreen_widget)
            self.fullscreen_widget = None
        
    def new_result(self, final_state, result):
        for k, plt in self.channel_plots.items():
            if k not in result:
                continue
            plt.append(result[k][1:])
            
        # Let the clamp decide which triggered regions of the data to extract
        # for pulse plots
        self.clamp_param.new_result(result)
        
        # update the schematic
        self.neuronview.update_state(final_state)

    def use_calculated_erev(self):
        print("Using caclculated ERevs")
        chans = self.params.child('Ion Channels')
        ENa = self.params.child('Ions', 'Na')
        for ch in ["INa", "INa1"]:
            chans[f"soma.{ch:s}", 'Erev'] = ENa.param('Erev').value()
        Ek = self.params.child('Ions', 'K')
        for ch in ["IK", "IKf", "IKs"]:
            chans[f"soma.{ch:s}", 'Erev'] = Ek.param('Erev').value()
    
    def use_default_erev(self):
        print("Using default Erevs")
        chans = self.params.child('Ion Channels')
        ENa_revs = {"INa": 50, "INa1": 74}
        for ch in ["INa", "INa1"]:
            chans[f"soma.{ch:s}", 'Erev'] = ENa_revs[ch]
        EK_revs = {"IK": -74, "IKf": -90, "IKs": -90}
        for ch in ["IK", "IKf", "IKs"]:
            chans[f"soma.{ch:s}", 'Erev'] = EK_revs[ch]


    def load_preset(self, preset):
        if preset == 'Passive Membrane':
            self.params['Temp'] = 6.3
            self.params['Speed'] = 1.0
            self.params['Patch Clamp', 'Plot Current'] = False
            self.params['Patch Clamp', 'Plot Voltage'] = False
            chans = self.params.child('Ion Channels')
            chans['soma.Ileak'] = True
            chans['soma.Ileak', 'Erev'] = 0
            chans['soma.Ileak', "Gmax"] = 1*NU.nS
            chans['soma.INa'] = False
            chans['soma.IK'] = False
            chans['soma.IH'] = False
            chans['soma.INa1'] = False
            chans['soma.IKf'] = False
            chans['soma.IKs'] = False

        elif preset == 'HH Action Potential':
            self.params['Temp'] = 6.3
            self.params['Speed'] = 1.0
            chans = self.params.child('Ion Channels')
            chans['soma.Ileak'] = True
            chans['soma.Ileak', 'Erev'] = -55*NU.mV
            chans['soma.Ileak', "Gmax"] = 1*NU.nS
            chans['soma.INa'] = True
            chans['soma.IK'] = True
            chans['soma.IH'] = False
            chans['soma.INa1'] = False
            chans['soma.IKf'] = False
            chans['soma.IKs'] = False

        elif preset == 'LG Action Potential':
            self.params['Temp'] = 37
            self.params['Speed'] = 1.0
            chans = self.params.child('Ion Channels')
            chans['soma.Ileak'] = True
            chans['soma.Ileak', 'Erev'] = -70*NU.mV
            chans['soma.Ileak', 'Gmax'] = 2.5*NU.nS
            chans['soma.INa'] = False
            chans['soma.IK'] = False
            chans['soma.IH'] = False
            chans['soma.INa1'] = True
            chans['soma.IKf'] = True
            chans['soma.IKs'] = True
        else:
            raise ValueError("Preset is not one of the implemented values")
            
        self.params['Preset'] = preset

    def closeEvent(self, ev):
        self.runner.stop()
        self.proc.close()
        QtWidgets.QApplication.instance().quit()


class ScrollingPlot(pg.PlotWidget):
    def __init__(self, dt, npts, pen='w', **kwds):
        pg.PlotWidget.__init__(self, **kwds)
        self.showGrid(True, True)
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

# def main():


if __name__ == '__main__':
    import sys
    proc = mp.QtProcess(debug=False)
    win = DemoWindow(proc)
    if sys.flags.interactive == 0:
        app.exec()
        app.processEvents()
    # main()
