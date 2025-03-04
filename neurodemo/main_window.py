from __future__ import annotations
import typing
import sys
import numpy as np
import pyqtgraph as pg
import pyqtgraph.parametertree as pt

from . import qt
from . import ions
import neurodemo
import neurodemo.units as NU
from neurodemo.channelparam import ChannelParameter
from neurodemo.channelparam import IonConcentrations
from neurodemo.clampparam import ClampParameter
from neurodemo.neuronview import NeuronView

pg.setConfigOption('antialias', True)

# Disable obnoxious app nap on OSX 
# Many thanks to https://github.com/minrk/appnope
# if sys.platform == 'darwin':
    # v = [int(x) for x in platform.mac_ver()[0].split('.')]
    # if (v[0] == 10 and v[1] >= 9) or v[0] >= 11:
    #     import appnope
    #     appnope.nope()

app = pg.mkQApp()
app.setStyle("Fusion")  # necessary to remove double labels on mac os w/pyqtgraph until PR is done


class DemoWindow(qt.QWidget):
    def __init__(self, multiprocessing=False):
        if multiprocessing:
            # Enable running simulation in background process:
            import pyqtgraph.multiprocess as mp
            self.proc = mp.QtProcess(debug=False)
        else:
            self.proc = None

        self.app = pg.mkQApp()
        self.app.setStyle("fusion")
        self.app.setStyleSheet("QLabel{font-size: 11pt;} QText{font-size: 11pt;} {QWidget{font-size: 8pt;}")
        self.app.setStyleSheet("QTreeWidgetItem{font-size: 9pt;}") #  QText{font-size: 11pt;} {QWidget{font-size: 8pt;}")

        if self.proc is None:
            print(sys.platform, "running without mp")
            # do not use remote process:
            self.ndemo = neurodemo
        else:
            print(sys.platform, "running with mp")
            self.ndemo = self.proc._import('neurodemo')

        self.scrolling_plot_duration = 1.0 * NU.s
        self.result_buffer = ResultBuffer(max_duration=self.scrolling_plot_duration)

        self.dt = 20e-6 * NU.s
        self.integrator = 'solve_ivp'
        self.sim = self.ndemo.Sim(temp=6.3, dt=self.dt)
        if self.proc is not None:
            self.sim._setProxyOptions(deferGetattr=True)  # only if using remote process
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
        self.runner = self.ndemo.SimRunner(self.sim)
        self.runner.set_speed(0.5)

        # if using remote process (only on Windows):
        if self.proc is not None:
            self.runner.new_result.connect(mp.proxy(self.new_result, autoProxy=False, callSync='off'))
        else: # Darwin (macOS) and Linux:
            self.runner.new_result.connect(self.new_result) 

        # set up GUI
        qt.QWidget.__init__(self)

        self.fullscreen_widget = None
        self.resize(1400, 800)
        self.layout = qt.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
 
        self.splitter = qt.QSplitter(qt.Qt.Orientation.Horizontal)
        self.layout.addWidget(self.splitter, 0, 0,1, 1)
        self.ptree = pt.ParameterTree(showHeader=False)
        self.splitter.addWidget(self.ptree)
 
        self.ptree_stim = pt.ParameterTree(showHeader=False)
        self.splitter.addWidget(self.ptree_stim)
        
        self.plot_splitter = qt.QSplitter(qt.Qt.Orientation.Vertical)
        self.splitter.addWidget(self.plot_splitter)
        
        self.neuronview = NeuronView(self.neuron, mechanisms)
        self.plot_splitter.addWidget(self.neuronview)
        
        self.clamp_param = ClampParameter(self.clamp, self)
        self.ptree_stim.setParameters(self.clamp_param)
        self.clamp_param.plots_changed.connect(self.plots_changed)
        self.clamp_param.mode_changed.connect(self.mode_changed)

        self.channel_plots: typing.Dict[str, ScrollingPlot] = {}
        
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
            # These will use default Erev values defined in neuronsim.py, Section(SimObject) class, which
            # are INDEPENDENT of values defined in self.use_default_erev(). Consider replacing both with
            # calculated Erev.
            ch.plots_changed.connect(self.plots_changed)

        self.ion_concentrations = [IonConcentrations(ion) for ion in ions.all_ions]
        for ion in self.ion_concentrations:
            ion.updateTemperature(self.sim.temp)  # match temperature with an update

        self.params = pt.Parameter.create(name='Parameters', type='group', children=[
            dict(name='Preset', type='list', value='HH AP', limits=['Passive', 'HH AP', 'LG AP']),
            dict(name='Run/Stop', type='action', value=False),
            dict(name="Elapsed", type='float', value=0.0, suffix='s', siPrefix=True),
            dict(name="dt", type='float', value=20e-6, limits=[2e-6, 200e-6], suffix='s', siPrefix=True),
            dict(name="Method", type='list', value="solve_ivp", limits=['solve_ivp', 'odeint']),
            dict(name='Speed', type='float', value=self.runner.speed, limits=[0.001, 10], step=0.5, minStep=0.001, dec=True),
            dict(name="Plot Duration", type='float', value=1.0, limits=[0.1, 10], suffix='s', siPrefix=True, step=0.2),
            dict(name='Temp', type='float', value=self.sim.temp, limits=[0., 41.], suffix='C', step=1.0),
            dict(name='Capacitance', type='float', value=self.neuron.cap, limits=[0.1e-12, 1000.e-12], suffix='F', siPrefix=True, dec=True, children=[
                dict(name='Plot Current', type='bool', value=False),
            ]),
            dict(name='Ions', type='group', children=self.ion_concentrations),
            dict(name='Cell Schematic', type='bool', value=True, children=[
                dict(name='Show Circuit', type='bool', value=False),
            ]),
            # self.clamp_param,  # now in the adjacent window
            dict(name='Ion Channels', type='group', children=self.channel_params),
        ])

        self.use_calculated_erev()

        # Now that add_plot() sets x-axis limits, it must be called AFTER defining self.params,
        # rather than before, since it uses "Plot Duration" field in self.params.
        self.vm_plot = self.add_plot('soma.V', 'Membrane Potential', 'V')
        self.splitter.setSizes([300, 300, 800])

        self.ptree.setParameters(self.params)
        self.params.sigTreeStateChanged.connect(self.params_changed)
        # make Run/Stop button change color to indicate running state
        p = self.params.child("Run/Stop")
        rsbutton = list(p.items.keys())[0].button
        rsbutton.setCheckable(True)  # toggle
        rsbutton.setStyleSheet("QPushButton { background-color: #225522}"
                               "QPushButton:checked { background-color: #882222}" )

        # self.start()  # if autostart desired

        # Make Pulse Once and Pulse Sequence buttons green, since they are also capable of starting the simulator.
        # Unlike Run/Stop, these do not change color, since they are one-shot, and don't run continuously.
        for p in [self.clamp_param.child("Pulse", "Pulse Once"), self.clamp_param.child("Pulse", "Pulse Sequence")]:
            rsbutton = list(p.items.keys())[0].button
            rsbutton.setStyleSheet("QPushButton { background-color: #225522}")

        self.clamp_param['Plot Command'] = True

        # Set default heights of neuronview, membrane voltage, and other initial scrolling graphs.
        self.plot_splitter.setSizes([300, 500] + [200] * (len(self.plot_splitter.sizes()) - 2))

        # Any additional plots must be added after setting default heights
        self.pause_shortcut = qt.QShortcut(qt.QKeySequence('Space'), self)
        self.pause_shortcut.activated.connect(self.pause)
        self.slow_shortcut = qt.QShortcut(qt.QKeySequence('-'), self)
        self.slow_shortcut.activated.connect(self.slower)
        self.fast_shortcut = qt.QShortcut(qt.QKeySequence('='), self)
        self.fast_shortcut.activated.connect(self.faster)
        self.fullscreen_shortcut = qt.QShortcut(qt.QKeySequence('F11'), self)
        self.fullscreen_shortcut.activated.connect(self.fullscreen)

        #self.fullscreen_shortcut.setContext().Qt.ShortcutContext(qt.Qt.ApplicationShortcut)
        self.show()

    def params_changed(self, root, changes):
        for param, change, val in changes:
            path = self.params.childPath(param)
            if path[0] == "Run/Stop":
                if self.running() is True:
                    self.stop()
                else:
                    self.start()
            if change != 'value':
                continue

            elif param is self.params.child('Speed'):
                self.runner.set_speed(val)
            elif param is self.params.child('dt'):
                self.reset_dt(val)
            elif param is self.params.child("Method"):
                self.integrator = val
                self.sim.set_integrator(val)
            elif param is self.params.child('Plot Duration'):
                self.set_scrolling_plot_duration(val)
            elif param is self.params.child('Temp'):
                self.sim.temp = val
                # also update the ion channel values = specifically Erev
                for ion in self.ion_concentrations:
                    ion.updateTemperature(self.sim.temp)
            elif param is self.params.child('Capacitance'):
                self.neuron.cap = val
            elif param is self.params.child('Capacitance', 'Plot Current'):
                if val:
                    self.add_plot('soma.I', "Membrane Capactiance", 'I')
                else:
                    self.remove_plot('soma.I')
            elif param is self.params.child('Preset'):
                self.load_preset(val)
            elif param is self.params.child('Cell Schematic'):
                self.neuronview.setVisible(val)
            elif param is self.params.child('Cell Schematic', 'Show Circuit'):
                self.neuronview.show_circuit(val)
            elif param in [  # change in ion concentrations and erev...
                    self.params.child('Ions', "Na", "[C]in"),
                    self.params.child('Ions', "Na", "[C]out"),
                    self.params.child('Ions', "K", "[C]in"),
                    self.params.child('Ions', "K", "[C]out"),
                    self.params.child('Ions', "Cl", "[C]in"),
                    self.params.child('Ions', "Cl", "[C]out"),
                ]:
                self.use_calculated_erev()  # force update of erevs
            elif param in self.params.child('Ion Channels'):
                if not param.value():
                    # When turning ion off, also remove any plots associated with this ion
                    for p in param.children():
                        if p.name().startswith('Plot') and p.value():
                            # Turn off any plots that are currently ON
                            p.setValue(False)

    def plots_changed(self, param, channel, name, plot):
        key = channel.name + '.' + name
        if plot:
            self.add_plot(key, param.name(), name)
        else:
            self.remove_plot(key)

    def command_units(self):
        return 'V' if self.clamp_param['Mode'] == 'vc' else 'A'

    def mode_changed(self):
        key = self.clamp.name + '.cmd'
        if key in self.channel_plots:
            plt = self.channel_plots[key]
            plt.setLabels(left=(plt.label_name, self.command_units()))

    def add_plot(self, key, pname, name) -> ScrollingPlot:
        # decide on y range, label, and units for new plot
        yranges = {
            'V': (-100*NU.mV, 50*NU.mV),
            'cmd': None,
            'I': (-1*NU.nA, 1*NU.nA),
            'G': (0, 100*NU.nS),
            'OP': (0, 1),
            'm': (0, 1),
            'h': (0, 1),
            'n': (0, 1),
        }
        color = {'I': 'c', 'G': 'y', 'OP': 'g', 'V': 'w'}.get(name, 0.7)
        units = {'I': 'A', 'G': 'S', 'V': 'V', 'cmd': self.command_units()}
        label = pname + ' ' + name
        if name in units:
            label = (label, units[name])
            
        # create new scrolling plot
        plt = ScrollingPlot(dt=self.dt, npts=int(self.scrolling_plot_duration / self.dt),
                            labels={'left': label}, pen=color)
        plt.label_name = label[0]

        if hasattr(self, 'vm_plot'):
            plt.setXLink(self.vm_plot)
        else:
            plt.setXRange(-self.scrolling_plot_duration, 0)
        yrange = yranges.get(name, (0, 1))
        if yrange is None:
            plt.enableAutoRange(y=True)
        else:
            plt.setYRange(*yrange)

        # Add a vertical line
        plt.hover_line = pg.InfiniteLine(pos=0, angle=90, movable=False)
        plt.addItem(plt.hover_line, ignoreBounds=True)
        plt.hover_line.setVisible(False)

        # register this plot for later..
        self.channel_plots[key] = plt

        # Prevent user from zooming out beyond actual data time limits
        plt.setLimits(xMin=-self.params['Plot Duration'], xMax=0)
        
        # add new plot to splitter and resize all accordingly
        sizes = self.plot_splitter.sizes()
        self.plot_splitter.addWidget(plt)
        size = self.plot_splitter.height() / (len(sizes) + 1.)
        r = len(sizes) / (len(sizes)+1)
        sizes = [int(s * r) for s in sizes] + [int(size)]
        self.plot_splitter.setSizes(sizes)

        # Ask sequence plotter to update as well
        self.clamp_param.add_plot(key, label)

        # Track mouse over plot
        plt.plotItem.scene().sigMouseHover.connect(self.mouse_moved_over_plot)

        return plt
            
    def remove_plot(self, key):
        plt = self.channel_plots.pop(key)
        self.clamp_param.remove_plot(key)
        plt.plotItem.scene().sigMouseHover.disconnect(self.mouse_moved_over_plot)
        plt.setParent(None)
        plt.close()
        
    def mouse_moved_over_plot(self, items):
        # only process hover events while paused
        if self.running():
            return
        if len(items) == 0:
            # Can happen if mouse is dragged past valid time limit
            return
        item = items[0]
        widget = item.getViewWidget()
        globalPos = qt.QCursor.pos()
        localPos = widget.mapFromGlobal(globalPos)
        scenePos = item.mapFromDevice(localPos)
        if isinstance(item, pg.PlotItem):
            viewPos = item.vb.mapSceneToView(scenePos)
            self.set_hover_time(viewPos.x())

    def set_hover_time(self, t):
        """Move vertical lines to time *t* and update the schematic accordingly.
        """
        for plt in self.channel_plots.values():
            plt.hover_line.setVisible(True)
            plt.hover_line.setPos(t)
        state = self.result_buffer.get_state_at_time(t)
        if state is not None:
           self.clamp_param['Cursor values', 'Relative time'] = t
           if 'soma.V' in self.clamp_param.plot_keys:
               self.clamp_param['Cursor values', 'Memb voltage'] = state['soma.V']
           if 'soma.PatchClamp.cmd' in self.channel_plots.keys():
               # Command values are not in the state variable, so we grab them out of the actual plot object
               plt: ScrollingPlot = self.channel_plots['soma.PatchClamp.cmd']
               dc: pg.PlotDataItem = plt.data_curve
               [x, y] = [dc.xData, dc.yData]
               dx = x[1] - x[0]
               idx = int(t / dx)  # Note that t will be negative, since it represents time before present. This will make idx negative, so index will count from end
               if idx >= -len(y):
                   self.clamp_param['Cursor values', 'Command'] = y[idx]
           self.neuronview.update_state(state)

    def running(self):
        return self.runner.running()

    def start(self, stop_after_cmd=False):

        self.runner.start(stop_after_cmd=stop_after_cmd, blocksize=1000)
        for plt in self.channel_plots.values():
            plt.hover_line.setVisible(False)
        
    def stop(self):
        self.runner.stop()

    def reset_dt(self, val):
        was_running = self.running()
        if was_running:
            self.stop()
        self.dt = val

        self.clamp_param.set_dt(self.dt)
        self.sim.change_dt(self.dt)
        if was_running:
            self.start() # restart.
        # make sure to change dt elsewhere as well.
        self.set_scrolling_plot_dt(self.dt)

    def set_scrolling_plot_dt(self, val):
        was_running = self.running()
        if was_running:
            self.stop()
        for k in self.channel_plots.keys():
            self.channel_plots[k].set_dt(val)
        if was_running:
            self.start() # restart.

    def set_scrolling_plot_duration(self, val):
        self.scrolling_plot_duration = val
        # Update x-axis limits for membrane voltage plot, so we don't "lose" the traces by accident.
        self.vm_plot.setLimits(xMin=-val)
        for k in self.channel_plots.keys():
            self.channel_plots[k].set_duration(val)
            # Update range and x-limits on all other ScrollingPlot objects
            # Note: ScrollingPLot.setLimits() and .setXRange() are wrapped from pyqtgraph ViewBox
            self.channel_plots[k].setLimits(xMin=-val, xMax=0)
            self.channel_plots[k].setXRange(-val, 0)

    def pause(self):
        self.params['Run'] = not self.params['Run']

    def slower(self):
        self.params['Speed'] *= 0.5
        
    def faster(self):
        self.params['Speed'] = min(self.params['Speed'] * 2.0, 10)

    def fullscreen(self):
        if self.fullscreen_widget is None:
            w = qt.QApplication.focusWidget()
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
        
    def new_result(self, result):
        for k, plt in self.channel_plots.items():
            if k not in result:
                continue
            if isinstance(result[k], float):
                plt.append(result[k])
            else:
                plt.append(result[k][1:])
            
        # Let the clamp decide which triggered regions of the data to extract
        # for pulse plots
        self.clamp_param.new_result(result)

        self.params['Elapsed'] = result['t'][-1]

        if self.runner.stop_after_cmd:
            # We have elected to stop after command (either single pulse or sequence) is done
            if len(self.clamp_param.triggers) == 0 and self.runner.sim.cmd_done():
                # Stop after command queue and trigger queue are BOTH empty
                self.runner.stop()
        
        # update the schematic
        self.neuronview.update_state(result.get_final_state())

        # store a running buffer of results
        self.result_buffer.add(result)

    def _get_Eh(self):
        ENa = self.params.child('Ions', 'Na')
        ena = ENa.param('Erev').value()
        Ek = self.params.child('Ions', 'K')
        ek = Ek.param('Erev').value()
        x = 0.75  # this sets ration of gk to gna given ions
        eh = x*ek + (1.0-x)*ena  # gives -43 mV when ena is 50 and ek is -74
        return eh

    def _get_Eleak(self):
        ENa = self.params.child('Ions', 'Na')
        ena = ENa.param('Erev').value()
        Ek = self.params.child('Ions', 'K')
        ek = Ek.param('Erev').value()
        x = 0.84677419  # gives -55 mV when ena is 50 and ek is -74
        eleak = x*ek + (1.0-x)*ena  
        return eleak

    def use_calculated_erev(self):
        chans = self.params.child('Ion Channels')
        ENa = self.params.child('Ions', 'Na')
        ENa_erev = ENa.param('Erev').value()
        for ch in ["INa", "INa1"]:
            chans[f"soma.{ch:s}", 'Erev'] = ENa_erev
        Ek = self.params.child('Ions', 'K')
        EK_erev = Ek.param('Erev').value()
        for ch in ["IK", "IKf", "IKs"]:
            chans[f"soma.{ch:s}", 'Erev'] = EK_erev
        ECl = self.params.child('Ions', 'Cl')
        ECl_erev = ECl.param('Erev').value()
        # for ch in ["ICl"]:
        #     chans[f"soma.{ch:s}", 'Erev'] = ECl.param('Erev').value()
        Eleak_erev = self._get_Eleak()
        for ch in ["Ileak"]:
            chans[f"soma.{ch:s}", 'Erev'] = Eleak_erev
        Eh_erev = self._get_Eh()
        for ch in ["IH"]:
            chans[f"soma.{ch:s}", 'Erev'] = Eh_erev
        self.set_hh_erev(ENa_erev, EK_erev, Eleak_erev, Eh_erev)
        self.set_lg_erev(ENa_erev, EK_erev, EK_erev, -55*NU.mV)
    
    def set_hh_erev(self, ENa_erev=50*NU.mV, EK_erev=-74*NU.mV,
                    Eleak_erev=-55*NU.mV, Eh_erev=-43*NU.mV):
        """Set new reversal potentials for the HH currents

        Args:
            ENa_erev (float): new value for Na channel
            EK_erev (float): new value for K channel
            Eleak_erev (float): new value for leak channel
            Eh_erev (float): new value for dexh (IH)
        """
        self.hhna.set_erev(ENa_erev)
        self.hhk.set_erev(EK_erev)
        self.dexh.set_erev(Eh_erev)
        self.leak.set_erev(Eleak_erev)

    def set_lg_erev(self, ENa_erev=74*NU.mV, EKf_erev=-90*NU.mV,
                EKs_erev=-90*NU.mV, Eleak_erev=-70*NU.mV):
        """Set new reversal potentials for the LG currents

        Args:
            ENa_erev (float): new value for Na channel
            EKf_erev (float): new value for Kf channel
            EKs_erev (float): new value for Ks channel
            Eleak_erev (float): new value for leak channel
        """
        self.lgna.set_erev(ENa_erev)
        self.lgkf.set_erev(EKf_erev)
        self.lgks.set_erev(EKs_erev)
        self.leak.set_erev(Eleak_erev)


    def load_preset(self, preset):
        """Load preset configurations for the simulations.

        Args:
            preset (string): which preset values to select and load

        Raises:
            ValueError: if preset is not known.
        """
        if preset == 'Passive':
            self.params['Temp'] = 6.3
            self.params['Speed'] = 0.5
            self.clamp_param['Plot Current'] = False
            self.clamp_param['Plot Voltage'] = False
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
            self.set_ions_off()
            self.neuron.set_default_erev()

        elif preset == 'HH AP':
            self.params['Temp'] = 6.3
            self.params['Speed'] = 0.5
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
            self.set_ions_off()
            self.set_hh_erev()

        elif preset == 'LG AP':
            self.params['Temp'] = 37
            self.params['Speed'] = 0.5
            chans = self.params.child('Ion Channels')
            chans['soma.Ileak'] = True
            chans['soma.Ileak', 'Erev'] = -70*NU.mV
            chans['soma.Ileak', 'Gmax'] = 2.5*NU.nS
            chans['soma.INa'] = False
            chans['soma.IK'] = False
            chans['soma.IH'] = False
            chans['soma.INa1'] = True
            chans['soma.INa1', "Erev"] = 74 * NU.mV
            chans['soma.IKf'] = True
            chans['soma.IKf', "Erev"] = -90 * NU.mV
            chans['soma.IKs'] = True
            chans['soma.IKs', "Erev"] = -90 * NU.mV
            self.set_ions_off()
            self.set_lg_erev()
        else:
            raise ValueError("Preset is not one of the implemented values")
            
        self.params['Preset'] = preset

    def closeEvent(self, ev):
        self.runner.stop()
        # self.proc.close()
        qt.QApplication.instance().quit()


class ScrollingPlot(pg.PlotWidget):
    def __init__(self, dt, npts, pen='w', **kwds):
        pg.PlotWidget.__init__(self, **kwds)
        self.showGrid(True, True)
        self.data_curve = self.plot(pen=pen)
        self.data = np.array([], dtype=float)
        self.npts = npts
        self.dt = dt
        self.plot_duration = 1.0
    
    def set_dt(self, dt):
        self.dt = dt
        # update npts as well
        self.npts=int(self.plot_duration / self.dt)
        # print(self.plot_duration, self.npts, self.dt)

    def set_duration(self, dur):
        self.plot_duration = dur
        self.npts=int(self.plot_duration / self.dt)
        self.setXRange(-self.plot_duration, 0)
        # print(self.plot_duration, self.npts, self.dt, len(self.data))

    def append(self, data):
        # print("len data, len self.data: ", len(data), len(self.data))
        self.data = np.concatenate((self.data, data), axis=0)
        if len(self.data) >= self.npts:
            self.data = self.data[-self.npts:]
        t = np.arange(len(self.data)) * self.dt
        t -= t[-1]
        # print("appending npts: ", len(self.data), self.npts, self.dt, self.plot_duration)
        self.data_curve.setData(t, self.data)


class ResultBuffer:
    def __init__(self, max_duration=10):
        self.max_duration = max_duration
        self.results = []

    def add(self, result):
        self.results.append(result)

    def get_state_at_time(self, t):
        if len(self.results) == 0:
            return None
        if t < 0:
            t = self.results[-1]['t'][-1] + t
        for result in self.results:
            if result['t'][0] <= t <= result['t'][-1]:
                return result.get_state_at_time(t)
        return None
