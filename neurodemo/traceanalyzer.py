# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""

import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph.parametertree as pt
from lmfit import Model
from lmfit.models import ExponentialModel

class TraceAnalyzer(QtGui.QWidget):
    def __init__(self, seq_plotter):
        QtGui.QWidget.__init__(self)
        self.plotter = seq_plotter
        
        self.layout = QtGui.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        self.hsplitter = QtGui.QSplitter(QtCore.Qt.Orientation.Horizontal)
        self.layout.addWidget(self.hsplitter)
        
        self.ptree = pt.ParameterTree(showHeader=False)
        self.hsplitter.addWidget(self.ptree)
        
        self.table = pg.TableWidget()
        self.hsplitter.addWidget(self.table)
        self.table.verticalHeader().hide()
        
        self.analysis_plot = EvalPlotter()
        self.hsplitter.addWidget(self.analysis_plot)
        
        self.hsplitter.setSizes([300, 200, 400])
        
        self.clear()
        
        self.params = TraceAnalyzerGroup(name="Analyzers")
        self.params.need_update.connect(self.update_analyzers)
        self.ptree.setParameters(self.params)

    def clear(self):
        self.data = []
        self.table.clear()
        
    def add_data(self, t, data, info):
        self.params.set_inputs(data.dtype.names)
        self.data.append((t, data, info))
        self.update_analysis()

    def update_analyzers(self):
        for anal in self.params.children():
            self.plotter.plots[anal['Input']].addItem(anal.rgn, ignoreBounds=True)
        self.update_analysis()
        
    def update_analysis(self):
        fields = ['cmd'] + [anal.name() for anal in self.params.children()]
        data = np.empty(len(self.data), dtype=[(str(f), float) for f in fields])
        for i, rec in enumerate(self.data):
            t, d, info = rec
            data['cmd'][i] = info['amp']
            for analysis in self.params.children():
                data[analysis.name()][i] = analysis.process(t, d)
        self.table.setData(data)
        self.analysis_plot.update_data(data)
        

class TraceAnalyzerGroup(pt.parameterTypes.GroupParameter):
    need_update = QtCore.Signal()

    def __init__(self, **kwds):
        analyses = ['min', 'max', 'mean', 'expTauDecay',  'spikeCount', 'spike_latency']
        self.inputs = []
        pt.parameterTypes.GroupParameter.__init__(self, addText='Add analysis..', addList=analyses, **kwds)

    def addNew(self, typ):
        param = TraceAnalyzerParameter(name=typ, analysis_type=typ, inputs=self.inputs, autoIncrementName=True)
        self.addChild(param)
        param.need_update.connect(self.need_update)
        self.need_update.emit()
        
    def set_inputs(self, inputs):
        self.inputs = list(inputs)
        self.inputs.remove('t')
        for ch in self.children():
            ch.set_input_list(inputs)


class TraceAnalyzerParameter(pt.parameterTypes.GroupParameter):
    need_update = QtCore.Signal(object)  # self

    def __init__(self, **kwds):
        kwds.update({'removable': True, 'renamable': False})
        childs = [
            dict(name='Input', type='list', values=kwds.pop('inputs')),
            dict(name='Type', type='list', value=kwds.pop('analysis_type'),
                values=['mean', 'min', 'max', 'expTauDecay', 'spikeCount', 'spikeLatency']),
            dict(name='Start', type='float', value=0, suffix='s', siPrefix=True, step=5e-3),
            dict(name='End', type='float', value=10e-3, suffix='s', siPrefix=True, step=5e-3),
            dict(name='Threshold', type='float', value=-30e-3, suffix='V', siPrefix=True, step=5e-3, visible=False),
        ]
        kwds['children'] = childs + kwds.get('children', [])
        
        pt.parameterTypes.GroupParameter.__init__(self, **kwds)
        self.sigTreeStateChanged.connect(self.tree_changed)
        
        self.rgn = pg.LinearRegionItem([self['Start'], self['End']])
        self.rgn.sigRegionChanged.connect(self.region_changed)
    
    def tree_changed(self, root, changes):
        for param, change, val in changes:
            if change != 'value':
                continue
            if param is self.child('Start') or param is self.child('End'):
                self.rgn.sigRegionChanged.disconnect(self.region_changed)
                try:
                    self.rgn.setRegion([self['Start'], self['End']])
                finally:
                    self.rgn.sigRegionChanged.connect(self.region_changed)
            elif param is self.child('Type'):
                needs_threshold = val in ['spikeCount', 'spike_latency']
                self.child('Threshold').setOpts(visible=needs_threshold)
            self.need_update.emit(self)
        
    def region_changed(self):
        start, end = self.rgn.getRegion()
        self.sigTreeStateChanged.disconnect(self.tree_changed)
        try:
            self['Start'] = start
            self['End'] = end
        finally:
            self.sigTreeStateChanged.connect(self.tree_changed)
            
        self.need_update.emit(self)
            
    def set_input_list(self, inputs):
        self.child('Input').setLimits(inputs)

    def process(self, t, data):
        dt = t[1] - t[0]
        i1 = int(self['Start'] / dt)
        i2 = int(self['End'] / dt)
        data = data[self['Input']][i1:i2]
        t = t[i1:i2]
        typ = self['Type']
        if typ == 'mean':
            return data.mean()
        elif typ == 'min':
            return data.min()
        elif typ == 'max':
            return data.max()
        elif typ.startswith('spike'):
            spikes = np.argwhere((data[1:] > self['Threshold']) & (data[:-1] < self['Threshold']))[:,0]
            if typ == 'spikeCount':
                return len(spikes)
            elif typ == 'spike_latency':
                if len(spikes) == 0:
                    return np.nan
                else:
                    return spikes[0] * dt
        elif typ == 'expTauDecay':
            return self.measure_tauDecay(data, t)
        elif typ == 'expTauRise4':
            return(self.measure_tauRise4(data, t))
            
    def measure_tau_old(self, data, t):
        from scipy.optimize import curve_fit
        dt = t[1] - t[0]
        def expfn(x, yoff, amp, tau):
            return yoff + amp * np.exp(-x / tau)
        guess = (data[-1], data[0] - data[-1], t[-1] - t[0])
        fit = curve_fit(expfn, t-t[0], data, guess)
        return fit[0][2]

    def measure_tauDecay(self, data, t):
        model = ExponentialModel()
        pars = model.guess(data-data[-1], x=t-t[0])
        result = model.fit(data-data[-1], pars,  x=t-t[0])
        return result.params['decay'] # fit[0][2]            
        
    def measure_tauRise4(self, data, t):
        # this is not working quite right yet... 
        print('taurise4')
        def expfn4(x, amp, tau):
            return amp * ((1.0-np.exp(-x / tau))**4.0)
        emodel = Model(expfn4)
        params = emodel.make_params()
        d = data-data[0]
        tp = t-t[0]
        emodel.set_param_hint('tau', value=t[-1]-t[0], min=0.0001)
        emodel.set_param_hint('amp', value=data[-1]-data[0])
        result = emodel.fit(d[1:], params, x=tp[1:])
        print(result.params)
        return result.params['tau'] # fit[0][2]       

class EvalPlotter(QtGui.QWidget):
    def __init__(self):
        self.held_plots = []
        self.last_curve = None
        self.held_index = 0
        
        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        
        self.x_label = QtGui.QLabel('X data')
        self.y_label = QtGui.QLabel('Y data')
        self.layout.addWidget(self.x_label, 0, 0)
        self.layout.addWidget(self.y_label, 1, 0)
        
        self.x_code = QtGui.QLineEdit('cmd')
        self.y_code = QtGui.QLineEdit()
        self.layout.addWidget(self.x_code, 0, 1)
        self.layout.addWidget(self.y_code, 1, 1)
        
        self.x_units_label = QtGui.QLabel('units')
        self.y_units_label = QtGui.QLabel('units')
        self.layout.addWidget(self.x_units_label, 0, 2)
        self.layout.addWidget(self.y_units_label, 1, 2)
        
        self.x_units_text = QtGui.QLineEdit('A')
        self.y_units_text = QtGui.QLineEdit()
        self.layout.addWidget(self.x_units_text, 0, 3)
        self.layout.addWidget(self.y_units_text, 1, 3)

        self.plot = pg.PlotWidget()
        self.layout.addWidget(self.plot, 2, 0, 1, 4)

        self.hold_plot_btn = QtGui.QPushButton('Hold Plot')
        self.layout.addWidget(self.hold_plot_btn, 3, 0, 1, 2)
        self.clear_plot_btn = QtGui.QPushButton('Clear Plot')
        self.layout.addWidget(self.clear_plot_btn, 3, 2, 1, 2)
        
        self.x_code.editingFinished.connect(self.replot)
        self.y_code.editingFinished.connect(self.replot)
        self.x_units_text.editingFinished.connect(self.replot)
        self.y_units_text.editingFinished.connect(self.replot)
        self.hold_plot_btn.clicked.connect(self.hold_plot)
        self.clear_plot_btn.clicked.connect(self.clear_plot)
        
    def update_data(self, data):
        self.data = data
        self.replot()
        
    def replot(self):
        data = self.data
        ns = {}
        for k in data.dtype.names:
            ns[k.replace(' ', '_')] = data[k]
        xcode = str(self.x_code.text())
        ycode = str(self.y_code.text())
        if xcode == '' or ycode == '':
            return
        
        try:
            x = eval(xcode, ns)
        except:
            pg.debug.printExc('Error evaluating plot x values:')
            self.x_code.setStyleSheet("QLineEdit { border: 2px solid red; }")
            return
        else:
            self.x_code.setStyleSheet("")
            
        try:
            y = eval(ycode, ns)
        except:
            pg.debug.printExc('Error evaluating plot y values:')
            self.y_code.setStyleSheet("QLineEdit { border: 2px solid red; }")
            return
        else:
            self.y_code.setStyleSheet("")
            
        if self.last_curve is None:
            self.last_curve = self.plot.plot(x, y, symbol='o', symbolBrush=(self.held_index, 10))
        else:
            self.last_curve.setData(x, y)
            
        self.plot.setLabels(bottom=(xcode, self.x_units_text.text()),
                            left=(ycode, self.y_units_text.text()))
        
    def hold_plot(self):
        if self.last_curve is None:
            return
        self.held_plots.append(self.last_curve)
        self.held_index += 1
        self.last_curve = None

    def clear_plot(self):
        self.held_plots = []
        self.held_index = 0
        self.plot.clear()
        self.replot()
