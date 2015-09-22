# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""
from __future__ import division, unicode_literals
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph.parametertree as pt


class TraceAnalyzer(QtGui.QWidget):
    def __init__(self, seq_plotter):
        QtGui.QWidget.__init__(self)
        self.plotter = seq_plotter
        
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.vsplitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.hsplitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.vsplitter)
        self.vsplitter.addWidget(self.hsplitter)
        
        self.analysis_plot = EvalPlotter()
        self.vsplitter.addWidget(self.analysis_plot)
        
        self.ptree = pt.ParameterTree()
        self.hsplitter.addWidget(self.ptree)
        
        self.table = pg.TableWidget()
        self.hsplitter.addWidget(self.table)
        
        self.hsplitter.setSizes([200, 300, 300])
        
        self.clear()
        
        self.params = TraceAnalyzerGroup(name='analyzers')
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
            for anal in self.params.children():
                data[anal.name()][i] = anal.process(t, d)
        self.table.setData(data)
        self.analysis_plot.update_data(data)
        

class TraceAnalyzerGroup(pt.parameterTypes.GroupParameter):
    need_update = QtCore.Signal()

    def __init__(self, **kwds):
        analyses = ['min', 'max', 'mean', 'exp tau', 'spike count', 'spike latency']
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
        kwds.update({'removable': True, 'renamable': True})
        childs = [
            dict(name='Input', type='list', values=kwds.pop('inputs')),
            dict(name='Type', type='list', value=kwds.pop('analysis_type'), values=['mean', 'min', 'max']),
            dict(name='Start', type='float', value=0, suffix='s', siPrefix=True, step=5e-3),
            dict(name='End', type='float', value=10e-3, suffix='s', siPrefix=True, step=5e-3),
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
            else:
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
        typ = self['Type']
        if typ == 'mean':
            return data.mean()
        elif typ == 'min':
            return data.min()
        elif typ == 'max':
            return data.max()
        

class EvalPlotter(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.xCode = QtGui.QLineEdit('cmd')
        self.yCode = QtGui.QLineEdit()
        self.plot = pg.PlotWidget()
        self.layout.addWidget(self.xCode, 0, 0)
        self.layout.addWidget(self.yCode, 0, 1)
        self.layout.addWidget(self.plot, 1, 0, 1, 2)
        self.xCode.textChanged.connect(self.replot)
        self.yCode.textChanged.connect(self.replot)
        
    def update_data(self, data):
        self.data = data
        self.replot()
        
    def replot(self):
        data = self.data
        ns = {}
        for k in data.dtype.names:
            ns[k] = data[k]
        try:
            xcode = str(self.xCode.text())
            ycode = str(self.yCode.text())
            if xcode == '' or ycode == '':
                return
            x = eval(xcode, ns)
            y = eval(ycode, ns)
        except:
            pg.debug.printExc('Error evaluating plot fields:')
        self.plot.plot(x, y, symbol='o', clear=True)
        
    
