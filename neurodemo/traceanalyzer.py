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
        self.splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.splitter)
        
        self.ptree = pt.ParameterTree()
        self.splitter.addWidget(self.ptree)
        
        self.table = pg.TableWidget()
        self.splitter.addWidget(self.table)
        self.splitter.setSizes([200, 600])
        
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
            self.plotter.plots[anal['Input']].addItem(anal.rgn)
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
            
        

class TraceAnalyzerGroup(pt.parameterTypes.GroupParameter):
    #inputs_changed = QtCore.Signal()
    #analyzers_changed = QtCore.Signal()
    need_update = QtCore.Signal()

    def __init__(self, **kwds):
        analyses = ['min', 'max', 'mean', 'exp tau', 'spike count', 'spike latency']
        self.inputs = []
        pt.parameterTypes.GroupParameter.__init__(self, addText='Add analysis..', addList=analyses, **kwds)

    def addNew(self, typ):
        param = TraceAnalyzerParameter(name=typ, analysis_type=typ, inputs=self.inputs, autoIncrementName=True)
        self.addChild(param)
        param.need_update.connect(self.need_update)
        #param.input_changed.connect(self.inputs_changed)
        #self.analyzers_changed.emit()
        self.need_update.emit()
        
    def set_inputs(self, inputs):
        self.inputs = list(inputs)
        self.inputs.remove('t')
        for ch in self.children():
            ch.set_input_list(inputs)


class TraceAnalyzerParameter(pt.parameterTypes.GroupParameter):
    #input_changed = QtCore.Signal(object, object)  # self, input
    need_update = QtCore.Signal(object)  # self

    def __init__(self, **kwds):
        kwds.update({'removable': True, 'renamable': True})
        childs = [
            dict(name='Input', type='list', values=kwds.pop('inputs')),
            dict(name='Type', type='list', value=kwds.pop('analysis_type'), values=['mean', 'min', 'max']),
            dict(name='Start', type='float', value=0, suffix='s', siPrefix=True),
            dict(name='End', type='float', value=10e-3, suffix='s', siPrefix=True),
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
        


