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
    def __init__(self):
        QtGui.QWidget.__init__(self)
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

    def clear(self):
        self.data = np.empty(0, dtype=[(str('command'), float)]) # note: no unicode allowed in py2/numpy dtypes
        self.table.clear()
        
    def add_data(self, t, v,i, info):
        mode, amp, cmd, seq_ind, seq_len = info
        rec = np.array([(amp,)], dtype=self.data.dtype)
        self.data = np.append(self.data, rec)
        self.table.setData(self.data)


#class TraceAnalyzerGroup(pt.GroupParameter):
    #def __init__(self):
        #analyses = ['min', 'max', 'mean', 'exp tau', 'spike count', 'spike latency']
        #pt.GroupParameter.__init__(self, addText='Add analysis..', addList=analyses)

    #def addNew(self, typ):
        #pass
    

#class TraceAnalyzerParameter(pt.GroupParameter):
    #def __init__(self):
        #pass
