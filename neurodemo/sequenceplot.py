# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""
from __future__ import division, unicode_literals
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
#from .analysisplot import AnalysisPlot

# Alternate approach to analysis: a much more user-friendly interface. This is
# incomplete because it was too much effort compared to AnalysisPlot, and also
# robs students of the opportunity to get a little code exposure. 
from .traceanalyzer import TraceAnalyzer


class SequencePlotWindow(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.mode = 'ic'
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.hold_check = QtGui.QCheckBox("Hold data")
        self.hold_check.setChecked(True)
        self.layout.addWidget(self.hold_check, 0, 0)
        self.clear_btn = QtGui.QPushButton("Clear data")
        self.layout.addWidget(self.clear_btn, 0, 1)
        self.clear_btn.clicked.connect(self.clear_data)
        
        self.splitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.layout.addWidget(self.splitter, 1, 0, 1, 2)
        
        self.plot_layout = pg.GraphicsLayoutWidget()
        self.splitter.addWidget(self.plot_layout)
        #self.vplot = self.plot_layout.addPlot(0, 0, labels={'left': ('Membrane Voltage', 'V'), 'bottom': ('Time', 's')})
        #self.iplot = self.plot_layout.addPlot(1, 0, labels={'left': ('Pipette Current', 'A'), 'bottom': ('Time', 's')})
        #self.iplot.setXLink(self.vplot)
        self.plots = {}
        
        self.analyzer = TraceAnalyzer()
        self.splitter.addWidget(self.analyzer)
        
        #self.analyzer = AnalysisPlot()
        #self.splitter.addWidget(self.analyzer)

    def add_plot(self, key, label):
        plot = self.plot_layout.addPlot(labels={'left': label, 'bottom': ('Time', 's')})
        self.plot_layout.nextRow()
        self.plots[key] = plot

    def remove_plot(self, key):
        plot = self.plots.pop(key)
        self.plot_layout.removeItem(plot)
        plot.hide()
        plot.setParentItem(None)
        
    def plot(self, t, data, info):
        if not self.hold_check.isChecked():
            self.clear_data()
        if self.mode != info['mode']:
            self.mode = info['mode']
            self.clear_data()
        
        if info['seq_len'] == 0:
            pen = 'w'
        else:
            pen = (info['seq_ind'], info['seq_len'] * 4./3.)
        
        for k, plt in self.plots.items():
            plt.plot(t, data[k], pen=pen)
        
        try:
            self.analyzer.add_data(t, data, info)
        except:
            pg.debug.printExc('Error analyzing data:')
        self.show()
        
    def clear_data(self):
        for plt in self.plots.values():
            plt.clear()
        self.analyzer.clear()
