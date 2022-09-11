# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore

from .analysisplot import AnalysisPlot   # simpler code-based analyzer
from .traceanalyzer import TraceAnalyzer  # user friendly analyzer


class SequencePlotWindow(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.mode = 'ic'
        self.layout = QtGui.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        self.hold_check = QtGui.QCheckBox("Hold data")
        self.hold_check.setChecked(True)
        self.layout.addWidget(self.hold_check, 0, 0)
        self.clear_btn = QtGui.QPushButton("Clear data")
        self.layout.addWidget(self.clear_btn, 0, 1)
        self.clear_btn.clicked.connect(self.clear_data)
        self.plotsigns = {}
        self.splitter = QtGui.QSplitter(QtCore.Qt.Orientation.Vertical)
        self.layout.addWidget(self.splitter, 1, 0, 1, 2)
        
        self.plot_layout = pg.GraphicsLayoutWidget()
        self.splitter.addWidget(self.plot_layout)
        self.plots = {}
        
        self.analyzer = TraceAnalyzer(self)
        self.splitter.addWidget(self.analyzer)
        
        #self.analyzer = AnalysisPlot()
        #self.splitter.addWidget(self.analyzer)

    def add_plot(self, key, label):
        plot = self.plot_layout.addPlot(labels={'left': label, 'bottom': ('Time', 's')})
        if 'soma.V' in self.plots:
            plot.setXLink(self.plots['soma.V'])
            
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
        # check to see if mode has changed, and if so, clear the plot
        if self.mode != info['mode']:
            self.mode = info['mode']
            self.clear_data()
        
        if info['seq_len'] == 0:
            pen = 'w'
        else:
            pen = (info['seq_ind'], info['seq_len'] * 4./3.)
        
        for k, plt in self.plots.items():
            sign = 1.0
            if k in ["soma.IK.I", "soma.IKf.I", "soma.IKs.I", "soma.INa.I",
                "soma.IH.I", "soma.INa1.I"]:
                sign = -1.0   # flip sign of cation currents for display
            plt.plot(t, sign*data[k], pen=pen)
        
        try:
            self.analyzer.add_data(t, data, info)
        except:
            pg.debug.printExc('Error analyzing data:')
        self.show()
        
    def clear_data(self):
        for plt in self.plots.values():
            plt.clear()
        self.analyzer.clear()
    
    def plot_triggers(self, t, d):
        for k, plt in self.plots.items():
            plt.plot([t,t], [-d, d], pen=pg.mkPen(color="w", width=1.5))
