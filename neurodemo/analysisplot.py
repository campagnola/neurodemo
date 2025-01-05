# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""
from __future__ import division, unicode_literals
import sys
import numpy as np
import pyqtgraph as pg
from . import qt
import pyqtgraph.console
from pyqtgraph.metaarray import MetaArray
from .editor import Editor


analysis_code = '''
base = data['Time': 0.0:0.01]
pulse = data['Time': 0.04:0.06]
vbase = base['V'].mean(axis='Time')
vpulse = pulse['V'].mean(axis='Time')
ipulse = pulse['Ipip'].mean(axis='Time')
x = ipulse
y = vpulse - vbase
'''

class AnalysisPlot(qt.QSplitter):
    def __init__(self):
        qt.QSplitter.__init__(self, qt.Qt.Vertical)
        self.ns = {}
        
        self.plot = pg.PlotWidget()
        self.addWidget(self.plot)
        self.editor = Editor()
        self.addWidget(self.editor)
        self.editor.setText(analysis_code)
        
        self.output = qt.QTextEdit()
        self.output.setVisible(False)
        self.addWidget(self.output)
        
        self.btn_widget = qt.QWidget()
        self.addWidget(self.btn_widget)
        self.btn_layout = qt.QGridLayout()
        self.btn_widget.setLayout(self.btn_layout)
        
        self.replot_btn = qt.QPushButton('Replot')
        self.btn_layout.addWidget(self.replot_btn, 0, 0)
        self.replot_btn.clicked.connect(self.update_plot)
        
        self.console = pg.console.ConsoleWidget(namespace=self.ns)
        self.console_btn = qt.QPushButton('Console')
        self.console_btn.setCheckable(True)
        self.console_btn.toggled.connect(self.console.setVisible)
        self.btn_layout.addWidget(self.console_btn, 0, 1)
        
        self.clear()
        
    def add_data(self, t, data, info):
        v = data['soma.V']
        i = data['soma.PatchClamp.I']
        data = np.vstack([v, i])[:, np.newaxis, :]
        if self.ns['data'] is None:
            cols = [
                {'name': 'Vm', 'units': 'V'},
                {'name': 'Ipip', 'units': 'A'},
            ]
            self.ns['data'] = MetaArray(data, info=[
                {'name': 'Signal', 'cols': cols}, 
                {'name': 'Trial', 'values': [info['amp']], 'units': 'A'},
                {'name': 'Time', 'units': 's', 'values': t},
                {}
            ])
        else:
            data = np.concatenate([self.ns['data'].asarray(), data], axis=1)
            minfo = self.ns['data'].infoCopy()
            minfo[1]['values'] = np.append(minfo[1]['values'], info['amp'])
            self.ns['data'] = MetaArray(data, info=minfo)
        self.update_plot()
        
    def clear(self):
        self.ns.clear()
        self.ns['data'] = None
        self.plot.clear()
        
    def update_plot(self):
        code = self.editor.text()
        try:
            exec(code, self.ns)
            self.output.setVisible(False)
        except:
            pg.debug.printExc("Error in analysis code:")
            #sys.excepthook(*sys.exc_info())
            return
        if 'x' in self.ns or 'y' in self.ns:
            x = self.ns.get('x', None)
            y = self.ns.get('y', None)
            symbol = 'o' if len(y) < 100 else None
            self.plot.plot(x=x, y=y, symbol=symbol, clear=True)
            
            if isinstance(x, MetaArray):
                self.plot.setLabel('left', x._info[-1]['name'], x._info[-1]['units'])
            if isinstance(y, MetaArray):
                self.plot.setLabel('bottom', y._info[-1]['name'], y._info[-1]['units'])
        

