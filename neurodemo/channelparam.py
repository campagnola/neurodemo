# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""
import numpy as np
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph.parametertree as pt
import neurodemo.units as NU

class ChannelParameter(pt.parameterTypes.SimpleParameter):
    
    # emitted when a plot should be shown or hidden
    plots_changed = QtCore.Signal(object, object, object, object)  # self, channel, name, on/off
    
    def __init__(self, channel):
        self.channel = channel
        name = channel.name
        ch_params = [
            dict(name='Gmax', type='float', value=channel.gmax, suffix='S', siPrefix=True, step=0.1, dec=True),
            dict(name='Erev', type='float', value=channel.erev, suffix='V', siPrefix=True, step=5*NU.mV),
            dict(name='Plot I', type='bool', value=False),
            dict(name='Plot G', type='bool', value=False),
            dict(name='Plot OP', type='bool', value=False),
        ]
        for sv in channel.difeq_state():
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
            elif param is self.child('Gmax'):
                self.channel.gmax = val
            elif param is self.child('Erev'):
                self.channel.erev = val
            elif param.name().startswith('Plot'):
                self.plots_changed.emit(self, self.channel, param.name()[5:], val)

class IonConcentrations(pt.parameterTypes.SimpleParameter):
    def __init__(self, ion):
        self.ion = ion
        name = ion.name
        print(ion)
        ion_params = [
            dict(name='Use C', type='bool', value=False),
            dict(name='Cout', type='float', value=ion.Cout, suffix='mM', siPrefix=True, step=1),
            dict(name='Cin', type='float', value=ion.Cin, suffix='mM', siPrefix=True, step=1),
        ]

        ion_params.append(dict(name="Erev", type='float', value=self.Nernst(ion), readonly=True))
        pt.parameterTypes.SimpleParameter.__init__(self, name=name, type='bool', 
                                                   value=ion.enabled, children=ion_params,
                                                   expanded=False)
        self.sigTreeStateChanged.connect(self.treeChange)

    def treeChange(self, root, changes):
        for param, change, val in changes:
            if change != 'value':
                continue
            print('param: ', param, change, val)
            if param is self:
                self.ion.enabled = val
            elif param is self.child('Cout'):
                self.ion.Cout = val
                self.ion.Erev=self.Nernst(self.ion)
            elif param is self.child('Cin'):
                self.ion.Cin = val
                self.ion.Erev=self.Nernst(self.ion)

    def Nernst(self, ion):
        Er = (58/ion.valence)*np.log10(ion.Cout/ion.Cin)
        return Er
