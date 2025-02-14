import typing
import numpy as np
from . import qt
import pyqtgraph.parametertree as pt
import neurodemo.units as NU
from .ions import IonClass


class ChannelParameter(pt.parameterTypes.SimpleParameter):
    
    # emitted when a plot should be shown or hidden
    plots_changed = qt.Signal(object, object, object, object)  # self, channel, name, on/off
    
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
                self.channel.Erev = val
            elif param.name().startswith('Plot'):
                self.plots_changed.emit(self, self.channel, param.name()[5:], val)


class IonConcentrations(pt.parameterTypes.GroupParameter):
    """Group of parameters describing the concentrations (in and out)
    and reversal potential of an ion.
    """
    def __init__(self, ion: IonClass):
        self.ion = ion
        name = ion.name
        self.ion_params = [
            dict(name='[C]out', type='float', value=ion.Cout, suffix='mM', siPrefix=True, step=1),
            dict(name='[C]in', type='float', value=ion.Cin, suffix='mM', siPrefix=True, step=1),
        ]
        self.temperatureC = 37.
        self.ion.Erev = self.Nernst(ion)
        self.ion_params.append(dict(name="Erev", type='float', value=self.ion.Erev, suffix='V', siPrefix=True, readonly=True))
        pt.parameterTypes.GroupParameter.__init__(self, name=name, type='group',
                                                  children=self.ion_params,
                                                  expanded=False)
        self.sigTreeStateChanged.connect(self.treeChange)

    def treeChange(self, root, changes: list[tuple[pt.Parameter, str, typing.Any]]):
        for param, change, val in changes:
            if change != 'value':
                continue
            elif param is self.child('[C]out'):
                self.ion.Cout = val
                self['Erev'] = self.Nernst(self.ion)
            elif param is self.child('[C]in'):
                self.ion.Cin = val
                self['Erev'] = self.Nernst(self.ion)
    
    def updateTemperature(self, temp:float=37.):
        self.temperatureC = temp
        self['Erev']=self.Nernst(self.ion)

    def Nernst(self, ion):
        R = 8.135 # J/K/M
        F = 96840 # C/M
        RTF = R*(273.16+self.temperatureC)/F
        Er = (2.303*RTF/ion.valence)*np.log10(ion.Cout/ion.Cin)
        return Er
