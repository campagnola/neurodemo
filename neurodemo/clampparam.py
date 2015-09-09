# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""
from __future__ import division, unicode_literals
import numpy as np
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph.parametertree as pt
from .sequenceplot import SequencePlotWindow
from .units import *

class ClampParameter(pt.parameterTypes.SimpleParameter):
    # emitted when a plot should be shown or hidden
    plots_changed = QtCore.Signal(object, object, object, object)  # self, channel, name, on/off
    
    def __init__(self, clamp, dt):
        self.clamp = clamp
        self.dt = dt
        self.plot_win = SequencePlotWindow()
        self.triggers = []  # items are (trigger_time, pointer, trigger_buffer, (mode, amp, cmd, seq_ind, seq_len))
        pt.parameterTypes.SimpleParameter.__init__(self, name='Patch Clamp', type='bool', value=True, children=[
            dict(name='Mode', type='list', values={'Current Clamp': 'ic', 'Voltage Clamp': 'vc'}, value='ic'),
            dict(name='Holding', type='float', value=0, suffix='A', siPrefix=True, step=10*pA),
            dict(name='Pipette Capacitance', type='float', value=clamp.cpip, limits=[0.01*pF, None], suffix='F', siPrefix=True, dec=True, step=0.5),
            dict(name='Access Resistance', type='float', value=clamp.ra, limits=[10*kOhm, None], suffix='Ω', siPrefix=True, step=0.5, dec=True),
            dict(name='Plot Current', type='bool', value=False),
            dict(name='Pulse', type='group', children=[
                dict(name='Capture Results', type='bool', value=False),
                dict(name='Pulse Once', type='action'),
                dict(name='Amplitude', type='float', value=50*pA, suffix='A', siPrefix=True),
                dict(name='Pre-delay', type='float', value=20*ms, suffix='s', siPrefix=True, limits=[0, None]),
                dict(name='Duration', type='float', value=50*ms, suffix='s', siPrefix=True, limits=[0, None]),
                dict(name='Post-delay', type='float', value=50*ms, suffix='s', siPrefix=True, limits=[0, None]),
                dict(name='Pulse Sequence', type='action'),
                dict(name='Start Amplitude', type='float', value=-50*pA, suffix='A', siPrefix=True),
                dict(name='Stop Amplitude', type='float', value=50*pA, suffix='A', siPrefix=True),
                dict(name='Pulse Number', type='int', value=11, limits=[2,None]),
            ]),
        ])
        self.sigTreeStateChanged.connect(self.treeChange)
        self.child('Pulse', 'Pulse Once').sigActivated.connect(self.pulse_once)
        self.child('Pulse', 'Pulse Sequence').sigActivated.connect(self.pulse_sequence)

    def treeChange(self, root, changes):
        for param, change, val in changes:
            if change != 'value':
                continue
            if param is self:
                self.clamp.enabled = val
            elif param is self.child('Mode'):
                self.set_mode(val)
            elif param is self.child('Holding'):
                self.clamp.set_holding(self.mode(), val)
            elif param is self.child('Pipette Capacitance'):
                self.clamp.cpip = val
            elif param is self.child('Access Resistance'):
                self.clamp.ra = val
            elif param.name().startswith('Plot'):
                self.plots_changed.emit(self, self.clamp, 'I', val)

    def mode(self):
        return self['Mode']

    def set_mode(self, mode):
        self.clamp.set_mode(mode)
        suff = {'ic': 'A', 'vc': 'V'}[mode]
        amp, start, stop, step = {'ic': (-10*pA, -100*pA, 100*pA, 10*pA), 
                                  'vc': (-10*mV, -80*mV, 50*mV, 5*mV)}[mode]
        self.sigTreeStateChanged.disconnect(self.treeChange)
        try:
            self.child('Holding').setOpts(suffix=suff, value=self.clamp.holding[mode], step=step)
            self.child('Pulse', 'Amplitude').setOpts(suffix=suff, value=amp, step=step)
            self.child('Pulse', 'Start Amplitude').setOpts(suffix=suff, value=start, step=step)
            self.child('Pulse', 'Stop Amplitude').setOpts(suffix=suff, value=stop, step=step)
        finally:
            self.sigTreeStateChanged.connect(self.treeChange)
            
    def pulse_template(self):
        d1 = self['Pulse', 'Pre-delay']
        d2 = self['Pulse', 'Duration']
        d3 = self['Pulse', 'Post-delay']
        dur = d1 + d2 + d3
        npts = dur / self.dt
        cmd = np.empty(npts)
        i1 = d1 / self.dt
        i2 = i1 + d2 / self.dt
        cmd[:] = self['Holding']
        return cmd, i1, i2
        
    def pulse_once(self):
        cmd, i1, i2 = self.pulse_template()
        amp = self['Pulse', 'Amplitude']
        cmd[i1:i2] += amp
        t = self.clamp.queue_command(cmd, self.dt)
        if self['Pulse', 'Capture Results']:
            self.triggers.append([t, 0, np.empty((len(cmd), 2)), (self.mode(), amp, cmd, 0, 0)]) 
    
    def pulse_sequence(self):
        cmd, i1, i2 = self.pulse_template()
        cmds = []
        amps = np.linspace(self['Pulse', 'Start Amplitude'],
                           self['Pulse', 'Stop Amplitude'],
                           self['Pulse', 'Pulse Number'])
        for amp in amps:
            cmd2 = cmd.copy()
            cmd2[i1:i2] += amp
            cmds.append(cmd2)
        
        times = self.clamp.queue_commands(cmds, self.dt)
        for i, t in enumerate(times):
            self.triggers.append([t, 0, np.empty((len(cmd), 2)), (self.mode(), amps[i], cmds[i], i, len(amps))])
        
    def new_result(self, result):
        if len(self.triggers) == 0:
            return
        try:
            vm = result['soma.Vm'][:-1]
            ip = result['soma.PatchClamp.I'][:-1]
            t = result['t'][:-1]
        except KeyError:
            return
        self.plot_triggered(vm, ip, t)
        
    def plot_triggered(self, vm, ip, t):
        if len(self.triggers) == 0:
            return
        tt, ptr, data, info = self.triggers[0]
        if tt > t[-1]:
            # no triggers ready
            return
        
        # Copy data from result to trigger buffer
        i = max(0, np.round((tt - t[0]) / self.dt)) # index of trigger within new data
        npts = min(len(data)-ptr, len(vm)-i) # number of samples to copy from new data
        data[ptr:ptr+npts, 0] = vm[i:i+npts] 
        data[ptr:ptr+npts, 1] = ip[i:i+npts]
            
        ptr += npts
        if ptr >= data.shape[0]:
            # If the trigger buffer is full, plot and remove
            self.plot_win.plot(np.arange(data.shape[0])*self.dt, data[:,0], data[:,1], info)
            self.triggers.pop(0)
            if len(vm) > npts:
                # If there is data left over, try feeding it to the next trigger
                self.plot_triggered(vm[i+npts:], ip[i+npts:], t[i+npts:])
        else:
            # otherwise, update the pointer and wait for the next result
            self.triggers[0][1] = ptr
