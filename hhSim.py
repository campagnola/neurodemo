# -*- coding: utf-8 -*-
"""
Simple neuron simulator for Python.
Also simulates voltage clamp and current clamp with access resistance.

Luke Campagnola 2015
"""

from collections import OrderedDict
import numpy as np
import scipy.integrate
from PyQt4 import QtGui, QtCore


# Define a set of scaled unit symbols to make the code more clear
for unit in 'msVAΩFS':
    for pfx, val in [('p', -12), ('n', -9), ('u', -6), ('m', -3), ('c', -2), ('k', 3), ('M', 6), ('G', 9)]:
        locals()[pfx+unit] = 10**val


#Ik2 = 24.3

# alpha synapse
#Alpha_t0 = 500.  # msec
#Alpha_tau = 2.0
#gAlpha = 1e-3 * Area/cm**2
#EAlpha = -7e-3  # V



def IAlpha(Vm, t):
    if t < Alpha_t0:
        return 0.
    else:
        # g = gmax * (t - onset)/tau * exp(-(t - onset - tau)/tau)
        tn = t - Alpha_t0
        if tn > 10.0 * Alpha_tau:
            return 0.
        else:
            return gAlpha * (Vm - EAlpha)*(tn/Alpha_tau) * np.exp(-(tn-Alpha_tau)/Alpha_tau)


class Sim(object):
    """Simulator for a collection of objects that derive from SimObject
    """
    def __init__(self, objects=None):
        if objects is None:
            objects = []
        self._objects = objects
        self._all_objs = None
        self._time = 0.0
        self.temp = 6.3

    def add(self, obj):
        assert obj._sim is None
        obj._sim = self
        self._objects.append(obj)

    def all_objects(self):
        if self._all_objs is None:
            objs = []
            for o in self._objects:
                objs.extend(o.all_objects())
            self._all_objs = objs
        return self._all_objs
    
    def state(self):
        s = OrderedDict()
        for o in self.all_objects():
            for k,v in o.state().items():
                s[k] = v
        return s
    
    @property
    def state_vars(self):
        return list(self.state().keys())

    def run(self, dt=0.1*ms, dur=100*ms, **args):
        self._all_objs = None
        svars = []
        for o in self.all_objects():
            for k in o.state_vars:
                svars.append((o, k))
        self._simstate = SimState(svars)
        npts = int(dur/dt)
        t = np.linspace(0, dur, npts)

        # Run the simulation
        init_state = list(self.state().values())
        global result
        result, info = scipy.integrate.odeint(self.derivatives, init_state, t, (dt,),
                                              rtol=1e-6, atol=1e-6, hmax=5e-2, full_output=1, **args)
        p = 0
        for o in self.all_objects():
            nvar = len(o.state_vars)
            o.update_result(result[:, p:p+nvar])
            p += nvar
            
        self._last_run_time = t

    def derivatives(self, state, t, dt):
        objs = self.all_objects()
        #for o in objs:
            #nvars = len(o.state_vars)
            #substate = state[:nvars]
            #state = state[nvars:]
            #o.set_current_state(substate)
        self._simstate.state = state
        d = []
        for o in objs:
            d.extend(o.derivatives(self._simstate, t))
            
        return d
    

class SimState(object):
    def __init__(self, keys):
        self.keys = keys
        self.indexes = dict([(k,i) for i,k in enumerate(keys)])
        self.state = None
        
    def __getitem__(self, key):
        return self.state[self.indexes[key]]


class SimObject(object):
    """
    Base class for objects that participate in integration by providing a set
    of state variables and their derivatives.
    """
    def __init__(self, init_state):
        self._sim = None
        self._init_state = init_state
        self._state_vars = tuple(init_state.keys())
        self._current_state = init_state.copy()
        self._sub_objs = []
        self.records = []
        self._rec_dtype = [(sv, float) for sv in init_state.keys()]
    
    def state(self):
        return self._current_state

    def all_objects(self):
        objs = [self]
        for o in self._sub_objs:
            objs.extend(o.all_objects())
        return objs
        
    @property
    def state_vars(self):
        return self._state_vars

    def update_result(self, result):
        self._last_result = result
        self.set_current_state(result[-1])

    def set_current_state(self, state):
        for i,k in enumerate(self.state_vars):
            self._current_state[k] = state[i]

    def derivatives(self, state, t):
        """Return derivatives of all state variables.
        
        Must be reimplemented in subclasses.
        """
        raise NotImplementedError()
    
    @property
    def sim(self):
        return self._sim


class Mechanism(SimObject):
    """Base class for simulation objects that interact with a section's
    membrane--channels, electrodes, etc.
    """
    def __init__(self, init_state, section=None):
        SimObject.__init__(self, init_state)
        self._section = section
        
    def current(self):
        """Return the membrane current being passed by this mechanism.
        
        Must be implemented in subclasses.
        """
        raise NotImplementedError()

    @property
    def section(self):
        return self._section
    
    @property
    def sim(self):
        return self.section.sim


class Channel(Mechanism):
    """Base class for simple ion channels.
    """
    # precomputed rate constant tables
    rates = None
    
    @classmethod
    def compute_rates(cls):
        return
        
    def __init__(self, gbar, init_state):
        Mechanism.__init__(self, init_state)
        self.gbar = gbar
        self._g = None
        if self.rates is None:
            type(self).compute_rates()
              
    @property
    def g(self):
        if self._g is None:
            self._g = self.gbar * self.section.area
        return self._g

    def current(self, state, t):
        vm = state[self.section, 'Vm']
        op = self.open_probability(state)
        return -self.g * op * (vm - self.erev)

    @staticmethod
    def interpolate_rates(rates, val, minval, step):
        """Helper function for interpolating kinetic rates from precomputed
        tables.
        """
        i = (val - minval) / step
        i1 = int(i)
        i2 = i1 + 1
        s = i2 - i
        if i1 < 0:
            return rates[0]
        elif i2 >= len(rates):
            return rates[-1]
        else:
            return rates[i1] * s + rates[i2] * (1-s)


class Section(SimObject):
    def __init__(self, radius=10*um):
        self.area = 4 * 3.1415926 * radius**2
        self.cm = self.area * (1 * uF / cm**2)
        self.ek = -77*mV
        self.ena = 50*mV
        self.ecl = -70*mV
        init_state = OrderedDict([('Vm', -65*mV)])
        SimObject.__init__(self, init_state)
        self.mechanisms = []

    def add(self, mech):
        assert mech._section is None
        mech._section = self
        self.mechanisms.append(mech)
        self._sub_objs.append(mech)

    def derivatives(self, state, t):
        Im = 0
        for mech in self.mechanisms:
            Im += mech.current(state, t)
            
        dv = Im / self.cm
        return [dv]
        

class Leak(Mechanism):
    def __init__(self, g, erev=-55*mV):
        self.g = g
        self.erev = erev
        Mechanism.__init__(self, {})
        
    def current(self, state, t):
        vm = state[self.section, 'Vm']
        return self.g * (self.erev - vm)

    def derivatives(self, state, t):
        return []


class MultiClamp(Mechanism):
    def __init__(self, mode='ic', cmd=None, dt=None, ra=5*MΩ, cpip=3*pF):
        self.ra = ra
        self.cpip = cpip
        self.mode = mode
        self.cmd = cmd
        self.dt = dt
        self.gain = 50e-6  # arbitrary VC gain
        init_state = OrderedDict([('Ve', -65*mV)])
        Mechanism.__init__(self, init_state)

    def set_command(self, cmd, dt):
        self.cmd = cmd
        self.dt = dt

    def pipette_current(self):
        """Compute current through pipette from most recent simulation run.
        """
        ve = self._last_result[:,0]
        vm = self.section._last_result[:,0]
        return (ve-vm) / self.ra

    def current(self, state, t):
        # Compute current through tip of pipette at this timestep
        vm = state[self.section, 'Vm']
        ve = state[self, 'Ve']
        return (ve - vm) / self.ra
    
    def derivatives(self, state, t):
        ## Select between VC and CC
        cmd = self.cmd
        if cmd is None:
            cmd = 0
            mode = 'ic'
        else:
            # interpolate command -- sharp steps confuse the integrator.
            fInd = t / self.dt
            ind = min(len(cmd)-1, np.floor(fInd))
            ind2 = min(len(cmd)-1, ind+1)
            s = fInd - ind
            cmd = cmd[ind] * (1-s) + cmd[ind2] * s
            
        # determine current generated by voltage clamp 
        if self.mode == 'vc':
            ve = state[self, 'Ve']
            cmd = (cmd-ve) * self.gain
        
        # Compute change in electrode potential
        dve = (cmd - self.current(state, t)) / self.cpip
        return [dve]


class HHK(Channel):
    """Hodgkin-Huxley K channel.
    """
    @classmethod
    def compute_rates(cls):
        cls.rates_vmin = -100
        cls.rates_vstep = 0.1
        vm = np.arange(cls.rates_vmin, cls.rates_vmin+400, cls.rates_vstep)
        cls.rates = np.empty((len(vm), 2))
        cls.rates[:,0] = (0.1 - 0.01*vm) / (np.exp(1.0 - 0.1*vm) - 1.0)
        cls.rates[:,1] = 0.125 * np.exp(-vm / 80.)
        
    def __init__(self, gbar=12*mS/cm**2):
        init_state = OrderedDict([('n', 0.3)]) 
        Channel.__init__(self, gbar, init_state)
        self.shift = 0
        
    @property
    def erev(self):
        return self.section.ek
        
    def open_probability(self, state):
        return state[self, 'n']**4

    def derivatives(self, state, t):
        # temperature dependence of rate constants
        q10 = 3 ** ((self.sim.temp-6.3) / 10.)
        vm = state[self.section, 'Vm'] - self.shift
        
        vm = vm + 65e-3   ## gating parameter eqns assume resting is 0mV
        vm *= 1000.   ##  ..and that Vm is in mV
        
        n = state[self, 'n']
        
        # disabled for now -- does not seem to improve speed.
        #an, bn = self.interpolate_rates(self.rates, vm, self.rates_vmin, self.rates_vstep)
        
        an = (0.1 - 0.01*vm) / (np.exp(1.0 - 0.1*vm) - 1.0)
        bn = 0.125 * np.exp(-vm / 80.)
        dn = q10 * (an * (1.0 - n) - bn * n)
        return [dn*1e3]
                 

class HHNa(Channel):
    """Hodgkin-Huxley Na channel.
    """
    @classmethod
    def compute_rates(cls):
        cls.rates_vmin = -100
        cls.rates_vstep = 0.1
        vm = np.arange(cls.rates_vmin, cls.rates_vmin+400, cls.rates_vstep)
        cls.rates = np.empty((len(vm), 4))
        cls.rates[:,0] = (2.5-0.1*vm) / (np.exp(2.5-0.1*vm) - 1.0)
        cls.rates[:,1] = 4. * np.exp(-vm / 18.)
        cls.rates[:,2] = 0.07 * np.exp(-vm / 20.)
        cls.rates[:,3] = 1.0 / (np.exp(3.0 - 0.1 * vm) + 1.0)
        
    def __init__(self, gbar=40*mS/cm**2):
        init_state = OrderedDict([('m', 0.05), ('h', 0.6)]) 
        Channel.__init__(self, gbar, init_state)
        self.shift = 0
        
    @property
    def erev(self):
        return self.section.ena
        
    def open_probability(self, state):
        return state[self, 'm']**3 * state[self, 'h']

    def derivatives(self, state, t):
        # temperature dependence of rate constants
        q10 = 3 ** ((self.sim.temp-6.3) / 10.)
        vm = state[self.section, 'Vm'] - self.shift
        m = state[self, 'm']
        h = state[self, 'h']

        vm = vm + 65e-3   ## gating parameter eqns assume resting is 0mV
        vm *= 1000.   ##  ..and that Vm is in mV
        
        # disabled for now -- does not seem to improve speed.
        #am, bm, ah, bh = self.interpolate_rates(self.rates, vm, self.rates_vmin, self.rates_vstep)
        
        am = (2.5-0.1*vm) / (np.exp(2.5-0.1*vm) - 1.0)
        bm = 4. * np.exp(-vm / 18.)
        dm = q10 * (am * (1.0 - m) - bm * m)
        
        ah = 0.07 * np.exp(-vm / 20.)
        bh = 1.0 / (np.exp(3.0 - 0.1 * vm) + 1.0)
        dh = q10 * (ah * (1.0 - h) - bh * h)

        return [dm*1e3, dh*1e3]
                 

class IH(Mechanism):
    """Ih from Destexhe 1993
    """
    def __init__(self, gbar=30*mS/cm**2):
        init_state = OrderedDict([('f', 0), ('s', 0)]) 
        Mechanism.__init__(self, init_state)
        self.gbar = gbar
        self.erev = -43*mV
        self._g = None
        self.shift = 0
        
    @property
    def g(self):
        if self._g is None:
            self._g = self.gbar * self.section.area
        return self._g
    
    def current(self, state, t):
        vm = state[self.section, 'Vm']
        f = state[self, 'f']
        s = state[self, 's']
        return -self.g * f * s * (vm - self.erev)
    
    def derivatives(self, state, t):
        vm = state[self.section, 'Vm'] - self.shift
        f = state[self, 'f']
        s = state[self, 's']
        
        #vm = vm + 65e-3   ## gating parameter eqns assume resting is 0mV
        vm *= 1000.   ##  ..and that Vm is in mV
        Hinf = 1.0 / (1.0 + np.exp((vm + 68.9) / 6.5))
        tauF = np.exp((vm + 158.6)/11.2) / (1.0 + np.exp((vm + 75.)/5.5))
        tauS = np.exp((vm + 183.6) / 15.24)
        df = (Hinf - f) / tauF
        ds = (Hinf - s) / tauS
        return [df*1e3, ds*1e3]


class LGNa(Mechanism):
    """Cortical sodium channel (Lewis & Gerstner 2002, p.124)
    """
    def __init__(self, gbar=112.5*mS/cm**2):
        init_state = OrderedDict([('m', 0.05), ('h', 0.6)]) 
        Mechanism.__init__(self, init_state)
        self.gbar = gbar
        self.erev = 74*mV
        self._g = None
        
    @property
    def g(self):
        if self._g is None:
            self._g = self.gbar * self.section.area
        return self._g

    def current(self, state, t):
        vm = state[self.section, 'Vm']
        m = state[self, 'm']
        h = state[self, 'h']
        return -self.g * m**3 * h * (vm - self.section.ena)

    def derivatives(self, state, t):
        # temperature dependence of rate constants
        q10 = 3 ** ((self.sim.temp-6.3) / 10.)
        vm = state[self.section, 'Vm'] - self.shift
        m = state[self, 'm']
        h = state[self, 'h']

        vm = vm + 65e-3   ## gating parameter eqns assume resting is 0mV
        vm *= 1000.   ##  ..and that Vm is in mV
        
        # disabled for now -- does not seem to improve speed.
        #am, bm, ah, bh = self.interpolate_rates(self.rates, vm, self.rates_vmin, self.rates_vstep)
        
        am = (2.5-0.1*vm) / (np.exp(2.5-0.1*vm) - 1.0)
        bm = 4. * np.exp(-vm / 18.)
        dm = q10 * (am * (1.0 - m) - bm * m)
        
        ah = 0.07 * np.exp(-vm / 20.)
        bh = 1.0 / (np.exp(3.0 - 0.1 * vm) + 1.0)
        dh = q10 * (ah * (1.0 - h) - bh * h)

        return [dm*1e3, dh*1e3]




def run(sim, dt=1e-4, mode='ic', cmd=None, dur=None):
    """
    Return array of Vm or Im values.        
    """
    dt = dt * 1e3  ## convert s -> ms
    
    if dur is None and cmd is not None:
        dur = dt*len(cmd)
    
    result = neuron.run(cmd=cmd, mode=mode, dt=dt, dur=dur)
    
    if mode == 'ic':
        out = result#[:,2] + np.random.normal(size=len(data), scale=0.3e-3)
    elif mode == 'vc':
        out = result#[:,1] + np.random.normal(size=len(data), scale=3.e-12)
    
    return out


if __name__ == '__main__':
    import pyqtgraph as pg
    from pyqtgraph.Qt import QtGui
    pg.setConfigOption('antialias', True)
    app = QtGui.QApplication([])
    
    sim = Sim()
    neuron = Section()
    neuron.add(Leak(g=0.1e-3 * neuron.area/cm**2))
    neuron.add(HHK())
    neuron.add(HHNa())
    #neuron.add(IH())
    clamp = MultiClamp(mode='ic')
    neuron.add(clamp)
    sim.add(neuron)
    
    win = pg.GraphicsWindow()
    win.resize(1000, 600)
    win.setWindowTitle('Testing hhSim.py')
    p1 = win.addPlot(title='IC', labels={'left': ('Vm', 'V')})
    p2 = win.addPlot(labels={'left': ('Ipip', 'A')}, row=1, col=0)
    win.ci.layout.setRowFixedHeight(1, 150)
    
    dur = 100 * ms
    dt = 1e-5
    npts = int(dur / dt)
    x1 = int(20*ms / dt)
    x2 = int(80*ms / dt)
    x = np.linspace(-200, 200, 11) * pA
    #x = [-200*pA]
    cmd = np.zeros((len(x), npts)) #*-65e-3
    data = np.zeros((len(x), npts, 9))
    for i, v in enumerate(x):
        print('V: ', v)
        cmd[i, x1:x2] = v
        clamp.set_command(cmd[i], dt)
        #data[i] = run(neuron, mode='ic', dt=dt, cmd=cmd[i])
        sim.run(dt=dt, dur=dur)
        data = neuron._last_result[:,0]
        t = sim._last_run_time
        p1.plot(t, data, pen=(i, 15))
        p2.plot(t, cmd[i], pen=(i, 15))

    #win2 = pg.GraphicsWindow()
    #for i, k in enumerate(('t', 'Im', 'Ve', 'Vm', 'm', 'h', 'n', 'f', 's')):
        #p = win2.addPlot(labels={'left': k})
        #win2.nextRow()
        #p.plot(data[0,:,i])

    import sys
    #if sys.flags.interactive == 0:
        #QtGui.QApplication.instance().exec_()
