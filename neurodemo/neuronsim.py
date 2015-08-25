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
from .units import *

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
    def __init__(self, objects=None, temp=37.0, dt=10*us):
        if objects is None:
            objects = []
        self._objects = objects
        self._all_objs = None
        self._time = 0.0
        self.temp = temp
        self.dt = dt

    def add(self, obj):
        assert obj._sim is None
        obj._sim = self
        self._objects.append(obj)
        return obj

    def all_objects(self):
        if self._all_objs is None:
            objs = []
            print('---------------')
            for o in self._objects:
                if not o.enabled:
                    continue
                objs.extend(o.all_objects())
            for o in objs:
                print('   ' + str(o))
            self._all_objs = objs
        return self._all_objs
    
    def state(self):
        s = OrderedDict()
        for o in self.all_objects():
            for k,v in o.state().items():
                s[(o,k)] = v
        return s
    
    @property
    def state_vars(self):
        return list(self.state().keys())

    @property
    def time(self):
        return self._time

    def run(self, samples=1000, **kwds):
        """Run the simulation until a number of *samples* have been acquired.
        
        Extra keyword arguments are passed to `scipy.integrate.odeint()`.
        """
        opts = {'rtol': 1e-6, 'atol': 1e-6, 'hmax': self.dt*50, 'full_output': 1}
        opts.update(kwds)
        
        self._all_objs = None
        svars = []
        if len(self.all_objects()) == 0:
            raise RuntimeError("No objects added to simulation.")
            
        for o in self.all_objects():
            for k in o.state_vars:
                svars.append((o, k))
        self._simstate = SimState(svars)
        t = np.arange(1, samples+1) * self.dt + self._time

        # Run the simulation
        init_state = list(self.state().values())
        result, info = scipy.integrate.odeint(self.derivatives, init_state, t, **opts)
        p = 0
        for o in self.all_objects():
            nvar = len(o.state_vars)
            o.update_result(result[:, p:p+nvar])
            p += nvar
            
        self._time = t[-1]
        self._last_run_time = t
        return SimState(svars, result.T, t=t) 

    def derivatives(self, state, t):
        objs = self.all_objects()
        self._simstate.state = state
        self._simstate.extra['t'] = t
        d = []
        for o in objs:
            d.extend(o.derivatives(self._simstate, t))
            
        return d
    

class SimState(object):
    """Contains the state of all variables in the simulation.
    
    During simulation runs, this is used to carry information about all
    variables at the current timepoint. After the simulation finishes, this is
    used to carry all state variable data collected during the simulation.
    """
    def __init__(self, keys, state=None, **extra):
        self.keys = keys
        self.indexes = {}
        for i, k in enumerate(keys):
            self.indexes[k] = i
            self.indexes[(k[0].name, k[1])] = i
        self.state = state
        self.extra = extra
        
    def __getitem__(self, key):
        if key in self.indexes:
            return self.state[self.indexes[key]]
        else:
            return self.extra[key]


class SimObject(object):
    """
    Base class for objects that participate in integration by providing a set
    of state variables and their derivatives.
    """
    instance_count = 0
    
    def __init__(self, init_state, name=None):
        self._sim = None
        if name is None:
            i = self.instance_count
            type(self).instance_count = i + 1
            name = type(self).__name__ + '%d' % i
        self.name = name
        self.enabled = True
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
            if not o.enabled:
                continue
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
    def __init__(self, init_state, section=None, **kwds):
        SimObject.__init__(self, init_state, **kwds)
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
        
    def __init__(self, gbar, init_state, **kwds):
        Mechanism.__init__(self, init_state, **kwds)
        self.gbar = gbar
        self._g = None
        if self.rates is None:
            type(self).compute_rates()
              
    @property
    def g(self):
        if self._g is None:
            self._g = self.gbar * self.section.area
        return self._g

    def conductance(self, state):
        op = self.open_probability(state)
        return self.g * op

    def current(self, state, t=None):
        vm = state[self.section, 'Vm']
        g = self.conductance(state)
        return -g * (vm - self.erev)

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
    def __init__(self, radius=10*um, vm=-65*mV, **kwds):
        self.area = 4 * 3.1415926 * radius**2
        self.cm = self.area * (1 * uF / cm**2)
        self.ek = -77*mV
        self.ena = 50*mV
        self.ecl = -70*mV
        init_state = OrderedDict([('Vm', vm)])
        SimObject.__init__(self, init_state, **kwds)
        self.mechanisms = []

    def add(self, mech):
        assert mech._section is None
        mech._section = self
        self.mechanisms.append(mech)
        self._sub_objs.append(mech)
        return mech

    def derivatives(self, state, t):
        Im = 0
        for mech in self.mechanisms:
            Im += mech.current(state, t)
            
        dv = Im / self.cm
        return [dv]
        

class MultiClamp(Mechanism):
    def __init__(self, mode='ic', cmd=None, dt=None, ra=5*MΩ, cpip=3*pF, **kwds):
        self.ra = ra
        self.cpip = cpip
        self.mode = mode
        self.cmd = None
        if cmd is not None:
            self.set_command(cmd, dt)
        self.gain = 50e-6  # arbitrary VC gain
        init_state = OrderedDict([('Ve', -65*mV)])
        Mechanism.__init__(self, init_state, **kwds)

    def set_command(self, cmd, dt, start=0):
        self.cmd = cmd
        self.dt = dt
        self.cmd_ptr = 0
        self.start = start

    def pipette_current(self):
        """Compute current through pipette from most recent simulation run.
        """
        ve = self._last_result[:,0]
        vm = self.section._last_result[:,0]
        return (ve-vm) / self.ra

    def current(self, state, t=None):
        # Compute current through tip of pipette at this timestep
        vm = state[self.section, 'Vm']
        ve = state[self, 'Ve']
        return (ve - vm) / self.ra
    
    def derivatives(self, state, t):
        t = t - self.start
        
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


class Leak(Channel):
    def __init__(self, gbar=0.1*mS/cm**2, erev=-55*mV, **kwds):
        Channel.__init__(self, gbar, {}, **kwds)
        self.erev = erev

    def open_probability(self, state):
        return 1

    def derivatives(self, state, t):
        return []


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
        
    def __init__(self, gbar=12*mS/cm**2, **kwds):
        init_state = OrderedDict([('n', 0.3)]) 
        Channel.__init__(self, gbar, init_state, **kwds)
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
        
    def __init__(self, gbar=40*mS/cm**2, **kwds):
        init_state = OrderedDict([('m', 0.05), ('h', 0.6)]) 
        Channel.__init__(self, gbar, init_state, **kwds)
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
                 

class IH(Channel):
    """Ih from Destexhe 1993
    """
    def __init__(self, gbar=30*mS/cm**2, **kwds):
        init_state = OrderedDict([('f', 0), ('s', 0)]) 
        Channel.__init__(self, gbar, init_state, **kwds)
        self.erev = -43*mV
        self.shift = 0
        
    def open_probability(self, state):
        return state[self, 'f'] * [self, 's']
    
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


class LGNa(Channel):
    """Cortical sodium channel (Lewis & Gerstner 2002, p.124)
    """
    def __init__(self, gbar=112.5*mS/cm**2, **kwds):
        init_state = OrderedDict([('m', 0.019), ('h', 0.876)]) 
        Channel.__init__(self, gbar, init_state, **kwds)
        self.erev = 74*mV
        
    def open_probability(self, state):
        return state[self, 'm']**3 * state[self, 'h']

    def derivatives(self, state, t):
        # temperature dependence of rate constants
        # TODO: not sure about the base temp:
        q10 = 3 ** ((self.sim.temp - 37.) / 10.)
        
        vm = state[self.section, 'Vm']
        m = state[self, 'm']
        h = state[self, 'h']

        #vm = vm + 65e-3   ## gating parameter eqns assume resting is 0mV
        vm *= 1000.   ##  ..and that Vm is in mV
        
        am = (-3020 + 40 * vm)  / (1.0 - np.exp(-(vm - 75.5) / 13.5))
        bm = 1.2262 / np.exp(vm / 42.248)
        mtau = 1 / (am + bm)
        minf = am * mtau
        dm = q10 * (minf - m) / mtau
        
        ah = 0.0035 / np.exp(vm / 24.186)
        # note: bh as originally written causes integration failures; we use
        # an equivalent expression that behaves nicely under floating point stress.
        #bh = (0.8712 + 0.017 * vm) / (1.0 - np.exp(-(51.25 + vm) / 5.2))
        bh = 0.017 * (51.25 + vm) / (1.0 - np.exp(-(51.25 + vm) / 5.2))
        htau = 1 / (ah + bh)
        hinf = ah * htau
        dh = q10 * (hinf - h) / htau

        return [dm*1e3, dh*1e3]


class LGKfast(Channel):
    """Cortical fast potassium channel (Lewis & Gerstner 2002, p.124)
    """
    def __init__(self, gbar=225*mS/cm**2, **kwds):
        init_state = OrderedDict([('n', 0.00024)]) 
        Channel.__init__(self, gbar, init_state, **kwds)
        self.erev = -90*mV
        
    def open_probability(self, state):
        return state[self, 'n']**2

    def derivatives(self, state, t):
        # temperature dependence of rate constants
        # TODO: not sure about the base temp:
        q10 = 3 ** ((self.sim.temp - 37.) / 10.)
        
        vm = state[self.section, 'Vm']
        n = state[self, 'n']

        #vm = vm + 65e-3   ## gating parameter eqns assume resting is 0mV
        vm *= 1000.   ##  ..and that Vm is in mV
        
        an = (vm - 95) / (1.0 - np.exp(-(vm - 95) / 11.8))
        bn = 0.025 / np.exp(vm / 22.22)
        ntau = 1 / (an + bn)
        ninf = an * ntau
        dn = q10 * (ninf - n) / ntau
        return [dn*1e3]


class LGKslow(Channel):
    """Cortical slow potassium channel (Lewis & Gerstner 2002, p.124)
    """
    def __init__(self, gbar=0.225*mS/cm**2, **kwds):
        init_state = OrderedDict([('n', 0.0005)]) 
        Channel.__init__(self, gbar, init_state, **kwds)
        self.erev = -90*mV
        
    def open_probability(self, state):
        return state[self, 'n']**4

    def derivatives(self, state, t):
        # temperature dependence of rate constants
        # TODO: not sure about the base temp:
        q10 = 3 ** ((self.sim.temp - 37.) / 10.)
        
        vm = state[self.section, 'Vm']
        n = state[self, 'n']

        #vm = vm + 65e-3   ## gating parameter eqns assume resting is 0mV
        vm *= 1000.   ##  ..and that Vm is in mV
        
        an = 0.014 * (vm + 44) / (1.0 - np.exp(-(44 + vm) / 2.3))
        bn = 0.0043 / np.exp((vm + 44) / 34)
        ntau = 1 / (an + bn)
        ninf = an * ntau
        dn = q10 * (ninf - n) / ntau

        return [dn*1e3]

