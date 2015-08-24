# -*- coding: utf-8 -*-
"""
Simple Hodgkin-Huxley simulator for Python. VERY slow.
Includes Ih from Destexhe 1993 [disabled]
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



Raccess = 10*MΩ
Cpip = 3*pF

Temp = 6.3

Radius = 20 * um
Area = (4 * 3.1415926 * Radius**2)
C = Area * uF / cm**2
gNa = 40e-3 * Area/cm**2
gK =  12e-3 * Area/cm**2
gKShift = 0
gL =  0.1e-3 * Area/cm**2
gH = 30. * Area/cm**2  ## 10mS * cm^2

EK = -77e-3
ENa = 50e-3
EL = -55e-3
EH = -43e-3
Ik2 = 24.3

# alpha synapse
Alpha_t0 = 500.  # msec
Alpha_tau = 2.0
gAlpha = 1e-3 * Area/cm**2
EAlpha = -7e-3  # V

def IK(n, Vm):
    return gK * n**4 * (Vm - EK)
    
def INa(m, h, Vm):
    return gNa * m**3 * h * (Vm - ENa)

def IL(Vm):
    return gL * (Vm - EL)

def IH(Vm, f, s):
    return gH * f * s * (Vm - EH)

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
        self._time = 0.0

    def add(self, obj):
        self._objects.append(obj)

    def all_objects(self):
        objs = []
        for o in self._objects:
            objs.extend(o.all_objects())
        return objs
            
    def state(self):
        s = OrderedDict()
        for o in self.all_objects():
            for k,v in o.state().items():
                s[k] = v
        return s
    
    @property
    def state_vars(self):
        return list(self.state().keys())

    def run(self, dt=0.1, dur=100, **args):
        npts = int(dur/dt)
        t = np.linspace(0, dur, npts)
        result = np.empty((npts, len(self.state_vars) + 1))

        # Run the simulation
        init_state = list(self.state().values())
        (result[:,1:], info) = scipy.integrate.odeint(self.derivatives, init_state, t, (dt,),
                                                rtol=1e-6, atol=1e-6, hmax=5e-2, full_output=1, **args)
        
        p = 0
        for o in self.all_objects():
            nvar = len(o.state_vars)
            o._last_result = result[:, p:p+nvar]
            p += nvar
            
        self._last_run_time = t
        #t, Im, Ve, Vm, m, h, n, f, s = [result[:,i] for i in range(result.shape[1])]
        
        # Compute electrode current sans pipette capacitance current
        #result[:,1] = (Ve-Vm) / Raccess

        # update state so next run starts where we left off
        #self.state = result[-1, 2:]

    def derivatives(self, state, t, dt):
        d = []
        for o in self.all_objects():
            nvars = len(o.state_vars)
            substate = state[:nvars]
            state = state[nvars:]
            d.extend(o._get_derivs(substate, t, dt))
        return d
    


class SimObject(object):
    """
    Base class for objects that participate in integration by providing a set
    of state variables and their derivatives.
    """
    def __init__(self, init_state):
        self._init_state = init_state
        self._last_state = init_state.copy()
        self._sub_objs = []
        self.records = []
        self._rec_dtype = [(sv, float) for sv in init_state.keys()]
    
    def state(self):
        return self._last_state.copy()

    def all_objects(self):
        objs = [self]
        for o in self._sub_objs:
            objs.extend(o.all_objects())
        return objs
        
    @property
    def state_vars(self):
        return list(self._init_state.keys())

    def _get_derivs(self, state, t, dt):
        for i,k in enumerate(self._last_state):
            self._last_state[k] = state[i]
        return self.derivatives(state, t, dt)
    
    def derivatives(self, state, t, dt):
        """Return derivatives of all state variables.
        
        Must be reimplemented in subclasses.
        """
        raise NotImplementedError()


class Mechanism(SimObject):
    def __init__(self, init_state, section=None):
        SimObject.__init__(self, init_state)
        self._section = section
        
    def current(self):
        """Return the membrane current being passed by this mechanism.
        
        Must be implemented in subclasses.
        """
        raise NotImplementedError()

    @property
    def vm(self):
        return self._section._last_state['Vm']


class Section(SimObject):
    def __init__(self):
        #self.state = [-65e-3, -65e-3, 0.05, 0.6, 0.3, 0.0, 0.0]
        init_state = OrderedDict([('Vm', -65*mV)])
        SimObject.__init__(self, init_state)
        self.mechanisms = []

    def add(self, mech):
        assert mech._section is None
        mech._section = self
        self.mechanisms.append(mech)
        self._sub_objs.append(mech)

    def derivatives(self, state, t, dt):
        Im = 0
        for mech in self.mechanisms:
            Im += mech.current()
            
        dv = 1e-3 * Im / C    # 1e-3 is because t is expressed in ms
        return [dv]
        
    def old_derivs(self, y, t, mode, cmd, dt):
        ## y is a vector [Ve, Vm, m, h, n, f, s], function returns derivatives of each variable
        ## with respect to time.
        ## t is current time expressed in ms.
        ## mode is 'ic' or 'vc'
        ## cmd is an array of ic or vc command values
        ## dt gives the time step of the command data.
        (Ve, Vm, m, h, n, f, s) = y
        ## Select between VC and CC
        if cmd is None:
            cmd = 0
            mode = 'ic'
        else:
            # interpolate command -- sharp steps confuse the integrator.
            fInd = t/dt
            ind = min(len(cmd)-1, np.floor(fInd))
            ind2 = min(len(cmd)-1, ind+1)
            s = fInd - ind
            cmd = cmd[ind] * (1-s) + cmd[ind2] * s
            
        # determine current generated by voltage clamp 
        if mode == 'vc':
            G = 50e-6 # arbitrary VC gain
            cmd = (cmd-Ve) * G
            
        # limit current
        #cmd = np.clip(cmd, -2*nA, 2*nA)
            
        # Compute current through tip of pipette
        Iaccess = (Ve-Vm) / Raccess
        
        # Compute change in electrode potential
        dve = 1e-3 * (cmd - Iaccess) / Cpip    # 1e-3 is because t is expressed in ms
        
        # Compute change in membrane potential
        #Im = (Iaccess - INa(m, h, Vm) - IK(n, Vm) - IL(Vm) - IH(Vm, f, s) - IAlpha(Vm, t))
        Im = Iaccess - INa(m, h, Vm) - IK(n, Vm) - IL(Vm) 
        dv = 1e-3 * Im / C    # 1e-3 is because t is expressed in ms
        
        ## Compute changes in gating parameters
        Vm = Vm + 65e-3   ## gating parameter eqns assume resting is 0mV
        Vm *= 1000.   ##  ..and that Vm is in mV
        
        # temperature dependence of rate constants
        q10 = 3 ** ((Temp-6.3) / 10.)
        
        am = (2.5-0.1*Vm) / (np.exp(2.5-0.1*Vm) - 1.0)
        bm = 4. * np.exp(-Vm / 18.)
        dm = q10 * (am * (1.0 - m) - bm * m)
        #am = 0.1 * (-Vm-40) / (np.exp((-Vm-40) / 10.) - 1.0)
        #bm = 4.0 * np.exp(-(Vm+65.0) / 18.0)
        #mtau = 1. / (q10 * (am + bm))
        #minf = am / (am + bm)
        #dm = (minf - m) / mtau
        
        ah = 0.07 * np.exp(-Vm / 20.)
        bh = 1.0 / (np.exp(3.0 - 0.1 * Vm) + 1.0)
        dh = q10 * (ah * (1.0 - h) - bh * h)
        #ah = .07 * np.exp(-(Vm+65.) / 20.)
        #bh = 1. / (np.exp(-(Vm+35.) / 10.) + 1.)
        #htau = 1. / (q10 * (ah + bh))
        #hinf = ah / (ah + bh)
        #dh = (hinf - h) / htau
        
        an = (0.1 - 0.01*(Vm-gKShift)) / (np.exp(1.0 - 0.1*(Vm-gKShift)) - 1.0)
        bn = 0.125 * np.exp(-Vm / 80.)
        dn = q10 * (an * (1.0 - n) - bn * n)
        #an = .01 * (-Vm-55) / (np.exp((-Vm-55) / 10.) - 1.0)
        #bn = .125 * np.exp(-(Vm+65.) / 80.)
        #ntau = 1. / (q10 * (an + bn))
        #ninf = an / (an + bn)
        #dn = (ninf - n) / ntau
        
        # Ih is disabled--very slow.
        #Hinf = 1.0 / (1.0 + np.exp((Vm + 68.9) / 6.5))
        #tauF = np.exp((Vm + 158.6)/11.2) / (1.0 + np.exp((Vm + 75.)/5.5))
        #tauS = np.exp((Vm + 183.6) / 15.24)
        #df = (Hinf - f) / tauF
        #ds = (Hinf - s) / tauS
        df = 0
        ds = 0
        
        return [dve, dv, dm, dh, dn, df, ds]

    def old_run(self, mode='ic', cmd=None, dt=0.1, dur=100, **args):
        npts = int(dur/dt)
        t = np.linspace(0, dur, npts)
        result = np.empty((npts, 9))

        # Run the simulation
        (result[:,2:], info) = scipy.integrate.odeint(self.derivatives, self.state, t, (mode, cmd, dt),
                                                rtol=1e-6, atol=1e-6, hmax=5e-2, full_output=1, **args)
        result[:,0] = t
        t, Im, Ve, Vm, m, h, n, f, s = [result[:,i] for i in range(result.shape[1])]
        
        # Compute electrode current sans pipette capacitance current
        result[:,1] = (Ve-Vm) / Raccess

        # update state so next run starts where we left off
        #self.state = result[-1, 2:]
        
        
        
        return result  ## result is array with dims: [npts, (time, Ie, Ve, Vm, Im, m, h, n, f, s)]


class Leak(Mechanism):
    def __init__(self, g, erev=-55*mV):
        self.g = g
        self.erev = erev
        Mechanism.__init__(self, {})
        
    def current(self):
        return self.g * (self.vm - self.erev)

    def derivatives(self, state, t, dt):
        return []


class MultiClamp(Mechanism):
    def __init__(self, mode='ic', cmd=None, ra=5*MΩ, cpip=3*pF):
        self.ra = ra
        self.cpip = cpip
        self.mode = mode
        self.cmd = cmd
        self.gain = 50e-6  # arbitrary VC gain
        init_state = OrderedDict([('Ve', -65*mV)])
        Mechanism.__init__(self, init_state)

    def set_command(self, cmd):
        self.cmd = cmd

    def current(self):
        # Compute current through tip of pipette
        return (self.state()['Ve'] - self.vm) / self.ra

    def derivatives(self, state, t, dt):
        ## Select between VC and CC
        cmd = self.cmd
        if cmd is None:
            cmd = 0
            mode = 'ic'
        else:
            # interpolate command -- sharp steps confuse the integrator.
            fInd = t/dt
            ind = min(len(cmd)-1, np.floor(fInd))
            ind2 = min(len(cmd)-1, ind+1)
            s = fInd - ind
            cmd = cmd[ind] * (1-s) + cmd[ind2] * s
            
        # determine current generated by voltage clamp 
        if self.mode == 'vc':
            ve = self.state()['Ve']
            cmd = (cmd-ve) * self.gain
        
        # Compute change in electrode potential
        dve = 1e-3 * (cmd - self.current()) / self.cpip    # 1e-3 is because t is expressed in ms
        return [dve]


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
    neuron.add(Leak(g=0.1e-3 * Area/cm**2))
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
    cmd = np.zeros((len(x), npts)) #*-65e-3
    data = np.zeros((len(x), npts, 9))
    for i, v in enumerate(x):
        print('V: ', v)
        cmd[i, x1:x2] = v
        clamp.set_command(cmd[i])
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
    if sys.flags.interactive == 0:
        QtGui.QApplication.instance().exec_()
