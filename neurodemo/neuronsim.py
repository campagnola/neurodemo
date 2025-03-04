from collections import OrderedDict
import numpy as np
import scipy.integrate
import neurodemo.units as NU
import warnings


class Sim(object):
    """Simulator for a collection of objects that derive from SimObject"""

    def __init__(self, objects=None, temp=37.0, dt=10., integrator:str='solve_ivp'):
        if objects is None:
            objects = []
        self._objects = objects
        self._all_objs = None
        self._time = 0.0
        self.temp = temp
        self.dt = dt
        self.integrator = integrator

    def set_integrator(self, integrator:str):
        if integrator in ["odeint", "solve_ivp"]:
            self.integrator = integrator

    def change_dt(self, newdt:float=100.e-6):
        if newdt < 5e-6:
            newdt = 5e-6
        elif newdt > 1001e-6:
            newdt = 1000e-6
        self.dt = newdt * NU.s
    
    def add(self, obj):
        assert obj._sim is None
        obj._sim = self
        self._objects.append(obj)
        return obj

    def all_objects(self):
        """Ordered dictionary of all objects to be simulated, keyed by their names."""
        if self._all_objs is None:
            objs = OrderedDict()
            for o in self._objects:
                if not o.enabled:
                    continue
                for k, v in o.all_objects().items():
                    if k in objs:
                        raise NameError(
                            'Multiple objects with same name "%s": %s, %s'
                            % (k, objs[k], v)
                        )
                    objs[k] = v
            self._all_objs = objs
        return self._all_objs

    @property
    def time(self):
        return self._time

    def run(self, blocksize:int=1000, **kwds):
        """Run the simulation until a number of *samples* have been acquired.

        Extra keyword arguments are passed to `scipy.integrate.odeint()`.
        """
        # print("Integrator: ", self.integrator)
        # reset all_objs cache in case some part of the sim has changed
        self._all_objs = None
        all_objs = self.all_objects().values()

        # check that there is something to simulate
        if len(all_objs) == 0:
            raise RuntimeError("No objects added to simulation.")

        # Collect / prepare state variables for integration
        init_state = []
        difeq_vars = []
        dep_vars = {}
        for o in all_objs:
            pfx = o.name + "."
            for k, v in o.difeq_state().items():
                difeq_vars.append(pfx + k)
                init_state.append(v)
            for k, v in o.dep_state_vars.items():
                dep_vars[pfx + k] = v
        self._simstate = SimState(difeq_vars, dep_vars)
        t = np.arange(0, blocksize) * self.dt + self._time
        # print("\nstarting run at:", self._time)
        opts = {"rtol": 1e-6, "atol": 1e-8, "hmax": 5e-4, "full_output": 1}
        opts.update(kwds)
        # Run the simulation

        if self.integrator == 'odeint':
            result, info = scipy.integrate.odeint(self.derivatives, init_state, t, tfirst=True, **opts)
            p = 0
            for o in all_objs:
                nvar = len(o.difeq_state())
                # print("odeint state: ", p, nvar, result[-1, p:p+nvar])
                o.update_state(result[-1, p : p + nvar])
                p += nvar
            self._time = t[-1]
            # print(f"   {self.integrator:s}  final state = {str(result.T[:, -1]):s}")
            # print("   start, finished at : ", t[0],t[-1])
            # print("   np.min(result.T): ", np.min(result.T), np.max(result.T))
            return SimState(difeq_vars, dep_vars, result.T, integrator=self.integrator, t=t)

        elif self.integrator == 'solve_ivp':
            """Notes:
            Different integrators were tried. 
            LSODA works ok; RK's seem to stall on AP. 
            Radau and BDF are very fast, but take some large steps, so the calculation of where the 
            pulse array is in the derivative via get_cmd are not always correct -
            the algorithm might step a long time into the future, invalidating all the
            commands in the queue, which are removed. Thus, only the first cmd
            is executed.
            Probably should not pop the queue in get_cmd until we are certain at THIS
            level (or maybe in runner?) that the trigger arrays are actually finished. 
            """
        # in future we will need to implement this instead:
            result = scipy.integrate.solve_ivp(
                self.derivatives,
                t_span=(t[0], t[-1]),
                t_eval=t,
                y0=init_state,
                method="LSODA",  # runs ok with LSODA
 
                dense_output=False,
                # args=dep_vars,
                rtol = opts['rtol'], #**opts,
                atol = opts['atol'],
                max_step = opts['hmax'],
            )
            # Update current state variables
            p = 0
            for i, o in enumerate(all_objs):
                nvar = len(o.difeq_state())
                # print("solve ivp state: ", p, nvar, result.y[p:p+nvar, -1])
                o.update_state(result.y[p:p+nvar, -1])
                p += nvar
            self._time = t[-1]
            # print(f"\n   {self.integrator:s}  {str(result.y[:, -1]):s}")
            # print("   start, finished at : ", t[0],t[-1])
            # print("    np.min(result.y): ", np.min(result.y), np.max(result.y))
            return SimState(difeq_vars, dep_vars, result.y, integrator=self.integrator, t=t)

    def derivatives(self, t, state):
        objs = self.all_objects().values()
        self._simstate.state = state
        self._simstate.extra["t"] = t
        d = []
        for o in objs:
            d.extend(o.derivatives(self._simstate))
        return d

    def state(self):
        """Return dictionary of all dependent and independent state
        variables.
        """
        state = {}
        for o in self.all_objects():
            for k, v in o.state(self._simstate).items():
                state[k] = v
        return state

    def cmd_done(self):
        if 'soma.PatchClamp' in self.all_objects():
            return len(self.all_objects()['soma.PatchClamp'].cmd_queue) == 0
        return False


class MissingCurrentException(Exception):
    pass


class SimState(object):
    """Contains the state of all diff. eq. variables in the simulation.

    During simulation runs, this is used to carry information about all
    variables at the current timepoint. After the simulation finishes, this is
    used to carry all state variable data collected during the simulation.

    Parameters
    ==========
        difeq_vars: list
            Names of all diff. eq. state variables
        dep_vars: dict
            Name:function pairs for all dependent variables that may be computed
        difeq_state: list
            Initial values for all dif. eq. state variables
        extra:
            Extra name:value pairs that may be accessed from this object
    """

    def __init__(self, difeq_vars, dep_vars=None, difeq_state=None, integrator='odeint', **extra):
        self.difeq_vars = difeq_vars
        # record indexes of difeq vars for fast retrieval
        self.indexes = dict([(k, i) for i, k in enumerate(difeq_vars)])

        self.dep_vars = dep_vars
        self.state = difeq_state
        self.extra = extra
        self.integrator = integrator

    def set_state(self, difeq_state):
        self.state = difeq_state

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self.get_slice(key)
        # allow lookup by (object, var)
        if isinstance(key, tuple):
            key = key[0].name + "." + key[1]
        try:
            # try this first for speed
            return self.state[self.indexes[key]]
        except KeyError:
            if key in self.dep_vars:
                return self.dep_vars[key](self)
            elif key in self.extra.keys():
                # Key is usually 't' here
                return self.extra[key]
            else:
                raise MissingCurrentException

    def keys(self):
        return list(self.indexes.keys()) + list(self.dep_vars.keys()) + list(self.extra.keys())

    def __contains__(self, key):
        # allow lookup by (object, var)
        if isinstance(key, tuple):
            key = key[0].name + "." + key[1]
        return key in self.indexes or key in self.dep_vars or key in self.extra

    def __str__(self):
        rep = "SimState:\n"
        for i, k in enumerate(self.difeq_vars):
            rep += f"  {k} = {self.state[i][-1]}\n"
        return rep

    def get_final_state(self):
        """Return a dictionary of all diff. eq. state variables and dependent
        variables for all objects in the simulation.
        """
        return self.get_state_at_index(-1)

    def get_state_at_time(self, t):
        index = np.searchsorted(self['t'], t)
        return self.get_state_at_index(index)

    def get_state_at_index(self, index):
        state = {}
        s = self.copy()
        clip = not np.isscalar(self["t"])
        if clip:
            # only get results for the last timepoint
            # print('state: ', self.state)
            # if self.integrator == 'odeint':
            s.set_state(self.state[:, index])
            # else:
            #     s.set_state(self.state[:, -1])

        for k in self.difeq_vars:
            state[k] = s[k]
        for k in self.dep_vars:
            state[k] = s[k]
        for k, v in self.extra.items():
            if clip:
                state[k] = v[index]
            else:
                state[k] = v

        return state

    def get_slice(self, sl):
        kwds = {'difeq_state': self.state[:, sl]}
        for k,v in self.extra.items():
            kwds[k] = v[sl]
        return self.copy(**kwds)

    def copy(self, **kwds):
        default_kwds = {
            'difeq_vars': self.difeq_vars,
            'dep_vars': self.dep_vars,
            'difeq_state': self.state, 
            'integrator': self.integrator,
        }
        default_kwds.update(self.extra)
        default_kwds.update(kwds)
        return SimState(**default_kwds)


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
            if i == 0:
                name = self.type
            else:
                name = self.type + "%d" % i
        self._name = name
        self.enabled = True
        self._init_state = init_state.copy()  # in case we want to reset
        self._current_state = init_state.copy()
        self._sub_objs = []
        self.records = []
        self._rec_dtype = [(sv, float) for sv in init_state.keys()]

        # maps name:function for computing state vars that can be computed from
        # a SimState instance.
        self.dep_state_vars = {}

    @property
    def name(self):
        return self._name

    def all_objects(self):
        """SimObjects are organized in a hierarchy. This method returns an ordered
        dictionary of all enabled SimObjects in this branch of the hierarchy, beginning
        with self.
        """
        objs = OrderedDict()
        objs[self.name] = self
        for o in self._sub_objs:
            if not o.enabled:
                continue
            objs.update(o.all_objects())
        return objs

    def difeq_state(self):
        """An ordered dictionary of all variables required to solve the
        diff. eq. for this object.
        """
        return self._current_state

    def update_state(self, result):
        """Update diffeq state variables with their last simulated values.
        These will be used to initialize the solver when the next simulation
        begins.
        """
        for i, k in enumerate(self._current_state.keys()):
            # print("simobject: k, r: ", k, result[i])
            self._current_state[k] = result[i]

    def derivatives(self, state):
        """Return derivatives of all state variables.

        Must be reimplemented in subclasses. This is used by the ODE solver
        to integrate during the simulation; should be as fast as possible.
        """
        raise NotImplementedError()

    @property
    def sim(self):
        """The Sim instance in which this object is being used."""
        return self._sim


class Mechanism(SimObject):
    """Base class for simulation objects that interact with a section's
    membrane--channels, electrodes, etc.
    """

    def __init__(self, init_state, section=None, **kwds):
        SimObject.__init__(self, init_state, **kwds)
        self._name = kwds.pop("name", None)  # overwrite auto-generated name
        self._section = section
        self.dep_state_vars["I"] = self.current

    def current(self, state):
        """Return the membrane current being passed by this mechanism.

        Must be implemented in subclasses.
        """
        raise NotImplementedError()

    @property
    def name(self):
        if self._name is None:
            # pick a name that is unique to the section we live in

            # first collect all names
            names = []
            if self._section is None:
                return None
            for o in self._section.mechanisms:
                if isinstance(o, Mechanism) and o._name is None:
                    # skip to avoid recursion
                    continue
                names.append(o.name)

            # iterate until we find an unused name
            pfx = self._section.name + "."
            name = pfx + self.type
            i = 1
            while name in names:
                name = pfx + self.type + str(i)
                i += 1
            self._name = name
        return self._name

    @property
    def section(self):
        return self._section

    @property
    def sim(self):
        return self.section.sim


class Channel(Mechanism):
    """Base class for simple ion channels."""

    # precomputed rate constant tables
    rates = None

    # maximum open probability (to be redefined by subclasses)
    max_op = 1.0

    @classmethod
    def compute_rates(cls):
        return

    def __init__(self, gmax=None, gbar=None, init_state=None, **kwds):
        Mechanism.__init__(self, init_state, **kwds)
        self._gmax = gmax
        self._gbar = gbar

        if self.rates is None:
            type(self).compute_rates()
        self.dep_state_vars["G"] = self.conductance
        self.dep_state_vars["OP"] = self.open_probability

    @property
    def gmax(self):
        if self._gmax is not None:
            return self._gmax
        else:
            return self._gbar * self.section.area

    @gmax.setter
    def gmax(self, v):
        self._gmax = v
        self._gbar = None

    @property
    def gbar(self):
        if self._gbar is not None:
            return self._gbar
        else:
            return self._gmax / self.section.area

    @gbar.setter
    def gbar(self, v):
        self._gbar = v
        self._gmax = None

    def conductance(self, state):
        try:
            op = self.open_probability(state)
        except MissingCurrentException:
            return 0
        return self.gmax * op

    def current(self, state):
        vm = state[self.section, "V"]
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
            return rates[i1] * s + rates[i2] * (1 - s)


class Section(SimObject):
    type = "section"

    def __init__(self, radius=None, cap=10e-12 * NU.F, vm=-65 * NU.mV, **kwds):
        self.cap_bar = 1 * NU.uF / NU.cm**2
        if radius is None:
            self.cap = cap
            self.area = cap / self.cap_bar
        else:
            self.cap = self.area * self.cap_bar
            self.area = 4 * 3.1415926 * radius**2
        self.set_default_erev()
        init_state = OrderedDict([("V", vm)])
        SimObject.__init__(self, init_state, **kwds)
        self.dep_state_vars["I"] = self.current
        self.mechanisms = []

    def set_default_erev(self):
        self.ek = -77 * NU.mV
        self.ena = 50 * NU.mV
        self.ena1 = 74 * NU.mV
        self.ekf = -90 * NU.mV
        self.eks = -90 * NU.mV
        self.ecl = -70 * NU.mV
        self.eh = -43 * NU.mV
        self.eleak = -55 * NU.mV

    def add(self, mech):
        assert mech._section is None
        mech._section = self
        self.mechanisms.append(mech)
        self._sub_objs.append(mech)
        return mech

    def derivatives(self, state):
        Im = 0
        for mech in self.mechanisms:
            if not mech.enabled:
                continue
            Im += mech.current(state)

        dv = Im / self.cap
        return [dv]

    def current(self, state):
        """Return the current flowing across the membrane capacitance."""
        dv = self.derivatives(state)[0]
        return -self.cap * dv


class PatchClamp(Mechanism):
    type = "PatchClamp"

    def __init__(self, mode="ic", ra=0.1 * NU.MOhm, cpip=0.5e-12 * NU.F, **kwds):
        self.ra = ra
        self.cpip = cpip
        self._mode = mode
        self.cmd_queue = []
        self.cmd = []
        self.last_time = 0.0
        self.holding = {"ic": 0.0 * NU.pA, "vc": -65 * NU.mV}
        self.gain = 50e-6  # arbitrary VC gain
        init_state = OrderedDict([("V", -65 * NU.mV)])

        Mechanism.__init__(self, init_state, **kwds)
        self.dep_state_vars['cmd'] = self.get_cmd_from_state

    def queue_command(self, cmd, dt, start=None):
        """Execute a command as soon as possible.

        Return the time at which the command will begin.
        """
        assert cmd.ndim == 1 and cmd.shape[0] > 0
        if len(self.cmd_queue) == 0:
            next_start = self.last_time + dt
        else:
            last_start, last_dt, last_cmd = self.cmd_queue[-1]
            next_start = last_start + len(last_cmd) * last_dt

        if start is None:
            start = next_start
        else:
            if start < next_start:
                raise ValueError(
                    "Cannot start next command before %f; asked for %f."
                    % (next_start, start)
                )

        self.cmd_queue.append((start, dt, cmd))
        return start

    def queue_commands(self, cmds, dt):
        """Queue multiple commands for execution."""
        return [self.queue_command(c, dt) for c in cmds]

    @property
    def mode(self):
        return self._mode

    def clear_queue(self):
        self.cmd_queue = []

    def set_mode(self, mode):
        self._mode = mode
        self.clear_queue()

    def set_holding(self, mode, val):
        if mode not in self.holding:
            raise ValueError("Mode must be 'ic' or 'vc'")
        self.holding[mode] = val

    def current(self, state):
        # Compute current through tip of pipette at this timestep
        vm = state[self.section, "V"]
        ve = state[self, "V"]
        return (ve - vm) / self.ra

    def derivatives(self, state):
        t = state["t"]
        self.last_time = t
        ## Select between VC and CC
        cmd = self.get_cmd(t)

        # determine current generated by voltage clamp
        if self.mode == "vc":
            ve = state[self, "V"]
            cmd = (cmd - ve) * self.gain

        # Compute change in electrode potential
        dve = (cmd - self.current(state)) / self.cpip
        return [dve]

    def get_cmd_from_state(self, state):
        if isinstance(state['t'], np.ndarray):
            return [self.get_cmd(t) for t in state['t']]            
        else:
            return self.get_cmd(state['t'])

    def get_cmd(self, t: float):
        """Return command value at time *t*.

        Values are interpolated linearly between command points.
        """
        hold = self.holding[self.mode]
        while len(self.cmd_queue) > 0:
            (start, dt, data) = self.cmd_queue[0]
            i1 = int(np.floor((t - start) / dt))
            # if i1 >= 0:
            #     print(f"I!********************* {i1:8d}, {t:.4f}, {start:.3f}, {len(data):6f}")

            if i1 < -1:
                # before start of next command; return holding
                return hold
            elif i1 == -1:
                # interpolate from holding into start of next command
                v1 = hold
                vt1 = start - dt
                v2 = data[0]
                vt2 = start
                break
            elif i1 >= len(data):
                # this command has expired; remove and try next command
                self.cmd_queue.pop(0)  # NOTE: this does not work with some integration algorithms.
                continue
            else:
                v1 = data[i1]
                vt1 = start + i1 * dt
                if i1 + 1 < len(data):
                    # interpolate to next command point
                    v2 = data[i1 + 1]
                    vt2 = vt1 + dt
                    break
                else:
                    if len(self.cmd_queue) > 1 and vt1 + dt >= self.cmd_queue[1][0]:
                        # interpolate from command to next command array
                        v2 = self.cmd_queue[1][2][0]
                        vt2 = self.cmd_queue[1][0]
                    else:
                        # interpolate from command back to holding
                        v2 = hold
                        vt2 = vt1 + dt
                    break

        if len(self.cmd_queue) == 0:
            return hold

        s = (t - vt1) / (vt2 - vt1)
        return v1 * (1 - s) + v2 * s


class Leak(Channel):
    type = "Ileak"

    def __init__(self, gbar=0.1 * NU.mS / NU.cm**2, **kwds):
        Channel.__init__(self, gbar=gbar, init_state={}, **kwds)
    
    @property
    def erev(self):
        return self.section.eleak      
    
    def set_erev(self, erev):
        self.section.eleak = erev

    def open_probability(self, state):
        if state.state.ndim == 2:
            # need to return an array of the correct length..
            return np.ones(state.state.shape[1])
        else:
            return 1

    def derivatives(self, state):
        return []


class HHK(Channel):
    """Hodgkin-Huxley K channel."""

    type = "IK"

    max_op = 0.55

    @classmethod
    def compute_rates(cls):
        cls.rates_vmin = -100
        cls.rates_vstep = 0.1
        vm = np.arange(cls.rates_vmin, cls.rates_vmin + 400, cls.rates_vstep)
        cls.rates = np.empty((len(vm), 2))
        cls.rates[:, 0] = (0.1 - 0.01 * vm) / (np.exp(1.0 - 0.1 * vm) - 1.0)

        cls.rates[:, 1] = 0.125 * np.exp(-vm / 80.0)

    def __init__(self, gbar=12 * NU.mS / NU.cm**2, **kwds):
        init_state = OrderedDict([("n", 0.3)])
        Channel.__init__(self, gbar=gbar, init_state=init_state, **kwds)
        self.shift = 0
        self.lastn = init_state["n"]

    @property
    def erev(self):
        return self.section.ek

    def set_erev(self, erev):
        self.section.ek = erev

    def open_probability(self, state):
        return state[self, "n"] ** 4

    def check_state(self, state, gv, lastgv):
        n = state[self, gv]
        if np.isnan(n):
            n = lastgv
        if n > 1.0:
            n = 1.0
        elif n < 0:
            n = 0.
        return n, n

    def derivatives(self, state):
        # temperature dependence of rate constants
        q10 = 3 ** ((self.sim.temp - 6.3) / 10.0)
        vm = state[self.section, "V"] - self.shift

        vm = vm + 65e-3  ## gating parameter eqns assume resting is 0mV
        vm *= 1000.0  ##  ..and that Vm is in mV
        n = state[self, "n"]
        # n, self.lastn = self.check_state(state, "n", self.lastn)

        # disabled for now -- does not seem to improve speed.
        # an, bn = self.interpolate_rates(self.rates, vm, self.rates_vmin, self.rates_vstep)

        an = (0.1 - 0.01 * vm) / (np.exp(1.0 - 0.1 * vm) - 1.0)
        bn = 0.125 * np.exp(-vm / 80.0)
        dn = q10 * (an * (1.0 - n) - bn * n)
        return [dn * 1e3]


class HHNa(Channel):
    """Hodgkin-Huxley Na channel."""

    type = "INa"

    max_op = 0.2

    @classmethod
    def compute_rates(cls):
        cls.rates_vmin = -100
        cls.rates_vstep = 0.1
        vm = np.arange(cls.rates_vmin, cls.rates_vmin + 400, cls.rates_vstep)
        cls.rates = np.empty((len(vm), 4))
        cls.rates[:, 0] = (2.5 - 0.1 * vm) / (np.exp(2.5 - 0.1 * vm) - 1.0)
        cls.rates[:, 1] = 4.0 * np.exp(-vm / 18.0)
        cls.rates[:, 2] = 0.07 * np.exp(-vm / 20.0)
        cls.rates[:, 3] = 1.0 / (np.exp(3.0 - 0.1 * vm) + 1.0)

    def __init__(self, gbar=40 * NU.mS / NU.cm**2, **kwds):
        init_state = OrderedDict([("m", 0.05), ("h", 0.6)])
        Channel.__init__(self, gbar=gbar, init_state=init_state, **kwds)
        self.shift = 0
        self.lastm = init_state["m"]
        self.lasth = init_state["h"]

    @property
    def erev(self):
        return self.section.ena

    def set_erev(self, erev):
        self.section.ena = erev

    def open_probability(self, state):
        return state[self, "m"] ** 3 * state[self, "h"]

    def check_state(self, state, gv, lastgv):
        n = state[self, gv]
        if np.isnan(n):
            n = lastgv
        if n > 1.0:
            n = 1.0
        elif n < 0:
            n = 0.
        return n, n

    def derivatives(self, state):
        # temperature dependence of rate constants
        q10 = 3 ** ((self.sim.temp - 6.3) / 10.0)
        vm = state[self.section, "V"] - self.shift
        m = state[self, "m"]
        h = state[self, "h"]
        # m, self.lastm = self.check_state(state, "m", self.lastm)
        # h, self.lasth = self.check_state(state, "h", self.lasth) # state[self, "h"]

        vm = vm + 65e-3  ## gating parameter eqns assume resting is 0mV
        vm *= 1000.0  ##  ..and that Vm is in mV

        # disabled for now -- does not seem to improve speed.
        # am, bm, ah, bh = self.interpolate_rates(self.rates, vm, self.rates_vmin, self.rates_vstep)

        with warnings.catch_warnings(record=True) as w:
            am = (2.5 - 0.1 * vm) / (np.exp(2.5 - 0.1 * vm) - 1.0)
            bm = 4.0 * np.exp(-vm / 18.0)
            dm = q10 * (am * (1.0 - m) - bm * m)

            ah = 0.07 * np.exp(-vm / 20.0)
            bh = 1.0 / (np.exp(3.0 - 0.1 * vm) + 1.0)
            dh = q10 * (ah * (1.0 - h) - bh * h)
            # if abs(vm) > 1:
            #     print("Vm: ", vm)
            # if len(w) > 0:
            #     print(f"Vm: {vm:f}   expval: {2.5 - 0.1 * vm:f}")

        return [dm * 1e3, dh * 1e3]


class IH(Channel):
    """Ih from Destexhe 1993"""

    type = "IH"

    max_op = 0.3

    def __init__(self, gbar=30 * NU.mS / NU.cm**2, **kwds):
        init_state = OrderedDict([("f", 0), ("s", 0)])
        Channel.__init__(self, gbar=gbar, init_state=init_state, **kwds)
        # self.erev = -43 * NU.mV
        self.shift = 0
        self.lastf = init_state["f"]
        self.lasts = init_state["s"]

    @property
    def erev(self):
        return self.section.eh

    def set_erev(self, erev):
        self.section.eh = erev

    def open_probability(self, state):
        return state[self, "f"] * state[self, "s"]

    def check_state(self, state, gv, lastgv):
        n = state[self, gv]
        if np.isnan(n):
            n = lastgv
        if n > 1.0:
            n = 1.0
        elif n < 0:
            n = 0.
        return n, n

    def derivatives(self, state):
        vm = state[self.section, "V"] - self.shift
        # f, self.lastf = self.check_state(state, "f", self.lastf) # [self, "f"]
        # s, self.lasts = self.check_state(state, "s", self.lasts) # state[self, "s"]
        f = state[self, "f"]
        s = state[self, "s"]
        # vm = vm + 65e-3   ## gating parameter eqns assume resting is 0mV
        vm *= 1000.0  ##  ..and that Vm is in mV
        Hinf = 1.0 / (1.0 + np.exp((vm + 68.9) / 6.5))
        tauF = np.exp((vm + 158.6) / 11.2) / (1.0 + np.exp((vm + 75.0) / 5.5))
        tauS = np.exp((vm + 183.6) / 15.24)
        df = (Hinf - f) / tauF
        ds = (Hinf - s) / tauS
        return [df * 1e3, ds * 1e3]


class LGNa(Channel):
    """Cortical sodium channel: Ersir et al., 1999;
       Jolivet, Lewis and Gerstner, 2004
       For fast-spiking interneuron.
    """

    type = "INa1"

    def __init__(self, gbar=112.5 * NU.mS / NU.cm**2, **kwds):
        init_state = OrderedDict([("m", 0.019), ("h", 0.876)])
        Channel.__init__(self, gbar=gbar, init_state=init_state, **kwds)


    @property
    def erev(self):
        return self.section.ena1

    def set_erev(self, erev):
        self.section.ena1 = erev

    def open_probability(self, state):
        return state[self, "m"] ** 3 * state[self, "h"]

    def derivatives(self, state):
        # temperature dependence of rate constants
        # TODO: not sure about the base temp:
        q10 = 3 ** ((self.sim.temp - 37.0) / 10.0)

        vm = state[self.section, "V"]
        m = state[self, "m"]
        h = state[self, "h"]

        # vm = vm + 65e-3   ## gating parameter eqns assume resting is 0mV
        vm *= 1000.0  ##  ..and that Vm is in mV

        am = (-3020 + 40 * vm) / (1.0 - np.exp(-(vm - 75.5) / 13.5))
        bm = 1.2262 / np.exp(vm / 42.248)
        mtau = 1 / (am + bm)
        minf = am * mtau
        dm = q10 * (minf - m) / mtau

        ah = 0.0035 / np.exp(vm / 24.186)
        # note: bh as originally written causes integration failures; we use
        # an equivalent expression that behaves nicely under floating point stress.
        # bh = (0.8712 + 0.017 * vm) / (1.0 - np.exp(-(51.25 + vm) / 5.2))
        bh = 0.017 * (51.25 + vm) / (1.0 - np.exp(-(51.25 + vm) / 5.2))
        htau = 1.0 / (ah + bh)
        hinf = ah * htau
        dh = q10 * (hinf - h) / htau

        return [dm * 1e3, dh * 1e3]


class LGKfast(Channel):
    """Cortical potassium channel (fast): Ersir et al., 1999;
       Jolivet, Lewis and Gerstner, 2004
       For fast-spiking interneuron.
    """

    type = "IKf"

    def __init__(self, gbar=225 * NU.mS / NU.cm**2, **kwds):
        init_state = OrderedDict([("n", 0.00024)])
        Channel.__init__(self, gbar=gbar, init_state=init_state, **kwds)

    @property
    def erev(self):
        return self.section.ekf   

    def set_erev(self, erev):
        self.section.ekf = erev

    def open_probability(self, state):
        return state[self, "n"] ** 2

    def derivatives(self, state):
        # temperature dependence of rate constants
        # TODO: not sure about the base temp:
        q10 = 3 ** ((self.sim.temp - 37.0) / 10.0)

        vm = state[self.section, "V"]
        n = state[self, "n"]

        # vm = vm + 65e-3   ## gating parameter eqns assume resting is 0mV
        vm *= 1000.0  ##  ..and that Vm is in mV

        an = (vm - 95) / (1.0 - np.exp(-(vm - 95) / 11.8))
        bn = 0.025 / np.exp(vm / 22.22)
        ntau = 1 / (an + bn)
        ninf = an * ntau
        dn = q10 * (ninf - n) / ntau
        return [dn * 1e3]


class LGKslow(Channel):
    """Cortical potassium channel (slow): Ersir et al., 1999;
       Jolivet, Lewis and Gerstner, 2004
       For fast-spiking interneuron.
    """

    type = "IKs"

    def __init__(self, gbar=0.225 * NU.mS / NU.cm**2, **kwds):
        init_state = OrderedDict([("n", 0.0005)])
        Channel.__init__(self, gbar=gbar, init_state=init_state, **kwds)

    @property
    def erev(self):
        return self.section.eks   

    def set_erev(self, erev):
        self.section.eks = erev

    def open_probability(self, state):
        return state[self, "n"] ** 4

    def derivatives(self, state):
        # temperature dependence of rate constants
        # TODO: not sure about the base temp:
        q10 = 3 ** ((self.sim.temp - 37.0) / 10.0)

        vm = state[self.section, "V"]
        n = state[self, "n"]

        # vm = vm + 65e-3   ## gating parameter eqns assume resting is 0mV
        vm *= 1000.0  ##  ..and that Vm is in mV

        an = 0.014 * (vm + 44) / (1.0 - np.exp(-(44 + vm) / 2.3))
        bn = 0.0043 / np.exp((vm + 44) / 34)
        ntau = 1 / (an + bn)
        ninf = an * ntau
        dn = q10 * (ninf - n) / ntau

        return [dn * 1e3]


# alpha synapse
Area = 1e-12 * NU.cm**2
Alpha_t0 = 500.0 * NU.ms  # msec
Alpha_tau = 2.0 * NU.ms
gAlpha = 1e-3 * Area / NU.cm**2
EAlpha = -7 * NU.mV  # mV


def IAlpha(Vm, t):
    if t < Alpha_t0:
        return 0.0
    else:
        # g = gmax * (t - onset)/tau * exp(-(t - onset - tau)/tau)
        tn = t - Alpha_t0
        if tn > 10.0 * Alpha_tau:
            return 0.0
        else:
            return (
                gAlpha
                * (Vm - EAlpha)
                * (tn / Alpha_tau)
                * np.exp(-(tn - Alpha_tau) / Alpha_tau)
            )
