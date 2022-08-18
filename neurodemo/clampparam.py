# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""

import numpy as np
from pyqtgraph.Qt import QtCore
import pyqtgraph.parametertree as pt
from .sequenceplot import SequencePlotWindow
import neurodemo.units as NU


class ClampParameter(pt.parameterTypes.SimpleParameter):
    # emitted when a plot should be shown or hidden
    plots_changed = QtCore.Signal(
        object, object, object, object
    )  # self, channel, name, on/off

    def __init__(self, clamp, sim):
        self.clamp = clamp
        self.sim = sim
        self.dt = sim.dt
        self.plot_win = SequencePlotWindow()

        self.triggers = (
            []
        )  # items are (trigger_time, pointer, trigger_buffer, (mode, amp, cmd, seq_ind, seq_len))
        self.plot_keys = []
        pt.parameterTypes.SimpleParameter.__init__(
            self,
            name="Patch Clamp",
            type="bool",
            value=True,
            children=[
                dict(
                    name="Mode",
                    type="list",
                    values={"CC": "ic", "VC": "vc"},
                    value="ic",
                ),
                dict(
                    name="Holding",
                    type="float",
                    value=0,
                    suffix="A",
                    siPrefix=True,
                    step=10 * NU.pA,
                ),
                # dict(name='Ideal', type='bool', value=True),
                dict(
                    #name="Pipette Capacitance",
                    name="Pipette Cap",
                    type="float",
                    value=clamp.cpip,
                    limits=[0.01 * NU.pF, None],
                    suffix="F",
                    siPrefix=True,
                    dec=True,
                    step=0.5,
                ),
                dict(
                    #name="Access Resistance",
                    name="Access Res",
                    type="float",
                    value=clamp.ra,
                    limits=[10 * NU.kOhm, None],
                    suffix="Ω",
                    siPrefix=True,
                    step=0.5,
                    dec=True,
                ),
                dict(name="Plot Current", type="bool", value=False),
                dict(name="Plot Voltage", type="bool", value=False),
                dict(
                    name="Pulse",
                    type="group",
                    children=[
                        dict(name="Capture Results", type="bool", value=False),
                        dict(name="Pulse Once", type="action"),
                        dict(
                            name="Hold-duration",
                            type="float",
                            value=10 * NU.ms,
                            suffix="s",
                            siPrefix=True,
                            limits=[0, None],
                            step=5e-3,
                        ),
                        dict(
                            name="Pre-duration",
                            type="float",
                            value=10 * NU.ms,
                            suffix="s",
                            siPrefix=True,
                            limits=[0, None],
                            step=5e-3,
                        ),
                        dict(
                            name="Pre-amplitude",
                            type="float",
                            value=0 * NU.pA,
                            suffix="A",
                            siPrefix=True,
                            dec=True,
                        ),
                        dict(
                            name="Duration",
                            type="float",
                            value=50 * NU.ms,
                            suffix="s",
                            siPrefix=True,
                            limits=[0, None],
                            step=5e-3,
                        ),
                        dict(
                            name="Amplitude",
                            type="float",
                            value=50 * NU.pA,
                            suffix="A",
                            siPrefix=True,
                            dec=True,
                        ),
                        dict(
                            name="Post-duration",
                            type="float",
                            value=50 * NU.ms,
                            suffix="s",
                            siPrefix=True,
                            limits=[0, None],
                            step=5e-3,
                        ),
                        dict(
                            name="Post-amplitude",
                            type="float",
                            value=0 * NU.pA,
                            suffix="A",
                            siPrefix=True,
                            dec=True,
                        ),
                        dict(name="Pulse Sequence", type="action"),
                        dict(name="Sequence Pulse", 
                            type="list",
                            values={"Pre": 1, "Pulse": 2, "Post": 3},
                            value=2,
                            ),
                        dict(
                            name="Start Amplitude",
                            type="float",
                            value=-50 * NU.pA,
                            suffix="A",
                            siPrefix=True,
                            dec=True,
                        ),
                        dict(
                            name="Stop Amplitude",
                            type="float",
                            value=50 * NU.pA,
                            suffix="A",
                            siPrefix=True,
                            dec=True,
                        ),
                        dict(name="Post-holdtime", type="float",
                            value = 500*NU.ms,
                            suffix="s",
                            si_refix=True,
                            limits=[0.01, 2.0],
                        ),

                        dict(
                            name="Pulse Number", type="int", value=11, limits=[2, None]
                        ),
                    ],
                ),
            ],
        )
        self.sigTreeStateChanged.connect(self.treeChange)
        self.child("Pulse", "Pulse Once").sigActivated.connect(self.pulse_once)
        self.child("Pulse", "Pulse Sequence").sigActivated.connect(self.pulse_sequence)

    def treeChange(self, root, changes):
        for param, change, val in changes:
            if change != "value":
                continue
            if param is self:
                self.clamp.enabled = val
            elif param is self.child("Mode"):
                self.set_mode(val)
            elif param is self.child("Holding"):
                self.clamp.set_holding(self.mode(), val)
            elif param is self.child("Pipette Capacitance"):
                self.clamp.cpip = val
            elif param is self.child("Access Resistance"):
                self.clamp.ra = val
            elif param is self.child("Plot Current"):
                self.plots_changed.emit(self, self.clamp, "I", val)
            elif param is self.child("Plot Voltage"):
                self.plots_changed.emit(self, self.clamp, "V", val)
            # elif param is self.child('Ideal'):
            # self.child('Pipette Capacitance').setOpts(visible=not val)
            # self.child('Access Resistance').setOpts(visible=not val)
            # self.clamp.set_ideal(val)

    def mode(self):
        return self["Mode"]

    def set_mode(self, mode):
        self.clamp.set_mode(mode)
        suff = {"ic": "A", "vc": "V"}[mode]
        pre_amp, amp, post_amp, start, stop, step = {
            "ic": (0 * NU.pA, -10 * NU.pA, 0 * NU.pA,-100 * NU.pA, 100 * NU.pA, 10 * NU.pA),
            "vc": (0 * NU.mV, -10 * NU.mV, 0 * NU.mV, -40 * NU.mV, 100 * NU.mV, 5 * NU.mV),
        }[mode]
        self.sigTreeStateChanged.disconnect(self.treeChange)
        try:
            self.child("Holding").setOpts(
                suffix=suff, value=self.clamp.holding[mode], step=step
            )
            self.child("Pulse", "Pre-amplitude").setOpts(suffix=suff, value=pre_amp, step=step)
            self.child("Pulse", "Post-amplitude").setOpts(suffix=suff, value=post_amp, step=step)
            self.child("Pulse", "Amplitude").setOpts(suffix=suff, value=amp, step=step)
            self.child("Pulse", "Start Amplitude").setOpts(
                suffix=suff, value=start, step=step
            )
            self.child("Pulse", "Stop Amplitude").setOpts(
                suffix=suff, value=stop, step=step
            )
        finally:
            self.sigTreeStateChanged.connect(self.treeChange)

    def pulse_template(self):
        d0 = self["Pulse", "Hold-duration"]
        d1 = self["Pulse", "Pre-duration"]
        d2 = self["Pulse", "Duration"]
        d3 = self["Pulse", "Post-duration"]
        d4 = self["Pulse", "Post-holdtime"]
        dur = d0 + d1 + d2 + d3 + d4
        durs = [d0, d1, d2, d3, d4]
        idurs = [0]*len(durs)
        npts = int(dur / self.dt)
        # cmd = np.zeros(npts)
        it0 = 0
        for i, d in enumerate(durs):
            idurs[i] = it0 + int(d / self.dt)
            it0 += idurs[i]
        # i1 = i0 + int(d1 / self.dt)
        # i2 = i1 + int(d2 / self.dt)
        # i3 = i2 + int(d3 / self.dt)
        # i4 = i3 + int(d4 / self.dt)
        cmd = np.ones(npts)*self["Holding"]
        return cmd, idurs # i0, i1, i2, i3, i4

    def pulse_once(self):
        cmd, idurs = self.pulse_template()
        amp_pre = self["Pulse", "Pre-amplitude"]
        cmd[idurs[0]:idurs[1]] += amp_pre
        amp = self["Pulse", "Amplitude"]
        cmd[idurs[1]:idurs[2]] += amp
        amp_post = self["Pulse", "Post-amplitude"]
        cmd[idurs[2]:idurs[3]] += amp_post
        t = self.clamp.queue_command(cmd, self.dt)
        if self["Pulse", "Capture Results"]:
            info = {
                "mode": self.mode(),
                "amp": amp,
                "cmd": cmd,
                "seq_ind": 0,
                "seq_len": 0,
            }
            self.add_trigger(len(cmd), t, info)

    def pulse_sequence(self):
        cmd, idurs = self.pulse_template()
        cmds = []
        amps = np.linspace(
            self["Pulse", "Start Amplitude"],
            self["Pulse", "Stop Amplitude"],
            self["Pulse", "Pulse Number"],
        )
        for amp in amps:
            cmd2 = cmd.copy()
            if self['Pulse', "Sequence Pulse"] == 1:
                cmd2[idurs[0]:idurs[1]] += amp
                cmd2[idurs[1]:idurs[2]] += self["Pulse", "Amplitude"]
                cmd2[idurs[2]:idurs[3]] += self["Pulse", "Post-amplitude"]
            elif self['Pulse', "Sequence Pulse"] == 2:
                cmd2[idurs[0]:idurs[1]] += self["Pulse", "Pre-amplitude"]
                cmd2[idurs[1]:idurs[2]] += amp
                cmd2[idurs[2]:idurs[3]] += self["Pulse", "Post-amplitude"]
 
            elif self['Pulse', "Sequence Pulse"] == 3:
                cmd2[idurs[0]:idurs[1]] += self["Pulse", "Pre-amplitude"]
                cmd2[idurs[1]:idurs[2]] += self["Pulse", "Amplitude"]
                cmd2[idurs[2]:idurs[3]] += amp
            cmds.append(cmd2)

        times = self.clamp.queue_commands(cmds, self.dt)
        for i, t in enumerate(times):
            info = {
                "mode": self.mode(),
                "amp": amps[i],
                "cmd": cmds[i],
                "seq_ind": i,
                "seq_len": len(amps),
            }
            self.add_trigger(len(cmd), t, info)

    def add_trigger(self, n, t, info):
        buf = np.empty(n, dtype=[(str(k), float) for k in self.plot_keys + ["t"]])
        self.triggers.append([t, 0, buf, info])

    def add_plot(self, key, label):
        self.plot_keys.append(key)
        self.triggers = []
        self.plot_win.add_plot(key, label)

    def remove_plot(self, key):
        self.plot_keys.remove(key)
        self.triggers = []
        self.plot_win.remove_plot(key)

    def new_result(self, result):
        if len(self.triggers) == 0:
            return
        t = result["t"]
        tt, ptr, buf, info = self.triggers[0]
        if tt > t[-1]:
            # no triggers ready
            return

        # Copy data from result to trigger buffer
        i = max(
            0, int(np.round((tt - t[0]) / self.dt))
        )  # index of trigger within new data
        npts = min(
            len(buf) - ptr, len(t) - i
        )  # number of samples to copy from new data
        for k in self.plot_keys:  # self.plot_keys is a list, buf is a recarray
            buf[k][ptr : ptr + npts] = result[k][i : i + npts]

        ptr += npts
        if ptr >= buf.shape[0]:
            # If the trigger buffer is full, plot and remove
            self.plot_win.plot(np.arange(buf.shape[0]) * self.dt, buf, info)
            self.triggers.pop(0)
            if len(t) > npts:
                # If there is data left over, try feeding it to the next trigger
                result = dict([(k, result[k][i + npts :]) for k in result])
                self.new_result(result)
        else:
            # otherwise, update the pointer and wait for the next result
            self.triggers[0][1] = ptr
