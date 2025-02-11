from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    # TYPE_CHECKING is always false at runtime, preventing unnecessary runtime import of modules that
    # are only needed by IDE for type checking
    import neurodemo.main_window

from dataclasses import dataclass
import numpy as np

from . import qt
import pyqtgraph.parametertree as pt
from .sequenceplot import SequencePlotWindow
import neurodemo.units as NU

@dataclass
class Trigger:
    trigger_time: float
    curr_buff_ptr: int
    buf: object
    info: dict


class ClampParameter(pt.parameterTypes.SimpleParameter):
    # emitted when a plot should be shown or hidden
    plots_changed = qt.Signal(
        object, object, object, object
    )  # self, channel, name, on/off
    mode_changed = qt.Signal(object, object)  # self, mode

    def __init__(self, clamp: neurodemo.neuronsim.PatchClamp, sim: neurodemo.main_window.DemoWindow):
        self.clamp = clamp
        self.sim = sim
        self.dt = sim.dt
        self.dt_updated = True
        self.plot_win = SequencePlotWindow()

        self.triggers = []  # items are (trigger_time, pointer, trigger_buffer, (mode, amp, cmd, seq_ind, seq_len))
        self.result_buffer = []  # store a few recent results to ensure triggers are caught
        self.result_buffer_size = 5

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
                    limits={"CC": "ic", "VC": "vc"},
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
                    limits=[0.1*1e-12 * NU.F, 1e-9 * NU.F],
                    suffix="F",
                    siPrefix=True,
                    dec=True,
                    step=0.5e-12,
                ),
                dict(
                    #name="Access Resistance",
                    name="Access Res",
                    type="float",
                    value=clamp.ra,
                    limits=[10000 * NU.Ohm, None],
                    suffix="Ω",
                    siPrefix=True,
                    step=0.5,
                    dec=True,
                ),
                dict(name="Plot Current", type="bool", value=False),
                dict(name="Plot Voltage", type="bool", value=False),
                dict(name="Plot Command", type="bool", value=False),
                dict(
                    name="Pulse",
                    type="group",
                    children=[
                        dict(name="Capture Results", type="bool", value=False),
                        dict(name="Clear Pulses", type="action"),
                        dict(name="Pulse Once", type="action"),
                        dict(
                            name="Hold-duration",
                            type="float",
                            value=0.010 * NU.s,
                            suffix="s",
                            siPrefix=True,
                            limits=[0, None],
                            step=5e-3,
                        ),
                        dict(
                            name="Pre-duration",
                            type="float",
                            value=0.010 * NU.s,
                            suffix="s",
                            siPrefix=True,
                            limits=[0, None],
                            step=5e-3,
                        ),
                        dict(
                            name="Pre-amplitude",
                            type="float",
                            value=0 * NU.A,
                            suffix="A",
                            siPrefix=True,
                            dec=True,
                        ),
                        dict(
                            name="Duration",
                            type="float",
                            value=0.050 * NU.s,
                            suffix="s",
                            siPrefix=True,
                            limits=[0, None],
                            step=5e-3,
                        ),
                        dict(
                            name="Amplitude",
                            type="float",
                            value=50e-12 * NU.A,
                            suffix="A",
                            siPrefix=True,
                            dec=True,
                        ),
                        dict(
                            name="Post-duration",
                            type="float",
                            value=0.050 * NU.s,
                            suffix="s",
                            siPrefix=True,
                            limits=[0, None],
                            step=5e-3,
                        ),
                        dict(
                            name="Post-amplitude",
                            type="float",
                            value=0 * NU.A,
                            suffix="A",
                            siPrefix=True,
                            dec=True,
                        ),
                        dict(name="Pulse Sequence", type="action"),
                        dict(name="Sequence Pulse", 
                            type="list",
                            limits={"Pre": 1, "Pulse": 2, "Post": 3},
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
                            step=1.0,
                            minStep=1 * NU.pA,
                        ),
                        dict(name="Post-holdtime", type="float",
                            value = 50 * NU.ms,
                            suffix="s",
                            siPrefix=True,
                            limits=[0.001, 5.0],
                            dec=True,
                            step=1.0,
                            minStep=0.001,
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
        self.child("Pulse", "Clear Pulses").sigActivated.connect(self.clear_triggers)

    def set_dt(self, dt):
        self.dt = dt
        self.dt_updated = True
    
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
            elif param is self.child("Pipette Cap"):
                self.clamp.cpip = val
            elif param is self.child("Access Res"):
                self.clamp.ra = val
            elif param is self.child("Plot Current"):
                self.plots_changed.emit(self, self.clamp, "I", val)
            elif param is self.child("Plot Voltage"):
                self.plots_changed.emit(self, self.clamp, "V", val)
            elif param is self.child("Plot Command"):
                self.plots_changed.emit(self, self.clamp, "cmd", val)
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
        self.mode_changed.emit(self, mode)

    def pulse_template(self):
        d0 = self["Pulse", "Hold-duration"]
        d1 = self["Pulse", "Pre-duration"]
        d2 = self["Pulse", "Duration"]
        d3 = self["Pulse", "Post-duration"]
        d4 = self["Pulse", "Post-holdtime"]
        dur = d0 + d1 + d2 + d3 + d4
        durs = [d0, d1, d2, d3, d4]
        idurs = [0]*len(durs)  # indexes
        npts = int(dur / self.dt)  # points in stimulus
        # print(dur, self.dt)
        # print("template: npts= ", npts)
        it0 = 0
        t0 = 0.
        for i, d in enumerate(durs):
            idurs[i] = it0 + int(d / self.dt)
            it0 = idurs[i]
            t0 = durs[i]
        cmd = np.ones(npts)*self["Holding"]
        return cmd, idurs # i0, i1, i2, i3, i4

    def print_triggers(self):
        if len(self.triggers) == 0:
            print("No triggers set")
            return
        print("Triggers:\n")
        for tr in self.triggers:
            print(f"    {tr.trigger_time:9.5f}  {tr.info['amp']*1e12:8.1f} pA")

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

        if not self.sim.running():
            # If not already running, then start a time-limited run.
            self.sim.start(stop_after_cmd=True)

        # self.print_triggers()

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
        # note: it is possible for data to arrive before triggers have been added.
        for i, t in enumerate(times):
            info = {
                "mode": self.mode(),
                "amp": amps[i],
                "cmd": cmds[i],
                "seq_ind": i,
                "seq_len": len(amps),
            }
            self.add_trigger(len(cmd), t, info)
        # self.print_triggers()
        if not self.sim.running():
            # If not already running, then start a time-limited run
            self.sim.start(stop_after_cmd=True)


    def add_trigger(self, n, t, info):
        buf = np.empty(n, dtype=[(str(k), float) for k in self.plot_keys + ["t"]])
        self.triggers.append(Trigger(t, 0, buf, info))

    def clear_triggers(self):
        self.triggers = []
        self.clamp.clear_queue()
        
    def add_plot(self, key, label):
        self.plot_keys.append(key)
        self.triggers = []
        self.plot_win.add_plot(key, label)

    def remove_plot(self, key):
        self.plot_keys.remove(key)
        self.triggers = []
        self.plot_win.remove_plot(key)

    def get_oldest_result(self):
        # Find index of result_buffer having the lowest timestamp
        earliest_t = -1
        earliest_i = -1
        for i, r in enumerate(self.result_buffer):
            t = r["t"][0]
            if earliest_t == -1:
                earliest_t = t
                earliest_i = i
            else:
                if t < earliest_t:
                    earliest_t = t
                    earliest_i = i
        if earliest_i >= 0:
            return self.result_buffer.pop(int(earliest_i))
        else:
            raise ValueError("Cannot find oldest result in queue")

    def new_result(self, result):

        # store a few recent results to ensure triggers are handled on time
        self.result_buffer.append(result)  # add the result to the end of the buffer list
        if len(self.result_buffer) > self.result_buffer_size: 
            # find the OLDEST time in the buffer, and return it
            result = self.get_oldest_result()
        else:
            return # wait until the list is full to plot anything
        # print("new result, time = ", result["t"][0], result["t"][-1])
        # self.print_triggers()
        if len(self.triggers) == 0:  # no triggers - nothing to plot
            return

        time_arr = result["t"]
        TR = self.triggers[0] 
        if TR.trigger_time > time_arr[-1]:  # no trigger yet
            # print(f"*** Trigger detected at {TR.trigger_time:.4f} for time block: {time_arr[0]:.6f} - {time_arr[-1]:.6f}")
            # print(len(time_arr))
            return
        # Copy data from result to trigger buffer
        trigger_index = max(
            0, int(np.round((TR.trigger_time - time_arr[0]) / self.dt))
        )  # index of trigger within new data
        # number of points available is the smaller of the remainder of the current
        # part of the unused buffer, or the remainder of the time array
        npts = min(
            len(TR.buf) - TR.curr_buff_ptr, len(time_arr) - trigger_index
        )  # number of samples to copy from new data
        # print("Trigger index: ", trigger_index, "npts: ", npts, "TR.curr: ", TR.curr_buff_ptr)
        for k in self.plot_keys:  # self.plot_keys is a list, buf is a recarray
            if k in result.keys():
                TR.buf[k][TR.curr_buff_ptr : TR.curr_buff_ptr + npts] = result[k][trigger_index : trigger_index + npts]

        TR.curr_buff_ptr += npts
        # print('buff shape [0]: ', TR.buf.shape[0], time_arr[0])
        # print("n triggers: ", len(self.triggers))
        if TR.curr_buff_ptr >= TR.buf.shape[0]:
            # If the trigger buffer would run over on next call, plot the current buffer
            #  and remove it - TR.trigger_time
            self.plot_win.plot((np.arange(TR.buf.shape[0]) * self.dt), TR.buf, TR.info)
            self.triggers.pop(0)
            if len(time_arr) > npts:
                # If there is data left over, try feeding it to the next trigger
                result = result[trigger_index+npts:]
                self.new_result(result)

