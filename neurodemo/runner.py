# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""
from pyqtgraph.Qt import QtCore
from timeit import default_timer as def_timer
from datetime import timedelta


class SimRunner(QtCore.QObject):
    """Run a simulation continuously and emit signals whenever results are ready.
    
    Results are emitted with a dictionary containing all state variables for the
    final timepoint in the simulation, as well as the complete time-record for
    any variables that have been requested using add_request().
    """
    new_result = QtCore.Signal(object, object)  # final_state, requested_records
    
    def __init__(self, sim):
        QtCore.QObject.__init__(self)
        
        # dumps profiling data to prof.pstat
        # view with: python gprof2dot/gprof2dot.py -f pstats prof.pstat  | dot -Tpng -o prof.png && gwenview prof.png
        #from cProfile import Profile
        #self.prof = Profile()
        #import atexit
        #atexit.register(lambda: self.prof.dump_stats('prof.pstat'))
        #self.prof.enable()
        
        self.sim = sim
        self.speed = 1.0
        self.requests = []
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.run_once)
        self.counter = 0
        
    def add_request(self, key):
        self.requests.append(key)

    def remove_request(self, key):
        self.requests.remove(key)
        
    def start(self, blocksize=500, **kwds):
        self.starttime = def_timer()
        self.blocksize = blocksize
        self.run_args = kwds
        self.timer.start(20)  # determines the width of the display window/update interval
        
    def stop(self):
        self.timer.stop()
        
    def run_once(self):
        self.counter += 1
        now = def_timer()
        elapsed = now - self.starttime
        # print(f"run_once***  count={self.counter:d}  secs: {elapsed:f}")
        blocksize = int(max(2, self.blocksize * self.speed))
        # print(blocksize, self.blocksize, self.speed)

        result = self.sim.run(blocksize, **self.run_args)
        # print("    result retrieved")
        rec = {}
        for key in self.requests:
            try:
                data = result[key]
            except KeyError:
                #print("Key '%s' not found in sim result; skipping." % key)
                continue
            rec[key] = data
        # print("    rec filled")
        state = result.get_final_state()
        # print("    final state retrieved")
        # print(state, rec.keys())
        self.new_result.emit(state, rec)   # this pauses on x86_64 after some time
        # print("    new result emitted")
        
    def set_speed(self, speed):
        self.speed = speed
