# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""
from pyqtgraph.Qt import QtGui, QtCore


class SimRunner(QtCore.QObject):
    """Run a simulation continuously and emit signals whenever results are ready.
    """
    new_result = QtCore.Signal(object, object)  # last_record, current_state
    
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
        self.requests = {}
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.run_once)
        
    def add_request(self, key, req):
        self.requests[key] = req

    def remove_request(self, key):
        self.requests.pop(key)
        
    def start(self, blocksize=500, **kwds):
        self.blocksize = blocksize
        self.run_args = kwds
        self.timer.start(16)
        
    def stop(self):
        self.timer.stop()
        
    def run_once(self):
        blocksize = int(max(2, self.blocksize * self.speed))
        result = self.sim.run(blocksize, **self.run_args)
        
        rec = {}
        for key, req in self.requests.items():
            try:
                if isinstance(req, tuple):
                    if req[1] == 'I':
                        data = req[0].current(result)
                    elif req[1] == 'G':
                        data = req[0].conductance(result)
                    elif req[1] == 'OP':
                        data = req[0].open_probability(result)
                    else:
                        data = result[req]
                else:
                    data = result[req]
            except KeyError:
                print("Key '%s' not found in sim result; skipping." % key)
                continue
            rec[key] = data
        
        state = self.sim.get_current_state()
        
        self.new_result.emit(rec, state)
        
    def set_speed(self, speed):
        self.speed = speed
