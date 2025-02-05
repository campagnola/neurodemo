# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""
from . import qt
from timeit import default_timer as def_timer


class SimRunner(qt.QObject):
    """Run a simulation continuously and emit signals whenever results are ready.    
    """
    new_result = qt.Signal(object)
    
    def __init__(self, sim):
        qt.QObject.__init__(self)
        
        # dumps profiling data to prof.pstat
        # view with: python gprof2dot/gprof2dot.py -f pstats prof.pstat  | dot -Tpng -o prof.png && gwenview prof.png
        #from cProfile import Profile
        #self.prof = Profile()
        #import atexit
        #atexit.register(lambda: self.prof.dump_stats('prof.pstat'))
        #self.prof.enable()
        
        self.sim = sim
        self.speed = 1.0
        self.timer = qt.QTimer()
        self.timer.timeout.connect(self.run_once)
        self.counter = 0

    def start(self, blocksize=500, **kwds):
        self.starttime = def_timer()
        self.blocksize = blocksize

        self.run_args = {"stop_after_cmd": False}  # If True, then stop after all queued commands are exhausted
        self.run_args.update(**kwds)
        self.timer.start(20)  # Argument (in milliseconds) determines width of the display window/update interval
        
    def stop(self):
        self.timer.stop()
        
    def running(self):
        return self.timer.isActive()

    def run_once(self):
        self.counter += 1
        blocksize = int(max(2, self.blocksize * self.speed))
        result = self.sim.run(blocksize, **self.run_args)
        self.new_result.emit(result)

    def set_speed(self, speed):
        self.speed = speed
