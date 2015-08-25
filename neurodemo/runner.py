from pyqtgraph.Qt import QtGui, QtCore


class SimRunner(QtCore.QObject):
    """Run a simulation continuously and emit signals whenever results are ready.
    """
    new_result = QtCore.Signal(object)
    
    def __init__(self, sim):
        QtCore.QObject.__init__(self)
        self.sim = sim
        self.speed = 1.0
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.run_once)
        
    def start(self, blocksize=500, **kwds):
        self.blocksize = blocksize
        self.run_args = kwds
        self.timer.start(16)
        
    def stop(self):
        self.timer.stop()
        
    def run_once(self):
        blocksize = int(max(2, self.blocksize * self.speed))
        res = self.sim.run(blocksize, **self.run_args)
        self.new_result.emit(res)
        
    def set_speed(self, speed):
        self.speed = speed
