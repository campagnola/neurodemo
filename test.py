# -*- coding: utf-8 -*-
"""
NeuroDemo - Physiological neuron sandbox for educational purposes
Luke Campagnola 2015
"""
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
from neurodemo import *

pg.setConfigOption('antialias', True)
app = QtGui.QApplication([])

# HH Simulation
#sim = Sim(temp=6.3)
#neuron = Section()
#leak = neuron.add(Leak(gbar=0.1*mS/cm**2))
#hhk = neuron.add(HHK())
#hhna = neuron.add(HHNa())

# Lewis & Gerstner cortical neuron
sim = Sim(temp=37)
neuron = Section(vm=-70*mV)
lgna = neuron.add(LGNa())
lgkf = neuron.add(LGKfast())
lgks = neuron.add(LGKslow())
leak = neuron.add(Leak(gbar=0.25*mS/cm**2, erev=-70*mV))

clamp = PatchClamp(mode='ic')
neuron.add(clamp)
sim.add(neuron)

win = pg.GraphicsWindow()
win.resize(1000, 600)
win.setWindowTitle('Testing hhSim.py')
p1 = win.addPlot(title='IC', labels={'left': ('Vm', 'V')})
p2 = win.addPlot(labels={'left': ('Ipip', 'A')}, row=1, col=0)
win.ci.layout.setRowFixedHeight(1, 150)
p3 = win.addPlot(row=2, col=0)
win.ci.layout.setRowFixedHeight(2, 150)

dur = 100 * ms
npts = int(dur / sim.dt)
x1 = int(20*ms / sim.dt)
x2 = int(80*ms / sim.dt)
x = np.linspace(-200, 200, 11) * pA
#x = [0*pA]
cmd = np.zeros((len(x), npts)) #*-65e-3
data = np.zeros((len(x), npts, 9))
t = np.arange(npts) * sim.dt
for i, v in enumerate(x):
    print('V: ', v)
    cmd[i, x1:x2] = v
    clamp.queue_command(cmd[i], sim.dt)
    #data[i] = run(neuron, mode='ic', dt=dt, cmd=cmd[i])
    result = sim.run(npts)
    data = result[neuron, 'V']
    p1.plot(t, data, pen=(i, 15))
    p2.plot(t, cmd[i], pen=(i, 15))
    #p2.plot(t, clamp.current(result), pen=(i, 15))
    
    #p3.plot(t, leak.current(result), pen=(i, 15))
    #p3.setLabel('left', 'Ileak', 'A')

    #p3.plot(t, lgkf.current(result), pen=(i, 15))
    #p3.setLabel('left', 'LG Kfast Current', 'A')

    #p3.plot(t, lgks.open_probability(result), pen=(i, 15))
    #p3.setLabel('left', 'LG Kslow O.P.')

    #p3.plot(t, lgna.current(result), pen=(i, 15))
    #p3.setLabel('left', 'LG Na Current', 'A')

    p3.plot(t, lgna.open_probability(result), pen=(i, 15))
    p3.setLabel('left', 'LG Na O.P.')

    #p3.plot(t, hhna.open_probability(result), pen=(i, 15))
    #p3.setLabel('left', 'HH Na O.P.')

import sys
if sys.flags.interactive == 0:
    QtGui.QApplication.instance().exec_()
