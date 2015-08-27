from neurodemo import MultiClamp
import numpy as np


def test_multiclamp_cmd():
    global mc
    mc = MultiClamp()
    mc.set_holding('ic', -1.0)
    a = np.ones(100)
    a[::2] = 2
    dt = 1e-3
    assert mc.queue_command(a, dt) == dt
    try:
        mc.queue_command(a, dt, start=2e-3) 
        raise AssertionError('Should have raised ValueError.')
    except ValueError:
        pass
    
    assert np.allclose(mc.queue_command(a, 1e-3, start=.2), 0.2)
    assert np.allclose(mc.queue_command(a, 1e-3), 0.3)
    
    vals = [
        (0.0, -1.0),
        (0.5, 0.5),
        (1.0, 2.0),
        (1.5, 1.5),
        (2.0, 1.0),
        (100.0, 1.0),
        (100.5, 0.0),
        (101.0, -1.0),
        (150.0, -1.0),
        (199.0, -1.0),
        (199.5, 0.5),
        (200.0, 2.0),
        (299.0, 1.0),
        (299.5, 1.5),
        (300.0, 2.0),
        (399.5, 0.0),
        (400.0, -1.0),
        (500.0, -1.0),
    ]
    for x,y in vals:
        t = x * dt
        y1 = mc.get_cmd(t)
        if not np.allclose(y1, y):
            raise ValueError("Expected %f for time %f; got %f." % (y, t, y1))
    


if __name__ == '__main__':
    test_multiclamp_cmd()
