import sys
import pyqtgraph as pg
from . import qt
import pyqtgraph.parametertree as pt
#import qdarktheme
STYLE_SHEET = """
QPushButton {
    color: #8ab4f7;
    border: 1px solid #3f4042;
    padding: 4px 8px;
    border-radius: $radius{4px};
}
QPushButton:!window {
    background: transparent;
}
"""
junk = """
QPushButton:flat,
QPushButton:default {
    border: none;
    padding: 5px 9px;
}
QPushButton:default {
    color: #202124;
    background: #8ab4f7;
}
QPushButton:hover,
QPushButton:flat:hover {
    background: rgba(46.000, 70.000, 94.000, 0.333);
}
QPushButton:pressed,
QPushButton:flat:pressed,
QPushButton:checked:pressed,
QPushButton:flat:checked:pressed {
    background: rgba(46.000, 70.000, 94.000, 0.933);
}
QPushButton:checked,
QPushButton:flat:checked {
    background: rgba(46.000, 70.000, 94.000, 0.733);
}
QPushButton:default:hover {
    background: #7fa7e5;
}
QPushButton:default:pressed {
    background: #6d8bbe;
}
QPushButton:default:disabled {
    background: #53575b;
}
"""
import numpy as np

app = pg.mkQApp()
# app.setStyle("Fusion") 
#app.setStyle("Fusion")  # necessary to remove double labels on mac os
# note that the default:
# print("Current style: ", app.style().name())
# returns "macos".
#app.setStyleSheet(qdarktheme.load_stylesheet())
#app.setStyleSheet("QPushButton {color: #8ab4f7;border: 1px solid #3f4042;padding: 4px 8px;border-radius: $radius{4px};}")
#app.setStyleSheet("QPushButton:!window {background: opaque; border: 2px solid #8f8f91;} QPushButton:!open {color:#aaaaaa;} QPushButton:!closed {color:#aa0000;}")
# app.setStyleSheet("QPushButton:pressed {background-color: qt.qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,stop: 0 #dadbde, stop: 1 #f6f7fa);}")
# app.setStyleSheet("QPushButton:!window {background: opaque;}") 
class PGTest(qt.QWidget):
    def __init__(self):

        qt.QWidget.__init__(self)
        self.fullscreen_widget = None
        self.resize(320, 320)
        self.layout = qt.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        self.ptree = pg.parametertree.ParameterTree(showHeader=False)
        self.layout.addWidget(self.ptree)
        self.params = pt.Parameter.create(name='params', type='group', children=[
            dict(name='Run/Stop', type='action', value=False),
        ])
        # self.params.setOpts(name = "Forward/Backward")
       

        self.ptree.setParameters(self.params)
        self.params.sigTreeStateChanged.connect(self.params_changed)
        self.params.setOpts(name = "Forward/Backward !")  # sets group paramater name (default was "params")
        self.params.setOpts(title = "Up/Down >") # adds options dict with Title  (also in parameter name)
        self.show()
    
    def params_changed(self, root, changes):
        print(root)
        print(changes)
        pass

if __name__ == "__main__":
    win = PGTest()
    if sys.flags.interactive == 0:
        pg.exec()