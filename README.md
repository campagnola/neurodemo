Neuron Demonstration
====================

Luke Campagnola, 2015


This is an educational simulation of a simple neuron.

* Hodgkin & Huxley channels
* Lewis & Gerstner (2002) cortical channels
* Destexhe 1993 Ih channel
* Current/voltage clamp electrode with access resistance


Requirements
------------

* Python 3
* NumPy, SciPy
* PyQt4
* PyQtGraph


Installation
------------

1. Experienced Python users can install the requirements listed above, then download the neurodemo code and run `demo.py`. For everybody else, I recommend installing the [Anaconda Python distribution] (http://continuum.io/downloads). Click on "I want Python 3.4" and download the installation package for your platform. To keep the size down, you can instead install [miniconda] (http://conda.pydata.org/miniconda.html). When you run the installer, _take note of the location it is installing to_.

2. If you installed miniconda, then it is necessary to manually install numpy, scipy, and pyqt4 (if you installed Anaconda, then these packages are already installed and you can skip this step). This can be done from a terminal:

```
> conda install numpy scipy pyqt
```

3. Install pyqtgraph (this is required regardless of which python distribution you installed):

```
> pip install pyqtgraph
```

4. Download the [neurodemo source code] (https://github.com/campagnola/neurodemo/archive/master.zip). You can find this by going to http://github.com/campagnola/neurodemo and clicking the "Download ZIP" button on the right side of the page. Unzip the file and you are ready to begin!

Running the Demo
----------------

You can launch the demo from the terminal by navigating to the location where you extracted the source code (eg, `cd Downloads\neurodemo\`) and then typing `python demo.py`.

Windows: Open the folder where you extracted the neurodemo source code. Right click `demo.py` and select "Open with..". Navigate to the location where you installed anaconda or miniconda (you may need to click something like "browse" or "choose another application" depending on your version of windows). Select `python.exe`.

OSX:


Linux: 














