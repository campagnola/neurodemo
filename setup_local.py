from setuptools import setup, find_packages, Extension
# from Cython.Build import cythonize
import numpy

version = '1.1'

setup(name='neurodemo',
      version=version,
      description='Simulations as a teaching tool for neuron electrophysiology',
      url='http://github.com/lcampagnola/neurodemo',
      author='Luke Campagnola',
      author_email='',
      license='MIT',
      packages=find_packages(include=['neurodemo*']),

      python_requires='>=3.10',

      zip_safe=False,
      entry_points={
          'console_scripts': [
               'ndemo=demo:main',
               ],
      },
      classifiers = [
             "Programming Language :: Python :: 3.10+",
             "Development Status ::  Release",
             "Environment :: Console",
             "Intended Audience :: Graduate Students",
             "License :: MIT",
             "Operating System :: OS Independent",
             "Topic :: Software Development :: Tools :: Python Modules",
             "Topic :: Computational Modeling :: Neuroscience",
             ],
    )
