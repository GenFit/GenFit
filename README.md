# GenFit Package

[![Build Status](https://travis-ci.com/GenFit/GenFit.svg?branch=master)](https://travis-ci.com/GenFit/GenFit)

GenFit is an experiment-independent framework for track reconstruction in particle and nuclear physics. It consists of three modular components:

* Track fitting algorithms

  Currently, GenFit contains a Kalman Filter, a Deterministic Annealing Filter, and a General Broken Lines fitter. Other algorithm modules can be added easily.

* Track representations

  These modules can perform extrapolations of track parameters through material and magnetic fields. GenFit is distributed with a well-tested track representation.
  Existing track extrapolation codes can be interfaced in a very straightforward way in this framework, using their native geometry and magnetic field interfaces.

* Measurements
  
  The measurement dimensionality and the orientation of planar tracking detectors can be chosen freely. GenFit is especially useful for tracking systems which include detectors which do not measure the passage of particles on predefined planes, like TPCs or wire-based drift chambers. The concept of so-called virtual detector planes provides a simple mechanism to use these detector hits in a transparent way without any geometrical simplifications.

GenFit has been developed in the framework of the PANDA experiment at FAIR, Darmstadt, Germany. It is also used in the Belle II, Fopi, and GEM-TPC experiments.

GenFit is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the license (LGPLv3), or (at your option) any later version. A copy of the license is distributed with GenFit in the `LICENSE.md` file.

GenFit Homepage: (https://github.com/GenFit/GenFit)
