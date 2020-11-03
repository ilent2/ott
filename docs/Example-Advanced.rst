
.. _example-advanced:

#################
Advanced examples
#################

This section provides a very brief overview of the more advanced examples
currently included in the toolbox.  These examples are all contained
in the ``examples/`` directory.

Previous versions of the toolbox included other examples, most functionality
is still present in this version of the toolbox but examples have not
all been included.

.. contents:: Contents
   :depth: 3
   :local:
..

Landscape
=========

Calculates trapping landscapes.  These show how optical trap properties,
such as trap stiffness or maximum optical force, vary as a function of
particle properties (such as radius and refractive index).
The ``ottLandscape.m`` example calculates how maximum trap depth varies
with different size and refractive index spherical particles, example
output is shown in [TODO].

These types of calculations often involve looping over parameters.
In some cases it is more optimal to pre-compute translations or rotations,
such as (in the example) translations along the beam axis for different
sized particles.
This version of the toolbox also adds support for re-using VSWF data
(for faster T-matrix calculation and field calculation), this is a feature
which may be useful for certain types of trapping landscapes.
See :class:`+ott.+utils.VswfData` and its usage in the various beam and
T-matrix classes.

Dynamics
========

The ``ottDynamics.m`` example shows how the :class:`+ott.+tools.Dynamics`
class can be used to simulate particle dynamics.
This requires a beam, T-matrix and drag tensor (which can be provided
as a :class:`+ott.+particle.Particle` instance); and solves the Langevin
equation using a fixed time step method.

The example creates a spherical particle and generates a trace of the
position and a 2-D histogram of the XY position, example output is shown
in [TODO].

Wall Effects
============

This ``ottWallEffects.m`` example builds on the ``ottDynamics.m`` example
by using a Sphere-Wall drag tensor.
Currently there is not optical interaction between the particle/beam/wall,
this may be added in a future release.
The example generates visualisations of [TODO], example output is shown
in [TODO].

Non-spherical particles
=======================

This version of the toolbox is focussed on T-matrix methods for scattering
calculations.  Future versions of the toolbox aim to include DDA,
shape surface approximation, and geometric optics.  For now, the main
methods for modelling almost-arbitrary shaped particles are point-matching,
extended boundary conditions method (EBCM), and discrete dipole
approximation.

The ``ottNonSpherical.m`` example shows how EBCM, DDA and PM can be used
to model a cylindrical particle.  For certain sizes, all three methods
agree well, other sizes/aspect ratios, the methods start to disagree,
as shown in [TODO].
Calculating T-matrices using different methods is a good method for
validating the generated T-matrix/predicted force.

DDA Vaterite
============

The ``ottDdaVaterite.m`` example shows some of the more advanced features
of the DDA implementation: the implementation can be used to model a
inhomogeneous birefringent particle.
The example calculates T-matrices for vaterites with a sheaf-of-wheat
structure and calculates the torque about the z axis.  Example output
is shown in [TODO].

The DDA implementation also supports calculating only specific columns of
the T-matrix.  This can be useful when the T-matrix is illuminated by only
a beam with particular orders; or when calculating T-matrices in parallel.

Neural Networks for fast simulation
===================================

The ``ottNeuralNetwork.m`` example shows how a neural network can be
trained to rapidly predict optical force data.  This requires first
generating a large training data set which is used to train the network.
Once trained, the network can be used for rapid force calculations
(for inputs within the bounds of the training data) or easily shared
with other researchers.

For this particular problem, interpolation could also be used, however
the stored data sets are often much larger than the resulting neural
network.
Additionally, for problems with additional inputs/outputs (such as
predicting force and torque from orientation and position), direct
interpolation becomes a much more computationally expensive task.

Further reading
===============

These examples cover most of the advanced functionality in the toolbox.
For information about how different methods work, the reference section
provides some insight.
Otherwise, take a look at the :ref:`Further-Reading` section or the
previous version of the toolbox for additional examples.
Experimental (viz., incomplete) new features can also be found in version 2
of the toolbox (available on GitHub).

