
.. _example-beams:

###################################################
`ottBeams` Example : Creating and visualising beams
###################################################

This example demonstrates how to create and use different kinds of
optical tweezers toolbox (OTT) beams.  The code for this example can
be found in the examples directory, if OTT has been added to the path,
run::

   open examples/ottBeams.m

A live script is also available for this example::

   open examples/liveScripts/beams.mlx

This example assumes the toolbox has been added to the Matlab path,
for information, see :ref:`adding-ott-to-matlabs-path`.

.. contents:: Concents
   :depth: 3
   :local:
..

Creating a beam
===============

OTT provides two main methods for creating beams: using the classes in
the :mod:`+ott.+beam` package and using the classes in the
:mod:`+ott.+bsc` package.
The :mod:`+ott.+beam` package is intended to provide a easy to use, high-level
interface for interacting with beams.  Internally, the :mod:`+ott.+beam`
classes use the functions declared in :mod:`+ott.+bsc`, and in some cases,
using the :mod:`+ott.+bsc` classes directly can give better run-times.
However, for most use cases, the :mod:`+ott.+beam` classes should be adequate.
In this example, we will focus on the :mod:`+ott.+beam` package.

The :mod:`+ott.+beam` package provides several classes
representing common beams, for example, to create a Gaussian beam, call::

   beam = ott.beam.Gaussian();

This create the default Gaussian beam (see the
:class:`+ott.+beam.+Gaussian` documentation for the
exact specification).  We can change the beam properties such as the
medium speed and power using the beam's properties or methods,
for example::

  beam.power = 0.1;           % Set power [Watts]
  beam.speed = 3e8/1.33;      % Set speed [meters/second]

Some properties, such as wavelength, need to be set using beam methods.
Beams are not handle classes, as such, it is important to store the
resulting beam object.  The following sets the wavelength by keeping the
speed fixed::

  beam = beam.setWavelength(532e-9, 'fixedSpeed'); % Set wavelength [meters]

Alternatively, we can create the beam with the desired properties at
the start by passing a list of named arguments to the constructor::

  beam = ott.beam.Gaussian('power', 0.1, ...
      'index_medium', 1.33, 'omega', 2*pi*3e8/532e-9);

Many of the :mod:`+ott.+beam` classes provide alternative construction
methods for convenience, for example, the Gaussian and Laguerre--Gaussian
beams provide a ``FromNa`` method that accepts a numerical aperture value
as the first argument.

To view our beam, we can use one of the beam visualisation methods,
for example, the following creates a XY slice through the beam focus::

   beam.visNearfield();

This creates a new plot in the current figure window (or creates a new figure
if required).  Alternatively, we can request the resulting image data and
plot the results ourselves (this is usually done in conjunction with
specifying the image range)::

   xrange = linspace(-1, 1, 80)*1e-6;   % Range in meters
   yrange = xrange;
   im = beam.visNearfield('range', {xrange, yrange});
   figure();
   contour(xrange, yrange, im);
   xlabel('X Position [m]');
   ylabel('Y Position [m]');

Translations and rotations
==========================

Beams have a ``position`` and ``rotation`` property.  These properties
are applied to the beam whenever the beam is used (for example, when a
visualisation method is called or when ``getData`` is called).

The position property is a 3x1 numeric vector.  To shift the beam by
1 wavelength in the x direction, we can directly set the position property::

   beam.position = [1;0;0]*beam.wavelength;

Alternatively, we can use the
:meth:`+ott.+utils.TranslateHelper.translateXyz` method.  The translation
method applies the translation on top of any existing displacement and
returns a new copy of the beam, for example, to translate our previous beam
along the Y-direction, we could use::

   tbeam = beam.translateXyz([0;1;0]*beam.wavelength);

Rotations are stored as 3x3 rotation matrices.  As with the ``position``
property, we can also directly set the ``rotation`` property, however
it is often easier to use the ``rotate*`` methods from
:class:`+ott.+utils.RotateHelper`.
The following rotates the beam pi/2 radians about the Y axis.
When the beam is used, the rotation is applied before the translation::

   rbeam = beam.rotateY(pi/2);

Arrays of beams
===============

The toolbox supports three kinds of arrays: Coherent arrays, Incoherent
arrays, and generic arrays.  Incoherent/Coherent arrays represent beams
which can be represented by a finite set of sub-beams.  Generic arrays
are simply collections of multiple beams.

To create a generic array, simply use Matlab's array syntax, for example::

   beams = [ott.beam.Gaussian(), ...
       ott.beam.LaguerreGaussian('lmode', 10)];

Most operations can be applied to generic arrays.  The result is the
same as applying the operation to each element of the array.  For example,
to translate the array of beams::

   tbeams = beams.translateXyz([1;0;0]*beam.wavelength);

Or to set the position of each beam::

   [tbeams.position] = deal([1;0;0]*beam.wavelength);

Coherent and Incoherent arrays can be created using the
:class:`+ott.+beam.Coherent` and :class:`+ott.+beam.Incoherent`,
for example::

   cbeams = ott.beam.Coherent(beams);


Calculating the change in momentum between two beams
====================================================

The :class:`+ott.+beam.Beam` class provides methods for calculating the change
in momentum between two beams.  Although it is more common to calculate
the force acting on a particle (see the ``ottForce.m`` example), the following
shows how to calculate the change in momentum between two beams::

   beam1 = ott.beam.Gaussian();
   beam2 = beam1.rotateY(pi/2);
   force = beam1.force(beam2)

Creating a custom beam
======================

Although the toolbox has several different beams commonly used in
optical trapping (for a complete list, see the `beam` package
reference section), it is sometimes necessary to create a custom beam.
The most common scenario is when modelling an SLM or the experimentally
measured field at the back aperture of the focussing objective.  For this
task we can use the `PmParaxial` class (for more control over the fields
we could also use the `ott.bsc` classes).  The following example shows
how we could model a phase-only SLM illuminated by a Gaussian-like beam::

   % Generate coordinates for pattern
   x = linspace(-1, 1, 20);
   y = linspace(-1, 1, 20);
   [X, Y] = ndgrid(x, y);

   % Calculate incident field
   E0 = exp(-(X.^2 + Y.^2)./4);

   % Calculate SLM-like pattern
   kx = 2;
   phi = 2*pi*kx*x;

   % Calculate field at back aperture
   E = E0 .* exp(1i*phi);

   % Calculate beam
   beam = ott.beam.PmParaxial.InterpProfile(X, Y, E);

Further reading
===============

For the full range of beams currently inculded in the toolbox, refer to
the :ref:`beams-package` part of the reference section.
The example code used to generate the overview figure in the reference
section can be found in the ``examples/packageOverview/`` directory.
More advanced beam functionality can be implemented by directly using
the beam shape coefficient classes (the :mod:`+ott.+bsc` package).
For a example which uses both :class:`+ott.+beam.Beam` and
:class:`+ott.+bsc.Bsc`, see ``examples/ottLandscape.m``.

