
.. _example-particles:

#####################################################
`ottParticles` Example : Creating and using particles
#####################################################

This examples demonstrates how to create and use different kinds of
optical tweezers toolbox (OTT) particles.  The code for this example can
be found in the examples directory, if OTT has been added to the path,
run::

   open examples/ottParticles.m

A live script is also available for this example::

   open examples/liveScripts/particles.mlx

This example assumes the toolbox has been added to the Matlab path,
for information, see :ref:`adding-ott-to-matlabs-path`.

.. contents:: Concents
   :depth: 3
   :local:
..

Creating a particle with a simple geometry
==========================================

OTT particles encapsulate the geometry, optical scattering and fluid drag
information for a particle in an optical tweezers simulation.
There are several ways to create particles, the simplest involve using the
`FromShape` methods of the particle classes.  These methods use the
corresponding `FromShape` methods of the drag and T-matrix classes, however
care should be taken when using large or complex particles, as these methods
may not always give good approximations for the input geometry.
The following create a shape and encapsulates it in a particle::

   % Generate a geometric shape (units of meters)
   shape = ott.shape.Sphere(1e-6);

   % Create a particle
   particle = ott.particle.Particle.FromShape(shape);

We can visualise the particle using the ``surf`` method, this simply calls
the shape's corresponding ``surf`` method::

   particle.surf();

Creating a particle with custom properties
==========================================

Particle properties don't have to be related: this is useful when there
is no available method for modelling the drag or scattering of a particle,
or when we want to visualise the shape using a simpler geometry.

To create a particle with custom properties, we can use the `Fixed`
sub-class and provide our own T-matrix, drag and geometry.
This class simply stores the provided properties.

For instance, we can create a particle instance with the geometry
of a cylinder, T-matrix for a spheroid, and drag for a cylinder::

   % Create geometry
   wavelength = 1e-6;
   shape = ott.shape.Cylinder(wavelength, wavelength);

   % Create drag using FromShape
   drag = ott.drag.Stokes.FromShape(shape, 'viscosity', 0.001);

   % Create T-matrix for spheroid.
   % There are several T-matrix methods included in the toolbox, Smarties
   % works well for spheroidal particle.  This example creates a spheroid
   % with similar dimensions to our Cylinder.  T-matrix methods use distance
   % units with the dimensions of wavelength (particle/beams use meters).
   index_relative = 1.2;
   tmatrix = ott.tmatrix.Smarties(...
       shape.radius ./ wavelength, shape.height ./ wavelength, index_relative);

And then we can store these properties in a particle instance::

   particle = ott.particle.Fixed(shape, 'drag', drag, 'tmatrix', tmatrix);


Translations and rotations
==========================

Similar to beams, particles have ``rotation`` and ``position`` properties
which can be used to control the position/rotation of the particle.

Both rotation and position can be directly set, for example::

   particle.position = [1;0;0]*wavelength;

And properties can be adjusted using the translate/rotate methods.
As with beams, the rotate/translate methods return a copy of the object::

   new_particle = particle.rotateY(pi/2);

We can see the effect of these operations by generating a surface plot
with the ``surf`` method::

   new_particle.surf();

For additional example, see the rotation/position part of the
:ref:`example-beams` (``ottBeam.m``) example.

Calculate optical scattering
============================

To calculate how a beam is scattered by a particle, we can either use
the matrix multiplication operator or the beam scatter method.

For this example, lets use a Gaussian beam::

   beam = ott.beam.Gaussian();

And instead of only calculating the external fields, we can also calculate
the internal fields by passing ``'internal', true`` to the particle
constructor::

   shape = ott.shape.Sphere(1e-6);
   particle = ott.particle.Particle.FromShape(shape, ...
       'internal', true, 'index_relative', 1.2);

Not all T-matrix calculation methods support calculating internal fields.
The T-matrix method that is used depends on the geometry (for this case,
a sphere, the internal method should give fairly accurate results).

To calculate the force, we multiply the particle by the beam (note that
the order of operations is important)::

   sbeam = particle * beam;

Alternatively, we can use the ``scatter`` method (this also supports
applying position/rotations to the beam)::

   sbeam = beam.scatter(particle);

The scattered beam stores an instance of the particle and the incident
beam, allowing us to easily visualise the internal and external fields,
for example, the following outputs the fields shown in [TODO]::

   sbeam.visNearfield('axis', 'y');


Calculate forces
================

Force can be calculated either directly using the ``force`` method of the
scattered beam or using the ``force`` method of the incident beam.
With the scattered beam::

   force = sbeam.force();

And, with the incident beam::

   force = beam.force(particle);

With both methods, the resulting force has units of Newtons.
The incident beam method has the advantage that we can also specify
a 3xN array of positions or rotations to apply to the particle, these
are applied on top of any existing particle translations/rotations::

   positions = randn(3, 5)*1e-6;
   forces = beam.force(particle, 'position', positions);

Additional force calculation examples are provided in the
:ref:`example-force` (``ottForce.m``) example.

Further reading
===============

The :ref:`example-advanced` section shows how the particle class can be
used for various tasks including dynamics simulations.

T-matrices cannot be calculated for all shapes, but if the particle is
homogeneous and small enough to simulate, it should be possible to compare
the scattering by the T-matrix to the scattering directly with DDA, giving
an estimate for accuracy.

DDA should work with most shapes, but this hasn't been thoroughly tested,
if you find something interesting, let us know.
For inspiration with creating different particle shapes, take a look at
the :mod:`+ott.+shape` reference section and the shapes used in the
``examples/packageOverview`` examples.

