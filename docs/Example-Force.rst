
.. _example-force:

##################################################
`ottForce` Example : Calculate and visualise force
##################################################

This example shows how to calculate and visualise force on a sphere.
This example includes two visualisations: a 1-D force profile in the x/z
directions, and a 2-D force profile in the x-z plane.

The example assumes the user is already familiar with the beam/particle
classes.  For examples of the beam/particle functionality, see the
:ref:`example-beams` and :ref:`example-particles` sections of the
documentation.

The code for this example can be found in the examples directory,
if OTT has been added to the path, run::

   open examples/ottBeams.m

A live script is also available for this example::

   open examples/liveScripts/beams.mlx

This example assumes the toolbox has been added to the Matlab path,
for information, see :ref:`adding-ott-to-matlabs-path`.

.. contents:: Concents
   :depth: 3
   :local:
..

Setup the particle
==================

The first step is to setup the particle.
For this simulation we use a spherical particle with a 1 micron radius.
We describe the geometry and let the ``FromShape`` method choose an
appropriate T-matrix method (for spheres, this will be ...
:class:`+ott.+tmatrix.Mie`)::

   % Create geometry for shape
   radius = 1.0e-6;      % Sphere radius [m]
   shape = ott.shape.Sphere(radius);

   % Setup particle
   particle = ott.particle.Particle.FromShape(shape, 'index_relative', 1.2);

Setup the beam
==============

The next step is to calculate the incident beam.
The following creates a Gaussian beam tightly focussed by a microscope
objective.  The back aperture of the microscope objective creates a
hard edge, we model this by setting the ``truncation_angle`` parameter
of the Gaussian beam class::

   % Calculate truncation angle for modelling a microscope objective
   NA = 1.2;
   index_medium = 1.33;
   truncation_angle = asin(NA/index_medium);

   % All beam parameters have SI units (radians for angles)
   beam = ott.beam.Gaussian(...
       'waist', 1.0e-6, 'power', 0.01, 'truncation_angle', truncation_angle, ...
       'index_medium', index_medium, 'wavelength0', 1064e-9);

Generate beam visualisation
===========================

It is always good to check the beam looks how we would expect.
If the beam doesn't look correct then there may be something wrong with
the choice of parameters or some parameters may be outside the range
of values for which :class:`+ott.+beam.Gaussian` gives reasonable results
(in which case, consider using the underlying :mod:`+ott.+bsc` classes
directly).  We can visualise the above beam using::

   figure();
   subplot(1, 2, 1);
   beam.visFarfield('dir', 'neg');
   title('Far-field');
   subplot(1, 2, 2);
   beam.visNearfield();
   title('Nearfield');

Calculate forces
================

For spherical particles in beams with only a single focus, we are often
interested in the 1-dimensional position-force profiles.  These can
be calculated by simply using the beam's ``force`` method and specifying
a ``3xN`` array of locations.  For the axial force::

   % Calculate force along axis
   z = linspace(-1, 1, 100)*6*beam.wavelength;
   fz = beam.force(particle, 'position', [0;0;1].*z);

For the radial force, we are interested in the force near the equilibrium
(both in the axial and radial directions)
The following code uses the axial force profile and then
uses the :class:`+ott.+tools.FindTraps1d` class to estimate where the
axial equilibrium is.
The radial force profile is then calculated at points passing through
the axial equilibrium. If the particle/beam does not have an axial
equilibrium then the coordinate origin is used::

   % Find traps along axis
   traps = ott.tools.FindTraps1d.FromArray(z, fz(3, :));

   % Get the axial equilibrium (might have no trap, in which case use z=0)
   if isempty(traps)
     z0 = 0.0;
   else
     z0 = traps(1).position;
   end

   % Calculate radial force
   r = linspace(-1, 1, 100)*6*beam.wavelength;
   fr = beam.force(particle, 'position', [1;0;0].*r + [0;0;z0]);

The particle/beam functions use SI units and the resulting force/torque
have units of [Newtons] and [Newton meters].
The following generates plots of the force profiles::

   subplot(1, 2, 1);
   plot(z, fz);
   xlabel('z position [m]');
   ylabel('force [N]');
   legend({'X', 'Y', 'Z'});

   subplot(1, 2, 2);
   plot(r, fr);
   xlabel('x position [m]');
   ylabel('force [N]');
   legend({'X', 'Y', 'Z'});


If the beam has an orbital component (for example, if there is orbital
or spin angular momentum) then there may be a force component around
the beam axis.

Generate 2-D visualisation of force
===================================

When the trap is harmonic, the axial/radial force profiles are sufficient
to characterise the properties of the trap.  However, this is rarely
the case for complex beams or particles far from equilibrium.

An alternative method for visualising the optical force is to plot a
colour map or vector field of the optical force as a function of position.
The following plots such a visualisation (this may take some time
depending on how many positions are required for the visualisation)::

   r = linspace(-1, 1, 40)*6*beam.wavelength;
   z = linspace(-1, 1, 40)*6*beam.wavelength;

   [R, Z] = meshgrid(r, z);
   F = beam.force(particle, 'position', {R, 0*R, Z});

   % Calculate magnitude of force
   Fmag = vecnorm(F, 2, 3);

   % Generate colormap/quiver plot showing force
   imagesc(R, Z, Fmag);
   hold on;
   quiver(R, Z, Fmag(:, :, 1), Fmag(:, :, 3));
   hold off;
   xlabel('x position [m]');
   ylabel('y position [m]');
   cb = colorbar();
   yaxis(cb, 'Force [N]');

Further reading
===============

The main part of an optical tweezers simulation is force calculation,
however another important component is being able to simulate how the
particle moves, this requires calculating the drag and potentially thermal
motion, and simulating the particle for a while.
The toolbox has preliminary support for dynamics simulations using a fixed
time step (see the ``ottDynamics.m`` example and the
:class:`+ott.+tools.Dynamics` class).

