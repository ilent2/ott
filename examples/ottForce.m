% ottForce.m -- Calculate and visualise force on a sphere
%
% This example shows how to calculate and visualise force on a sphere.
% This example includes two visualisations: a 1-D force profile in the x/z
% directions, and a 2-D force profile in the x-z plane.
%
% The example assumes the user is already familiar with the beam/particle
% classes.  For examples of the beam/particle functionality, see their
% corresponding examples.
%
% For more information, see the corresponding example in the OTT documentation
% and the related functions in the reference sections.
% A interactive live script of this example is also available in the
% ``examples/liveScripts/`` directory.
%
% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is an example from the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../');

%% Setup particle
% The first step is to setup the particle.
% For this simulation we use a spherical particle with a 1 micron radius.
% We describe the geometry and let the `FromShape` method choose an
% appropriate T-matrix method (for spheres, this will be `ott.tmatrix.Mie`).

% Create geometry for shape
radius = 1.0e-6;      % Sphere radius [m]
shape = ott.shape.Sphere(radius);

% Setup particle
particle = ott.particle.Particle.FromShape(shape, 'index_relative', 1.2);

%% Setup beam
% The next step is to calculate the incident beam.
% The following creates a Gaussian beam tightly focussed by a microscope
% objective.  The back aperture of the microscope objective creates a
% hard edge, we model this by setting the `truncation_angle` parameter
% of the Gaussian beam class.

% Calculate truncation angle for modelling a microscope objective
NA = 1.2;
index_medium = 1.33;
truncation_angle = asin(NA/index_medium);

% All beam parameters have SI units (radians for angles)
beam = ott.beam.Gaussian(...
    'waist', 1.0e-6, 'power', 0.01, 'truncation_angle', truncation_angle, ...
    'index_medium', index_medium, 'wavelength0', 1064e-9);

%% Generate visualisation of beam
% It is always good to check the beam looks how we would expect.
% If the beam doesn't look correct then there may be something wrong with
% the choice of parameters or some parameters may be outside the range
% of values for which `ott.beam.Gaussian` gives reasonable results (in
% which case, consider using the `ott.bsc` classes directly).

figure();
subplot(1, 2, 1);
beam.visFarfield('dir', 'neg');
title('Far-field');
subplot(1, 2, 2);
beam.visNearfield();
title('Nearfield');

%% Calculate axial and radial force profiles
% For spherical particles in beams with only a single focus, we are often
% interested in the 1-dimensional position-force profiles.  These can
% be calculated by simply using the beam's `force` method and specifying
% a 3xN array of locations.
%
% The following code first calculates the axial force profile and then
% uses the `FindTraps1d` class to estimate where the axial equilibrium is.
% The radial force profile is then calculated at points passing through
% the axial equilibrium. If the particle/beam does not have an axial
% equilibrium then the coordinate origin is used.

% Calculate force along axis
z = linspace(-1, 1, 100)*6*beam.wavelength;
fz = beam.force(particle, 'position', [0;0;1].*z);

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

%% Generate plot of force along x and z
% The particle/beam functions use SI units and the resulting force/torque
% have units of [Newtons] and [Newton meters].
% The following generates plots of the force profiles.
% If the beam has an orbital component (for example, if there is orbital
% or spin angular momentum) then there may be a force component around
% the beam axis.

figure();

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

%% Generate 2-D visualisation of force
% When the trap is harmonic, the axial/radial force profiles are sufficient
% to characterise the properties of the trap.  However, this is rarely
% the case for complex beams or particles far from equilibrium.
%
% An alternative method for visualising the optical force is to plot a
% colour map or vector field of the optical force as a function of position.
% The following plots such a visualisation (this may take some time
% depending on how many positions are required for the visualisation).

r = linspace(-1, 1, 40)*6*beam.wavelength;
z = linspace(-1, 1, 40)*6*beam.wavelength;

[R, Z] = meshgrid(r, z);
F = beam.force(particle, 'position', {R, 0*R, Z});

% Calculate magnitude of force
Fmag = vecnorm(F, 2, 3);

% Generate colormap/quiver plot showing force
figure();
imagesc(R, Z, Fmag);
hold on;
quiver(R, Z, Fmag(:, :, 1), Fmag(:, :, 3));
hold off;
xlabel('x position [m]');
ylabel('y position [m]');
cb = colorbar();
yaxis(cb, 'Force [N]');

