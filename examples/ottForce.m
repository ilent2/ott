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

%% Setup beam
% The first step is to calculate the incident beam.
% The following creates a Gaussian beam tightly focussed by a microscope
% objective.  The back aperture of the microscope objective creates a
% hard edge, we model this by setting the `truncation_angle` parameter
% of the Gaussian beam class.

% Calculate truncation angle for modelling a microscope objective
NA_obj = 1.2;
index_medium = 1.33;
truncation_angle = asin(NA_obj/index_medium);

% In this example, our Gaussian beam convergence angle is also
% described by an NA value.  This value is larger than the truncation
% angle, corresponding to overfilling.
NA_beam = 1.25;

% All beam parameters have SI units (radians for angles)
beam = ott.beam.Gaussian.FromNa(NA_beam, ...
    'power', 0.01, 'truncation_angle', truncation_angle, ...
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

%% Setup particle
% The next step is to setup the particle.
% For this simulation we use a spherical particle with a 1 micron radius.
% We describe the geometry and let the `FromShape` method choose an
% appropriate T-matrix method (for spheres, this will be `ott.tmatrix.Mie`).

% Create geometry for shape
radius = 1.0e-6;      % Sphere radius [m]
shape = ott.shape.Sphere(radius);

% Setup particle
particle = ott.particle.Fixed.FromShape(shape, ...
    'index_medium', beam.index_medium, ...
    'index_particle', 1.2*beam.index_medium, ...  % Could also use 'index_relative'
    'wavelength0', beam.wavelength0);

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

% Calculate force on the particle along the axis
z = linspace(-1, 1, 60)*6*beam.wavelength;
fz = -beam.force(particle, 'position', [0;0;1].*z);

% Find traps along axis
traps = ott.tools.FindTraps1d.FromArray(z, fz(3, :));

% Get the axial equilibrium (might have no trap, in which case use z=0)
if isempty(traps)
  z0 = 0.0;
  warning('No axial equilibrium found within range');
else
  z0 = traps(1).position;
end

% Calculate radial force
r = linspace(-1, 1, 100)*6*beam.wavelength;
fr = -beam.force(particle, 'position', [1;0;0].*r + [0;0;z0]);

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

% Add equilibrium marker (if found)
if ~isempty(traps)
  hold on;
  traps.plot()
  hold off;
  legend({'X', 'Y', 'Z', 'Center', 'Min/Max'});
end

subplot(1, 2, 2);
plot(r, fr);
xlabel('x position [m]');
ylabel('force [N]');
legend({'X', 'Y', 'Z'});

%% Calculate force using far-field integral
% In some cases, it is also useful to calculate the force by integrating
% the scattered fields in the far-field.  This will often take much longer
% and will be less accurate than the above method, use with caution.
% The following takes about 35 seconds.

f0 = beam.intensityMoment();
fz = zeros(3, numel(z));

tic
for ii = 1:numel(z)
  
  sbeam = beam.scatter(particle, 'position', [0;0;1].*z(ii));
  fz(:, ii) = -(sbeam.intensityMoment() - f0)./beam.speed;
  
end
toc

figure();
plot(z, fz);
xlabel('z position [m]');
ylabel('force [N]');

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

tic
[R, Z] = meshgrid(r, z);
F = -beam.force(particle, 'position', {R, 0*R, Z});
toc

% Calculate magnitude of force
Fmag = squeeze(vecnorm(F, 2));
Fr = squeeze(F(1, :, :));
Fz = squeeze(F(3, :, :));

% Generate colormap/quiver plot showing force
figure();
imagesc(r, z, Fmag);
axis image;
hold on;
quiver(R, Z, Fr, Fz, 'w');
hold off;
xlabel('x position [m]');
ylabel('y position [m]');
cb = colorbar();
ylabel(cb, 'Force [N]');

