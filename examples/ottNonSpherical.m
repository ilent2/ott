% ottNonSpherical.m -- Calculate force/torque on a non-spherical particle
%
% This example shows how to generate T-matrices for non-spherical particles
% and calculate the force/torques.  Care should be taken when using any
% of the T-matrix methods with non-spherical particles, not all methods
% work in all circumstances (as shown bellow).
%
% Copyright 2020 Isaac Lenton.
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../');

%% Setup simulation parameters

% Particle relative refractive index
index_medium = 1.0;
index_relative = 1.2;

wavelength0 = 1.064e-6;
wavelength = wavelength0/index_medium;

scale = 0.1;

% DDA and EBCM approximately agree
shape = ott.shape.Cylinder('radius', 0.4*wavelength0, 'height', 0.8*wavelength0)*scale;

% This shape is not nice!!!  Give it a try!
% shape = ott.shape.Cylinder('radius', 0.4*wavelength0, 'height', 0.5*wavelength0)*scale;

%% Calculate T-matrices

% Make memory warnings errors to make bailing easier
warning('error', 'ott:tmatrix:dda:Dipole:memory_may_exceed_physical');

Nmax = 10;  % Used for DDA/EBCM, doesnt work for PM
NmaxPm = 4;

tic
% DDA may give better results by choosing a smaller dipole spacing.
Tdda = ott.tmatrix.Dda.FromShape(shape./wavelength, ...
  'index_relative', index_relative, 'spacing', 0.005, 'Nmax', Nmax);
Pdda = ott.particle.Fixed('tmatrix', Tdda, 'shape', shape);
toc

tic
Tpm = ott.tmatrix.Pointmatch.FromShape(shape./wavelength, ...
  'index_relative', index_relative, 'Nmax', NmaxPm);
Ppm = ott.particle.Fixed('tmatrix', Tpm, 'shape', shape);
toc

tic
Tebcm = ott.tmatrix.Ebcm.FromShape(shape./wavelength, ...
  'index_relative', index_relative, 'Nmax', Nmax);
Pebcm = ott.particle.Fixed('tmatrix', Tebcm, 'shape', shape);
toc

%% Calculate force/torque

beam = ott.beam.Gaussian(0.6e-6, ...
  'wavelength0', wavelength0);

z = linspace(-1, 1, 100)*4*wavelength;
theta = linspace(-pi, pi, 100);
R = ott.utils.roty(theta*180/pi);

% Calculate force
Z = [0;0;1].*z;
fzDda = -beam.force(Pdda, 'position', Z);
fzPm = -beam.force(Ppm, 'position', Z);
fzEbcm = -beam.force(Pebcm, 'position', Z);

% Calculate torque
tyDda = -beam.torque(Pdda, 'rotation', R);
tyPm = -beam.torque(Ppm, 'rotation', R);
tyEbcm = -beam.torque(Pebcm, 'rotation', R);

%% Generate a plot

figure();

subplot(1, 2, 1);
plot(z, [fzDda(3, :); fzPm(3, :); fzEbcm(3, :)]);
xlabel('Position [wavelengths]');
ylabel('Force [Q]');
legend({'DDA', 'PM', 'EBCM'});

subplot(1, 2, 2);
plot(theta, [tyDda(2, :); tyPm(2, :); tyEbcm(2, :)]);
xlabel('Angle [radians]');
ylabel('Torque [Q_t]');
legend({'DDA', 'PM', 'EBCM'});

