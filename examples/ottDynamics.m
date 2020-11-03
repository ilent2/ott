% ottDynamics.m -- Dynamics simulation of an optically trapped particle
%
% This simulation shows how the toolbox can be used to simulate the
% dynamics of an optically trapped particle.  Accuracy of these results
% depends on the accuracy of the optical scattering calculation, the drag
% calculation, and the simulation properties (most importantly, time step
% size).
%
% The example bellow shows a simulation of a spherical particle in a
% Gaussian beam.  This should be fairly accurate.  If you replace the
% beam or particle, it is important to check that the drag, T-matrix and
% beam are still what you would expect.  You may also need to reduce the
% simulation step size if the forces change significantly.
%
% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../');

%% Setup beam and particle

% Spherical particle, 1um radius
radius = 1.0e-6;
shape = ott.shape.Sphere(radius);

% Particle relative refractive index
index_relative = 1.1;

particle = ott.particle.Fixed.FromShape(shape, index_relative, ...
    'viscosity', 0.001);

beam = ott.beam.Gaussian.FromNa(1.2, 'index_medium', 1.33, ...
    'wavelength0', 1064e-9, 'power', 0.01);

%% Setup dynamics simulation

temperature = 300.0;      % Temperature [K]

dynamics = ott.tools.Dynamics(beam, particle, ...
    'time_step', 1.0e-4, 'temperature', temperature);

%% Run the simulation

x0 = [0;0;0];
x = dynamics.simulate(x0);

%% Generate a plot of the results

figure();
plot3(x(1, :), x(2, :), x(3, :));
xlabel('x [m]');
xlabel('y [m]');
xlabel('z [m]');
daspect([1, 1, 1]);

%% Generate 2-D histogram of position

xbins = linspace(-1, 1, 100)*1e-6;
ybins = linspace(-1, 1, 100)*1e-6;

figure();
histogram2(x(1, :), x(2, :), xbins, ybins);
xlabel('X [m]');
xlabel('Y [m]');

