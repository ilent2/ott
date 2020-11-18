% ottWallEffect.m -- Simulate a particle near a wall
%
% This simulation shows how a spherical particle near a plane wall
% could be simulated.  There is a fluid interaction between the particle
% and the wall but no optical interaction (i.e., a transparent wall).
%
% This example is almost exactly the same as ottDynamics.m except for the
% drag tensor (which is a sphere-wall drag tensor) and the visualisations.
%
% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../');

%% Setup beam and particle

% Spherical particle, 1um radius, 1um above a wall
radius = 0.5e-6;
shape = ott.shape.Sphere(radius);
plane = ott.shape.Plane();
plane.position = [0;0;-2*radius];
geometry = [shape, plane];

% Particle relative refractive index
index_relative = 1.1;

% Setup beam
beam = ott.beam.Gaussian.FromNa(1.2, 'index_medium', 1.33, ...
    'wavelength0', 1064e-9, 'power', 0.01);

% Setup particle: sphere T-matrix + plane/sphere geometry
particle = ott.particle.Fixed();
particle.tmatrix = ott.tmatrix.Tmatrix.FromShape(shape./beam.wavelength, ...
  'index_relative', index_relative);
particle.shape = geometry;

%% Setup dynamics simulation

temperature = 300.0;      % Temperature [K]

dragMethod = @(shape) ott.drag.StokesSphereWall.FromShape(shape, ...
  'viscosity', 0.001);

dynamics = ott.dynamics.WallEffect('beam', beam, 'particle', particle, ...
    'timeStep', 1.0e-4, 'temperature', temperature, ...
    'dragMethod', dragMethod);

%% Run the simulation

% Setup a figure to show the progress
figure();
ax = axes();
axis([-1,1,-1,1,-1,1]*1.5e-6);

x0 = [0;0;0];
totalTime = 1e-2;   % Simulation time [s]
outputRate = 0.1;   % How frequently the figure is updated [s]
[t, x] = dynamics.simulate(totalTime, 'position', x0, ...
  'plot_axes', ax, 'outputRate', outputRate);

%% Generate 2-D histogram of position

zbins = linspace(-1, 1, 20)*0.1e-6;
xbins = linspace(-1, 1, 20)*0.1e-6;

figure();
counts = histcounts2(x(3, :), x(1, :), zbins, xbins);
imagesc(zbins, xbins, counts);
xlabel('z [m]');
ylabel('x [m]');
axis image;
