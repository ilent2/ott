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

% Spherical particle, 1um radius
radius = 1.0e-6;
shape = ott.shape.Sphere(radius);

% Particle relative refractive index
index_relative = 1.1;

particle = ott.particle.Fixed.FromShape(shape, index_relative, ...
    'viscosity', 0.001);

% Replace the particle drag tensor with a sphere-wall drag tensor
% Hmm, which one should we choose, and should it change?
%
% Also, this behaves slightly differently??
particle.drag = ott.drag.Stokes... todo

beam = ott.beam.Gaussian.FromNa(1.2, 'index_medium', 1.33, ...
    'wavelength0', 1064e-9, 'power', 0.01);

%% Setup dynamics simulation

temperature = 300.0;      % Temperature [K]

dynamics = ott.tools.Dynamics(beam, particle, ...
    'time_step', 1.0e-4, 'temperature', temperature);

%% Run the simulation

x0 = [0;0;0];
x = dynamics.simulate(x0);

%% Generate visualisations of the results

% TODO: What visualisations do we want?
% How fast does it run, can we look at power spectrums/correlations?

