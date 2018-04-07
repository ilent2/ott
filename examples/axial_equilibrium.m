% ott.axial_equilibrium can be used to calculate the trapping position
% and spring constant along the z axis.  If the beam is rotated, the
% method can also be used for other axies too.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add the toolbox to the path (assuming we are in ott/examples)
addpath('../');

% Make warnings less obtrusive
ott_warning('once');
change_warnings('off');

% Specify refractive indices
n_medium = 1.34;
n_particle = 1.59;

% Specify the wavelength in freespace [m]
wavelength = 1064.0e-9;

% Specify the particle radius (sphere)
radius = 1.0*wavelength/n_medium;

%% Calculate the beam

beam = ott.BscPmGauss('angle_deg', 50, ...
    'polarisation', [ 1 0 ], 'power', 1.0);

%% Calculate the particle T-matrix

T = ott.Tmatrix.simple('sphere', radius, ...
    'n_medium', n_medium, ...
    'n_particle', n_particle, ...
    'wavelength0', wavelength);

%% Find the equilibrium and trap stiffness in the x and x directions

% Find the equilibrium in the z-direction
[z,kz] = ott.axial_equilibrium(T, beam)

% Translate the beam to the z-axis equilibrium
beam = beam.translateZ(z);

% Rotate the beam about the y axis (so the beam is aligned with the x axis)
beam = beam.rotateY(pi/2.0);

% Calculate the equilibrium in the x-direction
[x,kx] = ott.axial_equilibrium(T, beam);

