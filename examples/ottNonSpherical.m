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

% All three methods agree with this shape (approximately)
shape = ott.shape.Cylinder('radius', 1.0e-6, 'height', 1.5e-6);

% This shape is not nice!!!  Give it a try!
%shape = ott.shape.Cylinder('radius', 1.0e-6, 'height', 1.5e-6);

% Particle relative refractive index
index_relative = 1.1;

%% Calculate T-matrices

tic
Tdda = ott.tmatrix.Dda.FromShape(shape, index_relative);
toc

tic
Tpm = ott.tmatrix.Pointmatch.FromShape(shape, index_relative);
toc

tic
Tebcm = ott.tmatrix.Ebcm.FromShape(shape, index_relative);
toc

%% Calculate force/torque

beam = ott.beam.Gaussian();

z = linspace(-1, 1, 100)*6;
theta = linspace(-pi, pi, 100);
R = ott.utils.roty(theta);

% Calculate force
fzDda = beam.force(Tdda, 'position', [0;0;1].*z);
fzPm = beam.force(Tpm, 'position', [0;0;1].*z);
fzEbcm = beam.force(Tebcm, 'position', [0;0;1].*z);

% Calculate torque
tyDda = beam.torque(Tdda, 'rotation', R);
tyPm = beam.torque(Tpm, 'rotation', R);
tyEbcm = beam.torque(Tebcm, 'rotation', R);

%% Generate a plot

figure();

subplot(1, 2, 1);
plot(z, [fzDda(3, :); fzPm(3, :); fzEbcm(3, :)]);
xlabel('Position [wavelengths]');
ylabel('Force [Q]');
legend({'DDA', 'PM', 'EBCM'});

subplot(1, 2, 1);
plot(theta, [tyDda(2, :); tyPm(2, :); tyEbcm(2, :)]);
xlabel('Angle [radians]');
ylabel('Torque [Q_t]');
legend({'DDA', 'PM', 'EBCM'});

