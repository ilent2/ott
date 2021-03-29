% Generate graphical TOC for ott.beam package
%
% This file is an example from the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../../');

figure('position', [680 544 444 554]);

rows = 3;
cols = 3;
idx = 1;

%% Gaussian

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.Gaussian();
beam.visNearfield();
title('Gaussian');
set(gca(), 'XTick', [], 'YTick', []);

%% Laguerre-Gaussian

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.LaguerreGaussian('lmode', 1, 'pmode', 2);
beam.visNearfield();
title('Laguerre-Gaussian');
set(gca(), 'XTick', [], 'YTick', []);

%% Hermite-Gaussian

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.HermiteGaussian('mmode', 1, 'nmode', 2);
beam.visNearfield();
title('Hermite-Gaussian');
set(gca(), 'XTick', [], 'YTick', []);

%% Ince-Gaussian

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.InceGaussian('lmode', 3, 'porder', 5, 'parity', 'even');
beam.visNearfield();
title('Ince-Gaussian');
set(gca(), 'XTick', [], 'YTick', []);

%% Plane wave

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.PlaneWave();
beam.visNearfield('axis', 'y', 'range', [1,1]*2e-6, 'field', 'Re(Ex)');
title('Plane Wave');
set(gca(), 'XTick', [], 'YTick', []);

%% Annular beam (Bessel, Webber, Mathieu)

subplot(rows, cols, idx); idx = idx + 1;

beam = [ott.beam.Bessel('theta', pi/4), ...
        ott.beam.Webber('theta', pi/4), ...
        ott.beam.Mathieu('theta', pi/4)];
beam(1).visNearfield('axis', 'y', 'range', [1,1]*2e-6, 'field', 'Re(Ex)');
title(['Annular Beams' newline '(Bessel, Webber, Mathieu)']);
set(gca(), 'XTick', [], 'YTick', []);

%% Arrays of beams

subplot(rows, cols, idx); idx = idx + 1;

beam = [ott.beam.Gaussian(), ott.beam.LaguerreGaussian(...
    'lmode', 1, 'pmode', 2)];
beam(1).position = [1;1;0]*beam(1).wavelength;
beam(2).position = [-1;-1;0]*beam(1).wavelength;
beam = ott.beam.Array(beam, 'arrayType', 'coherent');
beam.visNearfield();
title(['Beams Arrays' newline '(Coherent/Incoherent)']);
set(gca(), 'XTick', [], 'YTick', []);

%% Paraxial point matching
% Useful for simulating SLMs or beams were we know the far-field

subplot(rows, cols, idx); idx = idx + 1;

% Generate coordinates for pattern
x = linspace(-1, 1, 20);
y = linspace(-1, 1, 20);
[X, Y] = ndgrid(x, y);
P = atan2(X, Y);

% Calculate incident field
E0 = exp(-(X.^2 + Y.^2)./4);

% Calculate SLM-like pattern
kx = 2;
phi = 2*pi*kx*X + 2*P;

% Calculate field at back aperture
E = E0 .* exp(1i*phi);

% Calculate beam
beam = ott.beam.PmParaxial.InterpProfile(X, Y, E, 'Nmax', 20);
beam.visNearfield('range', [1,1]*6e-6);
title('PmParaxial');
set(gca(), 'XTick', [], 'YTick', []);

%% Scattered beams

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.Gaussian();
shape = ott.shape.Sphere(1.2*beam.wavelength0);
particle = ott.particle.Fixed.FromShape(shape, ...
    'internal', true, 'index_relative', 1.2, ...
    'wavelength0', beam.wavelength0);
sbeam = beam.scatter(particle);
sbeam.visNearfield('axis', 'y');
title('Scattered beam');
set(gca(), 'XTick', [], 'YTick', []);

