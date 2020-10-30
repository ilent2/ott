% Generate graphical TOC for ott.beam package
%
% This file is an example from the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT was manually installed)
%addpath('../../');

figure();

rows = 3;
cols = 3;
idx = 1;

%% Gaussian

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.Gaussian();
beam.visNearfield();
title('Gaussian');

%% Laguerre-Gaussian

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.LaguerreGaussian('lmode', 1, 'pmode', 2);
beam.visNearfield();
title('Laguerre-Gaussian');

%% Hermite-Gaussian

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.HermiteGaussian('mmode', 1, 'nmode', 2);
beam.visNearfield();
title('Hermite-Gaussian');

%% Ince-Gaussian

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.InceGaussian('lmode', 3, 'porder', 5, 'parity', 'even');
beam.visNearfield();
title('Ince-Gaussian');

%% Plane wave

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.PlaneWave();
beam.visNearfield('axis', 'y', 'field', 'Re(Ex)');
title('Plane Wave');

%% Annular beam (Bessel, Webber, Mathieu)

subplot(rows, cols, idx); idx = idx + 1;

beam = [ott.beam.Bessel('theta', pi/4),
        ott.beam.Webber('theta', pi/4),
        ott.beam.Mathieu('theta', pi/4)];
beam(1).visNearfield('axis', 'y', 'field', 'Re(Ex)');
title(['Annular Beams' newline '(Bessel, Webber, Mathieu)']);

%% Arrays of beams

subplot(rows, cols, idx); idx = idx + 1;

beam = [ott.beam.Gaussian(), ott.beam.LaguerreGaussian(1, 2)];
beam(1).position = [1;1;0];
beam(2).position = [-1;-1;0];
beam = ott.beam.Coherent(beam);  % Can also use ott.beam.Incoherent
beam.visNearfield();
title(['Beams Arrays' newline '(Coherent, Incoherent)']);

%% Paraxial point matching
% Useful for simulating SLMs or beams were we know the far-field

phase = ones(20, 1) .* [linspace(0, 5*pi, 10), linspace(5*pi, 0, 10)];
beam = ott.beam.PmParaxial('phase', phase);
beam.visNearfield();
title('PmParaxial');

subplot(rows, cols, idx); idx = idx + 1;

%% Scattered beams

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.Gaussian();
shape = ott.shape.Sphere();
particle = ott.particle.Particle.FromShape(shape, ...
    'internal', true, 'relative_index', 1.2);
sbeam = particle * beam;
sbeam.visNearfield('axis', 'y');

