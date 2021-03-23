% Generate graphical TOC for ott.bsc package
%
% This file is an example from the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../../');

hf = figure();
hf.Position = [680 877 642 221];

rows = 1;
cols = 4;
idx = 1;

%% Simple BSC beams

subplot(rows, cols, idx); idx = idx + 1;

bsc = ott.bsc.Bsc([1;0;0], zeros(3, 1));
beam = ott.beam.BscBeam(bsc);
beam.visFarfieldSphere('type', '3dpolar', ...
    'field', 'E2', 'normalise', true);
title('Bsc');
set(gca(), 'XTick', [], 'YTick', []);

%% Plane wave basis

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.PlaneWave();
beam.visNearfield('axis', 'y', 'field', 'Re(Ex)', 'range', 1e-6*[1,1]);
title('PlaneWave');
set(gca(), 'XTick', [], 'YTick', []);

%% Annular basis

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.Bessel('theta', pi/4);
beam.visNearfield('axis', 'y', 'field', 'Re(Ex)', 'range', 1e-6*[1,1]);
title('Annular');
set(gca(), 'XTick', [], 'YTick', []);

%% Pointmatched beams

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.LaguerreGaussian('lmode', 1, 'pmode', 2);
beam.visNearfield();
title('Pointmatch');
set(gca(), 'XTick', [], 'YTick', []);
