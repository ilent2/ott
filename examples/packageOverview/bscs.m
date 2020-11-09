% Generate graphical TOC for ott.bsc package
%
% This file is an example from the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../../');

figure();

rows = 1;
cols = 4;
idx = 1;

%% Simple BSC beams

bsc = ott.bsc.Bsc([1;0;0], 0);
beam = ott.beam.BscBeam(bsc);
beam.visFarfieldSphere('type', '3dpolar');
title('Bsc');

%% Plane wave basis

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.PlaneWave();
beam.visNearfield('axis', 'y', 'field', 'Re(Ex)');
title('PlaneWave');

%% Annular basis

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.Bessel('theta', pi/4);
beam.visNearfield('axis', 'y', 'field', 'Re(Ex)');
title('Annular');

%% Pointmatched beams

subplot(rows, cols, idx); idx = idx + 1;

beam = ott.beam.LaguerreGaussian('lmode', 1, 'pmode', 2);
beam.visNearfield();
title('Pointmatch');
