% Generate graphical TOC for particle representations
%
% This file is an example from the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT was manually installed)
%addpath('../../');

figure();

rows = 1;
cols = 2;
idx = 1;

%% Fixed particle

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Cube();
particle = ott.particle.Fixed.FromShape(shape);
particle.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title('Fixed');

%% Variable particle

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Sphere();
particle = ott.particle.Variable.FromShape(shape);

% Generate visualisation of shape at current size
particle.position = [-1;0;0];
particle.surf('surfOptions', {'EdgeColor', 'none'});

% Generate visualisation of shape at new size/ri
particle.position = [1;0;0];
particle.relative_index = 1.2;
particle.shape.radius = 1.5;
particle.surf('surfOptions', {'EdgeColor', 'none'});

camlight; lighting gouraud;
title('Variable');

