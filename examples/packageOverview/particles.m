% Generate graphical TOC for particle representations
%
% This file is an example from the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../../');

hf = figure();
hf.Position(4) = 247.2000;

rows = 1;
cols = 2;
idx = 1;

%% Fixed particle

subplot(rows, cols, idx); idx = idx + 1;
wavelength0 = 1064e-9;

shape = ott.shape.Cube(0.1e-6);
particle = ott.particle.Fixed.FromShape(shape, ...
  'index_relative', 1.2, 'wavelength0', wavelength0);
particle.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
axis([-1,1, -1,1, -1,1]*1.5);
view([27, 11]);
title('Fixed');

%% Variable particle

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Sphere(1e-6);
particle = ott.particle.Variable.FromShape('initial_shape', shape);

% Generate visualisation of shape at current size
particle.position = [-2;0;0]*1e-6;
particle.index_relative = 1.2;
particle.surf('surfOptions', {'EdgeColor', 'none'});

% Generate visualisation of shape at new size/ri
hold on;
particle.position = [2;0;0]*1e-6;
particle.index_relative = 1.5;
particle.shape.radius = 1.5e-6;
particle.surf('surfOptions', {'EdgeColor', 'none', 'FaceColor', 'blue'});
hold off;
view([0, 40]);

camlight; lighting gouraud;
title('Variable');

%% Figure formatting

% Add an arrow
annotation(hf,'arrow',[0.703214285714285 0.751785714285714],...
  [0.516619047619048 0.517619047619048]);
