% Generate graphical TOC for T-matrix methods
%
% This file is an example from the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../../');

figure();

rows = 2;
cols = 4;
idx = 1;

relative_index = 1.2;

%% Sphere (Mie)

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Sphere();
p = shape.surf();
[p.EdgeColor] = deal('none');
camlight; lighting gouraud;
title('Sphere (Mie)');

T = ott.tmatrix.Mie.FromShape(shape, ...
  'relative_index', relative_index);

%% Spheroid (Smarties)

subplot(rows, cols, idx); idx = idx + 1;

a = 1;
b = 1.5;
shape = ott.shape.Ellipsoid([a, a, b]);
p = shape.surf();
[p.EdgeColor] = deal('none');
camlight; lighting gouraud;
title('Spheroid (Smarties)');

T = ott.tmatrix.Smarties.FromShape(shape, ...
  'relative_index', relative_index);

%% Cylinder (Pointmatch, DDA, EBCM)

subplot(rows, cols, idx); idx = idx + 1;

radius = 0.1;
height = 0.2;
shape = ott.shape.Cylinder(radius, height);
shape.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title(['Cylinder' newline '(Pointmatch, Dda, Ebcm)']);

T = ott.tmatrix.Tmatrix.SmartCylinder(shape, ...
    'relative_index', relative_index);

%% Star shaped (Pointmatch)

h1 = subplot(rows, cols, idx); idx = idx + 1;

% Construct a random star shaped particle
theta = linspace(0, pi, 20);
phi = linspace(0, 2*pi, 21);
phi = phi(1:end-1);
[T, P] = meshgrid(theta, phi);
R = 0.1.*randn(size(T)) + 1.0;
R(:, 1) = R(1, 1);
R(:, end) = R(1, end);
T(end+1, :) = T(1, :);
R(end+1, :) = R(1, :);
P(end+1, :) = P(1, :);
[X, Y, Z] = ott.utils.rtp2xyz(R, T, P);
X = reshape(X, size(R));
Y = reshape(Y, size(R));
Z = reshape(Z, size(R));

shape = ott.shape.PatchMesh.FromSurfMatrix(X, Y, Z);
shape.starShaped = true;  % Not computed, must be set explicitly

% Also show another example shape
shape(2) = ott.shape.Cube();
shape(1).position = [-1;-1;0];
shape(2).position = [1;1;0];
shape.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title(['Star shaped' newline '(Pointmatch)']);

% Create T-matrix for only one shape
T = ott.tmatrix.Pointmatch.FromShape(shape(1).*0.1, ...
    'relative_index', relative_index);

%% Rotationally symmetric (EBCM)

h2 = subplot(rows, cols, idx); idx = idx + 1;

% A couple of example shapes
shape = [ott.shape.AxisymFunc.Pill(), ...
        ott.shape.AxisymInterp.Bicone()];
shape(1).position = [-1;-1;0];
shape(2).position = [1;1;0];
shape.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title(['Axially Symmetric' newline '(Ebcm)']);

% Create T-matrix for only one shape
T = ott.tmatrix.Ebcm.FromShape(shape(1), 'relative_index', relative_index);

%% Layered spheres (MieLayered)

h3 = subplot(rows, cols, idx); idx = idx + 1;

shape = [ott.shape.Sphere(0.8), ott.shape.Sphere(1.5), ott.shape.Sphere(2.0)];
p = shape.surf('surfOptions', {'EdgeColor', 'none'});
p(2).FaceAlpha = 0.7;
p(3).FaceAlpha = 0.3;
camlight; lighting gouraud;
title(['Concentric/Layered Spheres' newline '(MieLayered)']);

% Pass in whole array and array of indices
r_indices = [1.5, 1.2, 1.1];
T = ott.tmatrix.MieLayered.FromShape(shape, 'relative_indices', r_indices);

%% Arbitrary shapes (Dda)

h4 = subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.StlLoader('monkey.stl');
shape.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title(['Arbitrary Shape' newline '(Dda)']);

T = ott.tmatrix.Dda.FromShape(shape./shape.maxRadius.*0.1, ...
    'relative_index', relative_index);

%% Tweak positioning

f = gcf();
f.Position(3) = 560;
f.Position(4) = 364;
set(h1, 'Position', [0.7498 0.5826 0.1661 0.3446]);
set(h2, 'Position', [0.1457 0.0549 0.2061 0.4732]);
set(h3, 'Position', [0.4433 0.1298 0.1332 0.3412]);
set(h4, 'Position', [0.7108 0.1254 0.1332 0.3412]);

