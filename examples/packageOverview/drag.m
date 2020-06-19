% Graphical overview of the drag package
% Generates two figures: the graphical TOC and the wall effect comparison

% Add toolbox to path
addpath('../../');

figure();

rows = 2;
cols = 3;
idx = 1;

%% Sphere

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.Sphere();
p = shape.surf();
[p.EdgeColor] = deal('none');
camlight; lighting gouraud;
title('StokesSphere');

dragTensor = ott.drag.Stokes.FromShape(shape);

%% Sphere near wall

subplot(rows, cols, idx); idx = idx + 1;

shape = [ott.shapes.Sphere(), ott.shapes.Plane('position', [0; 0; -1.2])];
p(1) = shape(1).surf();
hold on;
p(2) = shape(2).surf('scale', 3);
hold off;
[p.EdgeColor] = deal('none');
camlight; lighting gouraud;
title(['Sphere/Plane Wall Effects', newline, ...
  '(FaxenSphere, PadeSphere,' newline 'ChaouiSphere)']);

dragTensor = ott.drag.Stokes.FromShape(shape);

%% Cylinder

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.Cylinder(0.5, 3);
p = shape.surf();
[p.EdgeColor] = deal('none');
camlight; lighting gouraud;
title('StokesCylinder');

dragTensor = ott.drag.Stokes.FromShape(shape);

%% Star shaped particle

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

shape = ott.shapes.PatchMesh.FromSurfMatrix(X, Y, Z);
shape.starShaped = true;  % Not computed, must be set explicitly
p = shape.surf();
[p.EdgeColor] = deal('none');
camlight; lighting gouraud;
title(['Star-shaped particle' newline, ...
  '(StokesLambNn, StokesLambPm)']);

dragTensor = ott.drag.Stokes.FromShape(shape);

%% Eccentric spherers

h2 = subplot(rows, cols, idx); idx = idx + 1;

shape = [ott.shapes.Sphere(1.0), ott.shapes.Sphere(2.0)];
shape(1).position = [0;0;-0.5];
p = shape.surf();
p(2).FaceAlpha = 0.5;
[p.EdgeColor] = deal('none');
camlight; lighting gouraud;
title('EccentricSpheresNn');

dragTensor = ott.drag.Stokes.FromShape(shape);

%% Tweak positioning

set(h1, 'Position', [0.2264 0.1100 0.2067 0.3049]);
set(h2, 'Position', [0.5733 0.1100 0.2067 0.3049]);
