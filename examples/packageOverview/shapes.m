% Generate visualisations of different shape methods
%
% This file is an example from the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT was manually installed)
%addpath('../../');

f = figure();
f.Units = 'pixels';
f.Position(2) = f.Position(2) - (585 - f.Position(4));
f.Position(3) = 560;
f.Position(4) = 585.6000;

rows = 5;
cols = 4;
idx = 1;

%% Cube

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Cube();
shape = shape.rotateYz(pi/4, pi/4);
shape.surf();
title('Cube');

%% Rectangular prism

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.RectangularPrism();
shape = shape.rotateYz(pi/4, pi/4);
shape.surf();
title('RectangularPrism');

%% Sphere

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Sphere();
shape.surf();
title('Sphere');

%% Cylinder

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Cylinder();
shape.surf();
title('Cylinder');

%% Ellipsoid

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Ellipsoid();
shape.surf();
title('Ellipsoid');

%% Superellipsoid

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Superellipsoid();
shape.surf();
title('Superellipsoid');

%% Plane

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Plane('normal', [1; 1; 1]);
shape.surf();
title('Plane');

%% Slab

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Slab('normal', [1; 1; 1]);
shape.surf();
title('Slab');

%% Strata

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Strata('normal', [1; 1; 1]);
shape.surf();
title('Strata');

%% Empty

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.Empty();
shape = shape.rotateYz(pi/4, pi/4);
shape.surf();
title('Empty');

%% Union

subplot(rows, cols, idx); idx = idx + 1;

shape1 = ott.shape.Cube();
shape2 = ott.shape.Cube().rotateYz(pi/4, pi/4);
shape = shape1 | shape2;
shape.surf();
camlight; lighting gouraud;
title('Union');

%% Intersection

subplot(rows, cols, idx); idx = idx + 1;

shape1 = ott.shape.Cube();
shape2 = ott.shape.Cube().rotateYz(pi/4, pi/4);
shape = shape1 & shape2;
shape.surf();
camlight; lighting gouraud;
title('Intersection');

%% Inverse

subplot(rows, cols, idx); idx = idx + 1;

shape1 = ott.shape.Cube();
shape2 = ott.shape.Cube().rotateYz(pi/4, pi/4);
shape = shape1 & ~shape2;
shape.surf();
camlight; lighting gouraud;
title('Inverse');

%% STL Loader

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.StlLoader('monkey.stl');
shape.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title('StlLoader');

%% OBJ Loader

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.ObjLoader('monkey.obj');
shape = shape.rotateX(0.8.*pi/2).rotateZ(pi/2);
shape.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title('ObjLoader');

%% AxisymFunc: Pill

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.AxisymFunc.Pill();
shape.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title('AxisymFunc.Pill');

%% AxisymFunc: BiconcaveDisc

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.AxisymFunc.BiconcaveDisc();
shape.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title('AxisymFunc.BiconcaveDisc');

%% AxisymInterp: Bicone

h1 = subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.AxisymInterp.Bicone();
shape.surf();
title('AxisymInterp.Bicone');

%% AxisymInterp: ConeTippedCylinder

h2 = subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shape.AxisymInterp.ConeTippedCylinder();
shape.surf();
title('AxisymInterp.ConeTippedCylinder');

%% Figure tweaks

set(h1, 'Position', [0.4093 0.1100 0.1389 0.1106]);
set(h2, 'Position', [0.7137 0.1100 0.1389 0.1106]);
