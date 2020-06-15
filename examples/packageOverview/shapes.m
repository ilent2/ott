% Generate visualisations of different shape methods

% Add toolbox to path
addpath('../../');

figure('Units', 'pixels', 'Position', [488.2000 176.2000 560 585.6000]);

rows = 5;
cols = 4;
idx = 1;

%% Cube

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.Cube();
shape = shape.rotateYz(pi/4, pi/4);
shape.surf();
title('Cube');

%% Rectangular prism

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.RectangularPrism();
shape = shape.rotateYz(pi/4, pi/4);
shape.surf();
title('RectangularPrism');

%% Sphere

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.Sphere();
shape.surf();
title('Sphere');

%% Cylinder

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.Cylinder();
shape.surf();
title('Cylinder');

%% Ellipsoid

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.Ellipsoid();
shape.surf();
title('Ellipsoid');

%% Superellipsoid

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.Superellipsoid();
shape.surf();
title('Superellipsoid');

%% Plane

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.Plane('normal', [1; 1; 1]);
shape.surf();
title('Plane');

%% Slab

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.Slab('normal', [1; 1; 1]);
shape.surf();
title('Slab');

%% Strata

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.Strata('normal', [1; 1; 1]);
shape.surf();
title('Strata');

%% Empty

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.Empty();
shape = shape.rotateYz(pi/4, pi/4);
shape.surf();
title('Empty');

%% Union

subplot(rows, cols, idx); idx = idx + 1;

shape1 = ott.shapes.Cube();
shape2 = ott.shapes.Cube().rotateYz(pi/4, pi/4);
shape = shape1 | shape2;
shape.surf();
title('Union');

%% Intersection

subplot(rows, cols, idx); idx = idx + 1;

shape1 = ott.shapes.Cube();
shape2 = ott.shapes.Cube().rotateYz(pi/4, pi/4);
shape = shape1 & shape2;
shape.surf();
title('Intersection');

%% Inverse

subplot(rows, cols, idx); idx = idx + 1;

shape1 = ott.shapes.Cube();
shape2 = ott.shapes.Cube().rotateYz(pi/4, pi/4);
shape = shape1 & ~shape2;
shape.surf();
title('Inverse');

%% STL Loader

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.StlLoader('monkey.stl');
shape.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title('StlLoader');

%% OBJ Loader

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.ObjLoader('monkey.obj');
shape = shape.rotateX(0.8.*pi/2).rotateZ(pi/2);
shape.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title('ObjLoader');

%% AxisymFunc: Pill

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.AxisymFunc.Pill();
shape.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title('AxisymFunc.Pill');

%% AxisymFunc: BiconcaveDisc

subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.AxisymFunc.BiconcaveDisc();
shape.surf('surfOptions', {'EdgeColor', 'none'});
camlight; lighting gouraud;
title('AxisymFunc.BiconcaveDisc');

%% AxisymInterp: Bicone

h1 = subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.AxisymInterp.Bicone();
shape.surf();
title('AxisymInterp.Bicone');

%% AxisymInterp: ConeTippedCylinder

h2 = subplot(rows, cols, idx); idx = idx + 1;

shape = ott.shapes.AxisymInterp.ConeTippedCylinder();
shape.surf();
title('AxisymInterp.ConeTippedCylinder');

%% Figure tweaks

set(h1, 'Position', [0.4093 0.1100 0.1389 0.1106]);
set(h2, 'Position', [0.7137 0.1100 0.1389 0.1106]);
