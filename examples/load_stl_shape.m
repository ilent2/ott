% Load a STL file, generate a voxel grid and scatter a plane wave off the shape

% Add ott to the path
addpath('../');

% Setup the figure for output
figure();

%% Load the STL object

% STL reader only supports binary STL files for now
% Change the filename to an STL file on your computer or download
% an STL file from online.
shape = ott.shapes.StlLoader('mesh.stl');

% Visualise the shape
h = subplot(1, 3, 1);
shape.surf('surfoptions', {'EdgeColor', 'none'});
axis equal;
view([-30, 25]);
title('Mesh');
xlabel('X');
ylabel('Y');
zlabel('Z');

%% Convert the shape to a voxel grid

spacing = 1.0;  % Adjust this for your mesh

subplot(1, 3, 2);
xyz = shape.voxels(spacing, 'visualise', true);
axis equal;
view([-30, 25]);
title('Voxels');
xlabel('X');
ylabel('Y');
zlabel('Z');

%% Calculate the T-matrix for this shape

scale = 0.1;   % Convert from mesh units to simulation units

Tmatrix = ott.TmatrixDda(scale*xyz, 'index_particle', 1.4, 'index_medium', 1.33, ...
  'wavelength0', 1.0);

%% Calculate torque as a function of axial angle

angles = linspace(0, 2*pi, 100);
rotation = ott.utils.roty(angles*180/pi);

beam = ott.BscPmGauss('NA', 1.02, 'power', 1.0, 'index_medium', 1.33, ...
  'polarisation', [1, 1i], 'wavelength0', 1.0);

[~, torque] = ott.forcetorque(beam, Tmatrix, 'rotation', rotation);

subplot(1, 3, 3);
plot(angles, torque);
xlabel('Y-Angle [rad]');
ylabel('Torque [a.u.]');
legend({'X', 'Y', 'Z'}, 'Location', 'SouthEast');
title('Torque');

%% Animate the mesh rotating (no translation)

beam = ott.BscPmGauss('NA', 1.02, 'power', 1.0, 'index_medium', 1.33, ...
  'polarisation', [1, 1i], 'wavelength0', 1.0);

x = [0;0;0];
Rx = eye(3);
dtheta = 2*pi/100;
dx = 0.0;  % no translation

figure();
shape.surf('rotation', Rx, 'surfoptions', {'EdgeColor', 'none'});
axis equal;
axis([-10, 10, -10, 10, -2, 10]);

for ii = 1:100
  [f, t] = ott.forcetorque(beam, Tmatrix, 'rotation', Rx, 'position', x);
  x = x + dx*f./vecnorm(f);
  Rx = ott.utils.rotation_matrix(dtheta*t./vecnorm(t))*Rx;
  
  shape.surf('rotation', Rx, 'surfoptions', {'EdgeColor', 'none'});
  axis equal;
  axis([-10, 10, -10, 10, -2, 10]);
  drawnow;
end
