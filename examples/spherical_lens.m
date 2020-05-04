% Demonstrates how the spherical lens class can be used

% Add OTT to the path
addpath('../');

% Cteate some lenses
lens = ott.scat.geometric.SphericalLens(3.0, 'position', [0;0;0]);

% Set-up an incident beam (parallel rays)
x = linspace(-1, 1, 5);
[xx, yy, zz] = meshgrid(x, x, -1);
xyz = [xx(:), yy(:), zz(:)].';
beam1 = ott.beam.Ray('origin', xyz, 'direction', [0;0.1;1]);

% Calculate rays scattered by lenses
beam2 = lens.scatter(beam1);

%% Generate a visualisation

figure();
lens.surf('surfoptions', {'FaceAlpha', 0.2});
hold on;
beam2.visualiseAllRays('show_polarisation', false);
hold off;
