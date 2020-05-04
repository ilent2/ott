% Demonstrates how the Thin lens class can be used

% Add OTT to the path
addpath('../');

% Cteate some lenses
L1 = ott.scat.geometric.ThinLens(1.0, 'position', [0;0;0]);
L2 = ott.scat.geometric.ThinLens(-0.8, 'position', [0;0;0.2]);

% Set-up an incident beam (parallel rays)
x = linspace(-1, 1, 20);
[xx, yy, zz] = meshgrid(x, x, -1);
xyz = [xx(:), yy(:), zz(:)].';
beam1 = ott.beam.Ray('origin', xyz, 'direction', [0;0;1]);

% Calculate rays scattered by lenses
beam2 = L1.scatter(beam1);
beam3 = L2.scatter(beam2);

%% Generate a visualisation

figure();
L1.surf('FaceAlpha', 0.2);
hold on;
L2.surf('FaceAlpha', 0.2);
beam3.visualiseAllRays();
hold off;
