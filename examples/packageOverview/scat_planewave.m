% Demonstration of scat.planewave methods

figure();

rows = 1;
cols = 3;
idx = 1;

%% Plane wave

subplot(rows, cols, idx); idx = idx + 1;

medium = ott.beam.medium.Vacuum.Unitary;
beam = ott.beam.PlaneWave('medium', medium);

beam.visualise();

%% scat.planewave.Plane method

subplot(rows, cols, idx); idx = idx + 1;

plane = ott.shapes.Plane();
media = Dielectric.FromIndex(1.2)./beam.medium;

particle = ott.scat.planewave.Plane(plane, media);

[rbeam, tbeam] = beam.scatter(particle);

% Change reflected beam to total-field beam
rbeam.type = 'total';

totBeam.visualise('mask', ~plane);
hold on;
tbeam.visualise('mask', plane);
hold off;

%% scat.planewave.Strata method

subplot(rows, cols, idx); idx = idx + 1;

strata = ott.shapes.Strata([0, 0.3, 0.8, 5]);
media = [Dielectric.FromIndex(1.2), ...
  Dielectric.FromIndex(1.5), Dielectric.FromIndex(1.2), ...
  Dielectric.FromIndex(1.0)]./beam.medium;

particle = ott.scat.planewave.Strata(strata, media);
[rbeam, tbeam] = beam.scatter(particle);

% Change reflected beam to total-field beam
rbeam.type = 'total';

totBeam.visualise('mask', ~plane);
hold on;
tbeam.visualise('mask', plane);
hold off;


