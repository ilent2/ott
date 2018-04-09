% Example showing how multiple beams can be added together
%
% There are a couple of ways to add beams together.
%
%   * Incoherently: the beam coefficients are not added together,
%       instead the force is calculated separately for each beam.
%
%   * Coherently with larger Nmax: the beam is calculated at the
%       origin, Nmax is then expanded to include the translation
%       distance and the beam coefficients are added.  The resulting
%       beam can then be translated inside the new Nmax region.
%
%   * Coherently with same Nmax: for each force calculation the
%       beams are translated, the translated beam coefficients are
%       added and the force calculated from the combined beam.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add the toolbox to the path (assuming we are in ott/examples)
addpath('../');

import ott.*
import ott.utils.*

% Make warnings less obtrusive
ott_warning('once');
change_warnings('off');

%% Create a particle to scatter the beams with

% Wavelength in medium/vacuum [m]
wavelength = 1064.0e-9;

T = ott.Tmatrix.simple('sphere', wavelength, 'index_medium', 1.0, ...
    'index_particle', 1.2, 'wavelength0', wavelength);

%% Create a simple gaussian beam
% We will create displaced copies of this beam

beam = ott.BscPmGauss('polarisation', [1 1i], 'angle_deg', 50, ...
    'index_medium', 1.0, 'wavelength0', wavelength);
beam = beam / beam.power();

% Displacement of beams [wavelength_medium]
displacement = 1.2*wavelength;

% Phase shift for 2nd beam
phase = exp(2*pi*1i/0.1);

% Range for force/displacement graph
x = linspace(-8, 8, 80)*wavelength;

%% Coherent beams with expanded Nmax

% Calculate new Nmax
Nmax = ott.utils.ka2nmax(ott.utils.nmax2ka(beam.Nmax) ...
    + displacement*T.k_medium);

% Change the Nmax and create the two beams
beam1 = beam;
beam1.Nmax = Nmax;
beam2 = beam1.translateXyz(displacement, 0, 0);
beam1 = beam1.translateXyz(-displacement, 0, 0);

% Add the beams
nbeam = beam1 + beam2 * phase;

% Calculate the force along the x-axis
fx1 = zeros(3, length(x));
for ii = 1:length(x)
  tbeam = nbeam.translateXyz(x(ii), 0, 0);
  sbeam = T * tbeam;
  fx1(:, ii) = ott.forcetorque(tbeam, sbeam);
end

%% Coherent beams with same Nmax

fx2 = zeros(3, length(x));

for ii = 1:length(x)

  % Translate and add the beams
  beam1 = beam.translateXyz(x(ii)+displacement, 0, 0);
  beam2 = beam.translateXyz(x(ii)-displacement, 0, 0);
  tbeam = beam1 + beam2 * phase;

  % Scatter the beam and calculate the force
  sbeam = T * tbeam;
  fx2(:, ii) = ott.forcetorque(tbeam, sbeam);
end

%% Incoherent beams

fx3 = zeros(3, length(x));

for ii = 1:length(x)

  % Translate beams
  beam1 = beam.translateXyz(x(ii)+displacement, 0, 0);
  beam2 = beam.translateXyz(x(ii)-displacement, 0, 0) * phase;

  % Scatter the beams
  sbeam1 = T * beam1;
  sbeam2 = T * beam2;

  % Calculate the force
  fx3(:, ii) = ott.forcetorque(beam1, sbeam1) ...
      + ott.forcetorque(beam2, sbeam2);

>>>>>>> Started work on the 1.4.0 examples
end

%% Generate a figure showing the force displacement graphs

x = x/wavelength;

figure(1);
plot(x, fx1(1, :), 'b', x, fx2(1, :), 'r*', x, fx3(1, :));
legend('Coherent (large Nmax)', 'Coherent (small Nmax)', 'Incoherent');
xlabel('x [\lambda]')
ylabel('Q_x')
title('Force displacement curves for sphere in multiple beams')

