% Example calculation of forces on a sphere in a optical beam
%
% This example includes a couple of example beams (Gaussian, LG and HG),
% combining the separate example scripts found in ottv1.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add the toolbox to the path (assuming we are in ott/examples)
addpath('../');

%% Describe the particle, beam and surrounding medium

% Make warnings less obtrusive
ott_warning('once');
change_warnings('off');

% Refractive index of particle and medium
n_medium = 1.33;
n_particle = 1.59;

% Wavelength of light in vacuum [m]
wavelength0 = 1064e-9;

% Calculate the wavelength in the medium
wavelength_medium = wavelength0/n_medium;

% Radius of particle
radius = 1.0*wavelength_medium;

% Set the beam type (must be one of lg, gaussian or hg)
% This is used bellow to select from the beam presets
beam_type = 'gaussian';

% Specify the numerical aparture of the beam/microscope
NA = 1.02;

%% Setup the T-matrix for the particle

% Create a T-matrix for a sphere
T = ott.Tmatrix.simple('sphere', radius, 'wavelength0', wavelength0, ...
    'index_medium', n_medium, 'index_particle', n_particle);

%% Setup the T-matrix for the beam

switch beam_type
  case 'gaussian'

    % Create a simple Gaussian beam with circular polarisation
    beam = ott.BscPmGauss('NA', NA, 'polarisation', [ 1 1i ], ...
        'index_medium', n_medium, 'wavelength0', wavelength0);

  case 'lg'

    % Create a LG03 beam with circular polarisation
    %
    % We could have specified the NA as before but this time we
    % calculate the beam angle ourselves and specify that
    beam_angle = asin(NA/n_medium);
    beam = ott.BscPmGauss('lg', [ 0 3 ], ...
        'polarisation', polarisation, 'angle', beam_angle, ...
        'index_medium', n_medium, 'wavelength0', wavelength0);

  case 'hg'

    % Create a HG23 beam with circular polarisation
    beam_angle = asin(NA/n_medium);
    beam = ott.BscPmGauss('hg', [ 2 3 ], ...
        'polarisation', polarisation, 'angle', beam_angle, ...
        'index_medium', n_medium, 'wavelength0', wavelength0);

  otherwise
    error('Unsupported beam type');
end

% Normalize the beam power
beam = beam / beam.power();

%% Generate force/position graphs

% Time the calculation
tic

%calculate the force along z
z = linspace(-8,8,80)*wavelength;
fz = zeros(3, length(z));
for nz = 1:length(z)
  tbeam = beam.translateZ(z(nz));
  sbeam = T * tbeam;
  fz(:, nz) = ott.forcetorque(tbeam, sbeam);
end

% Find the equilibrium along the z axis
zeq = ott.find_equilibrium(z, fz(3, :));
if isempty(zeq)
  warning('No axial equilibrium in range!')
  zeq=0;
end
zeq = zeq(1);

% Calculate force along x-axis (with z = zeq, if found)
r = linspace(-4,4,80)*wavelength;
fr = zeros(3, length(r));
for nr = 1:length(r)
  tbeam = beam.translateXyz([r(nr), 0.0, zeq]);
  sbeam = T * tbeam;
  fr(:, nr) = ott.forcetorque(tbeam, sbeam);
end

% Finish timing the calculation
toc

% Generate the plots

figure(1); plot(z/wavelength,fz(3, :));
xlabel('{\it z} (x\lambda)');
ylabel('{\it Q_z}');
aa = axis;
hold on;
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
hold off;

figure(2); plot(r/wavelength,fr(1, :));
xlabel('{\it r} (x\lambda)');
ylabel('{\it Q_r}');
aa = axis;
hold on;
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
hold off;
