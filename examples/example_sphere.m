% Example calculation of forces on a sphere in a optical beam
%
% This example includes a couple of example beams (Gaussian, LG and HG),
% combining the separate example scripts found in ottv1.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add the toolbox to the path (assuming we are in ott/examples)
addpath('../');

% Close open figures
close all;

%% Describe the particle, beam and surrounding medium

% Make warnings less obtrusive
ott.warning('once');
ott.change_warnings('off');

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

tic

% Create a T-matrix for a sphere
T = ott.Tmatrix.simple('sphere', radius, 'wavelength0', wavelength0, ...
    'index_medium', n_medium, 'index_particle', n_particle);

disp(['T-matrix calculation took ' num2str(toc) ' seconds']);

%% Calculate the beam shape coefficients

tic

switch beam_type
  case 'gaussian'

    % Create a simple Gaussian beam with circular polarisation
    % We could also calculate the beam angle ourseves and specify that.
    beam = ott.BscPmGauss('NA', NA, 'polarisation', [ 1 1i ], ...
        'index_medium', n_medium, 'wavelength0', wavelength0);

  case 'lg'

    % Create a LG03 beam with circular polarisation
    beam = ott.BscPmGauss('lg', [ 0 3 ], ...
        'polarisation', [ 1 1i ], 'NA', NA, ...
        'index_medium', n_medium, 'wavelength0', wavelength0);

  case 'hg'

    % Create a HG23 beam with circular polarisation
    beam = ott.BscPmGauss('hg', [ 2 3 ], ...
        'polarisation', [ 1 1i ], 'NA', NA, ...
        'index_medium', n_medium, 'wavelength0', wavelength0);

  otherwise
    error('Unsupported beam type');
end

% Normalize the beam power
beam.power = 1.0;

disp(['Beam calculation took ' num2str(toc) ' seconds']);

%% Generate force/position graphs

tic

% Calculate the force along z
z = [0;0;1]*linspace(-8,8,80)*wavelength_medium;
fz = ott.forcetorque(beam, T, 'position', z);

% Find the equilibrium along the z axis
zeq = ott.find_equilibrium(z(3, :), fz(3, :));
if isempty(zeq)
  warning('No axial equilibrium in range!')
  zeq=0;
end
zeq = zeq(1);

% Calculate force along x-axis (with z = zeq, if found)
r = [1;0;0]*linspace(-4,4,80)*wavelength_medium + [0;0;zeq];
fr = ott.forcetorque(beam, T, 'position', r);

disp(['Force calculation took ' num2str(toc) ' seconds']);

%% Generate the plots

figure(1); plot(z(3, :)/wavelength_medium,fz(3, :));
xlabel('{\it z} [\lambda_m]');
ylabel('{\it Q_z} [n_m P / c]');
aa = axis;
hold on;
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
hold off;

figure(2); plot(r(1, :)/wavelength_medium,fr(1, :));
xlabel('{\it r} [\lambda_m]');
ylabel('{\it Q_r} [n_m P / c]');
aa = axis;
hold on;
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
hold off;
