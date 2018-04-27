% Example of calculation of restoring forces on a layered sphere.
%
% Approximate of Figure 3a/b in Hu et al., Antireflection coating for 
% improved optical trapping, Journal of Applied Physics 103, 093119 (2008)
%
% How long should this take?
% Using the the parameters required to fairly faithfully reproduce figure 3
% took 800 seconds on a Core2 Duo 6600 w/ 6GB RAM. The version here uses a
% much coarser grid and takes 200 seconds.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add the toolbox to the path (assuming we are in ott/examples)
addpath('../');

% Close open figures
close all;

% Make warnings less obtrusive
ott.warning('once');
ott.change_warnings('off');

% Specify refractive indices
n_medium = 1.34;
n_particle = 1.8;

% Calculate the shell refractive index
n_shell = sqrt(n_medium*n_particle);

% Set the wavelength in freespace [m]
wavelength = 1064.0e-9;
k = 2*pi/wavelength;

% Radius of particle to use [freespace units]
radius = linspace(.25,2.5,50)/n_medium;

%% Generate beam

tic

% Gaussian beam with circular polarisation
beam = ott.BscPmGauss('NA', 1.02, 'polarisation', [ 1 1i ], ...
    'power', 1.0, 'index_medium', n_medium, 'wavelength0', wavelength);

disp(['Beam calculation took ' num2str(toc) ' seconds']);

%% Calculate trap depths for a bunch of different spheres

tic

% Allocate spcae for the arrays
z = zeros(1, 45);
r = zeros(1, 25);
fz = zeros(6, length(z));
fr = zeros(6, length(r));

axialrestoringforce = zeros(length(radius), 1);
axialrestoringforce1 = axialrestoringforce;
transverserestoringforce = zeros(length(radius), 1);
transverserestoringforce1 = transverserestoringforce;

% For each particle radius, calculate the axial and radial min trap depth
for ii=1:length(radius)

    disp(['Calculating radius ' num2str(ii) '/' ...
        num2str(length(radius)) ' (' num2str(radius(ii)) ')']);

    % Calculate grid [m]
    z = linspace(-.5,max([1,radius(ii)])*n_particle,45)*wavelength;
    r = linspace(0,max([1,radius(ii)])*1.33,25)*wavelength;

    % Calculate T-matrix for spherical particle with shell
    T = ott.TmatrixMie.simple('sphere', ...
        [radius(ii), radius(ii)+1/n_shell/4]*wavelength, ...
        'k_medium', k*n_medium, ...
        'k_particle', k*[n_particle, n_shell]);

    % Calculate T-matrix for spherical particle without shell
    T1 = ott.Tmatrix.simple('sphere', radius(ii)*wavelength, ...
        'k_medium', k*n_medium, 'k_particle', k*n_particle);

    % Calculate the force along the axial direction
    % The force arrays contain the x, y, z force for both particles
    fz = ott.forcetorque(beam, [T, T1], 'position', [0;0;1] * z);

    % Find the minimum axial restoring force
    axialrestoringforce(ii) = min(fz(3, :));
    axialrestoringforce1(ii) = min(fz(6, :));

    % Find z-axis equilibrium
    zeq = ott.find_equilibrium(z, fz(3, :));
    if isempty(zeq)
      ott.warning(['No zeq in range (radius = ' num2str(radius(ii)) ')']);
      zeq=0;
    else
      zeq = zeq(1);
    end

    % Calculate the radial force at z-axial equilibrium
    fr = ott.forcetorque(beam, [T, T1], ...
        'position', [1;0;0] * r + [0;0;zeq]);

    % Find the minimum transverse restoring force
    transverserestoringforce(ii) = min(fr(1, :));
    transverserestoringforce1(ii) = min(fr(4, :));
end

disp(['Simulation took ' num2str(toc()) ' seconds']);

%% Generate figures of the trap depth vs core radius

figure(1)
plot(radius*n_medium,-axialrestoringforce);
hold on
plot(radius*n_medium,-axialrestoringforce1,':');
hold off
legend('shell', 'solid', 'Location', 'NorthWest');
ylabel('axial force efficiency \it Q_z')
xlabel('core radius [{\it{\lambda_{medium}}}]')
xlim(([radius(1),radius(end)]*n_medium))
ylim([0,max(ylim)])
grid on

figure(2)
plot(radius*n_medium,-transverserestoringforce);
hold on
plot(radius*n_medium,-transverserestoringforce1,':');
hold off
legend('shell', 'solid', 'Location', 'NorthWest');
ylabel('transverse force efficiency \it Q_r')
xlabel('core radius [{\it{\lambda_{medium}}}]')
xlim(([radius(1),radius(end)]*n_medium))
ylim([0,max(ylim)])
grid on
