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

% Make warnings less obtrusive
ott_warning('once');
change_warnings('off');

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

% Gaussian beam with circular polarisation
beam = ott.BscPmGauss('NA', 1.02, 'polarisation' [ 1 i ]);
beam = beam / beam.power();

% Allocate spcae for the arrays
z = zeros(1, 45);
r = zeros(1, 25);
fz = zeros(3, length(z));
fr = zeros(3, length(r));
fz1 = fz;
fr1 = fr;

axialrestoringforce = zeros(length(radius), 1);
axialrestoringforce1 = axialrestoringforce;
transverserestoringforce = zeros(length(radius), 1);
transverserestoringforce1 = transverserestoringforce;

% For each particle radius, calculate the axial and radial min trap depth
for ii=1:length(radius)

    % Calculate grid [in medium units]
    z = linspace(-.5,max([1,radius(ii)])*n_particle,45)*n_medium;
    r = linspace(0,max([1,radius(ii)])*1.33,25)*n_medium;

    % Calculate T-matrix for spherical particle with shepp
    T = ott.Tmatrix.simple('sphere', [radius(ii), radius(ii)+1/n_shell/4], ...
        'k_medium', k*n_medium, ...
        'k_particle', k*[n_particle, n_shell]);

    % Calculate T-matrix for spherical particle without shell
    T1 = ott.Tmatrix.simple('sphere', radius(ii), ...
        'k_medium', k*n_medium, 'k_particle', k*n_particle);
    
    % Calculate the force along the axial direction
    for nz = 1:length(z)
      tbeam = beam.translateXyz(0, 0, z(nz));
      sbeam = T * tbeam;
      sbeam1 = T1 * tbeam;
      
      fz(:, nr) = ott.forcetorque(tbeam, sbeam);
      fz1(:, nr) = ott.forcetorque(tbeam, sbeam1);
    end

    % Find the minimum axial restoring force
    axialrestoringforce(ii) = min(fz(3, :));
    axialrestoringforce1(ii) = min(fz1(3, :));
    
    % Find z-axis equilibrium
    zeq = ott.find_equilibrium(z, fz(3, :));
    if isempty(zeq)
      warning(['No zeq in range (radius = ' num2str(radius(ii)) ')']);
      zeq=0;
    end

    % Calculate the radial force at z-axial equilibrium
    for nr = 1:length(r)
      tbeam = beam.translateXyz(r(nr), 0, zeq);
      sbeam = T * tbeam;
      sbeam1 = T1 * tbeam;
      
      fr(:, nr) = ott.forcetorque(tbeam, sbeam);
      fr1(:, nr) = ott.forcetorque(tbeam, sbeam1);
    end
    
    % Find the minimum transverse restoring force
    transverserestoringforce(ii) = min(fr(1, :));
    transverserestoringforce1(ii) = min(fr1(1, :));
end

%% Generate figures of the trap depth vs core radius

figure(1)
plot(radius*n_medium,-axialrestoringforce);
hold on
plot(radius*n_medium,-axialrestoringforce1,':');
hold off
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
ylabel('transverse force efficiency \it Q_r')
xlabel('core radius [{\it{\lambda_{medium}}}]')
xlim(([radius(1),radius(end)]*n_medium))
ylim([0,max(ylim)])
grid on
