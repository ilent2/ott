% Generate a T-matrix for a vaterite particle using DDA
%
% This example shows how DDA can be used to simulate a vaterite
% particle either as a homogeneous spherical particle or with a
% sheaf-of-wheat structure.
%
% The example calculates the T-matrix, calculates the torque when
% the particle is placed in a Gaussian beam and shows runtime information.
%
% This script takes a couple of minutes on a i7 with 16 GB ram.
% If this is too slow on your system, try reducing the number
% of vaterites and the vaterite size.  If you run out of memory,
% try the slower low memory DDA (see comment bellow).
%
% This simulation uses a rather large spacing (lambda/10) and is likely
% very inacurate, a spacing of lambda/20 or smaller would be preferable.
%
% Copyright 2020 Isaac Lenton.
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add OTT to the path (only required if OTT isn't already on the path)
% This assumes this script is in the ott/examples directory, if you
% moved the script you may need to modify this command
addpath('../');

% Make warnings less obtrusive
ott.warning('once');
ott.change_warnings('off');

%% Define properties for the particle and beam

% Refractive indices of water and vaterite
% Vaterite is positive uniaxial (n_e > n_o)
index_medium = 1.33;
index_o = 1.55;
index_e = 1.65;

% Vacuum wavelength [m]
wavelength0 = 1064e-9;

% Medium and particle wavelengths [m]
wavelength_medium = wavelength0 ./ index_medium;
wavelength_particle = wavelength0 ./ [index_o, index_e];

% Describe the particle properties
spacing = min(wavelength_particle) / 10;

%% Visualise the sheaf-of-wheat structur
%
% In this example we model the alignment of the vaterite crystal with
% a sheaf-of-wheat structure.  The definition is at the end of this file.
%
% The following code generates a image of the structure.

sc = 1e-6;
tol = 1e-16;   % Shift focus slightly away from zero to supress warning

x = linspace(-1, 1, 100)*sc;
y = linspace(-1, 1, 100)*sc;
[xx, yy] = meshgrid(x, y);
r = [0*xx(:), xx(:), yy(:)].';

v = sheafOfWheat(r+tol, 0.5*sc);

small = 1.0e-4*sc;
figure();
streamline(xx, yy, ...
  reshape(v(2, :), size(xx)), reshape(v(3, :), size(xx)), ...
  [linspace(-1, 1, 100), linspace(-1, 1, 100)]*sc, [zeros(1, 100)+small, zeros(1, 100)-small]);
axis image
title('Sheaf-of-wheat visualisation');
xlabel('X Position [m]');
ylabel('Z Position [m]');

%% Generate the T-matrix
%
% To calculates the T-matrix we use DDA with mirror symmetry and
% rotational symmetry optimisations.

% DDA has a low memory implementation but it is slower
% Enable this if you run out of ram for the larger particles
low_memory = false;

% The method supports using '\' and 'gmres'.  Set use iterative to
% true to use 'gmres' for solving the DDA system.  (Might be faster).
use_iterative = false;

radius = logspace(-8, -6.2147, 20);
% radius = 1.0e-6;   % Try with low_memory
% radius = logspace(-8, -6, 20);
times = zeros(size(radius));

% Setup variable for T-matrix (ensure it is free)
Ts = ott.TmatrixDda.empty(1, 0);

for ii = 1:length(times)

    % Pre-calculate voxel locations
    shape = ott.shapes.Sphere(radius(ii));
    voxels = shape.voxels(spacing, 'even_range', true);
    
    if numel(voxels) == 0
      times(ii) = nan;
      continue;
    end
    
    % Calculate polarizability unit
    upol = ott.utils.polarizability.LDR(spacing ./ wavelength_medium, ...
        [index_o; index_o; index_e] ./ index_medium);
    
    % Calculate sheaf-of-wheta orientations
    dirs = sheafOfWheat(voxels, 0.5*radius(ii));
    
    % Rotate polarizability and index_relative units
    polarizabilities = ott.utils.rotate_3x3tensor(diag(upol), 'dir', dirs);
      
    tic();

    Tmatrix = ott.TmatrixDda(voxels, ...
        'polarizability', polarizabilities, ...
        'index_medium', index_medium, ...
        'spacing', spacing, ...
        'z_rotational_symmetry', 0, ...
        'z_mirror_symmetry', true, ...
        'wavelength0', wavelength0, ...
        'low_memory', low_memory, ...
        'use_iterative', use_iterative);

    times(ii) = toc();
    
    Ts(end+1) = Tmatrix;
end

%% Generate plot of runtime

figure();
loglog(radius, times, '*--');
xlabel('Radius [m]');
ylabel('Time [s]');
title('DDA Run-time');

%% Generate plot of the torque
%
% This is probably not very acurate and not very large or tiny particles.
%
% If we do larger particles and improve the accuracy we should be able
% to approximate figure 3 from
%
%   Highly birefringent vaterite microspheres: production,
%   characterization and applications for optical micromanipulation.
%   S. Parkin, et al., Optics Express Vol. 17, Issue 24, pp. 21944-21955 (2009)
%   https://doi.org/10.1364/OE.17.021944
%
% Tips to try to reproduce this figure: specify the ms to use with
% TmatrixDda to match the ms in the beam.  Increase resolution.

beam = ott.BscPmGauss('NA', 1.1, 'index_medium', index_medium, ...
  'power', 1.0, 'wavelength0', wavelength0, 'polarisation', [1, -1i]);

[~, tz] = ott.forcetorque(beam, Ts, ...
    'position', [0;0;0], 'rotation', eye(3));

figure();
plot(radius(end-length(Ts)+1:end), tz(3:3:end, :));
xlabel('Radius [m]');
ylabel('Torque [a.u.]');

%% Sheaf-of-wheat definition

function v = sheafOfWheat(voxels, foci_distance)
% Generate the vectors for the sheaf-of-wheat structure
%
% - voxels (3xN numeric) -- voxel locations
% - foci_distance (numeric) -- distance from centre to focii
%
% Based on code from Vincent Loke (circa 2009)

if any(voxels(:) == 0)
  warning('Voxels should not be placed on planes through 0');
end

xneg = voxels(1, :) < 0;
yneg = voxels(2, :) < 0;
zneg = voxels(3, :) < 0;

voxels = abs(voxels);

rho = sqrt(voxels(1, :).^2 + voxels(2, :).^2);
z = voxels(3, :);

gam = asin((2^(1/2)*((foci_distance^4 - 2*foci_distance^2*rho.^2 ...
  + 2*foci_distance^2*z.^2 + rho.^4 + 2*rho.^2.*z.^2 + z.^4).^(1/2) ...
  + foci_distance^2 - rho.^2 - z.^2).^(1/2))/(2*foci_distance));

dzdrho = foci_distance/2*sin(gam).* ...
  ((rho./(foci_distance*cos(gam))).^2-1).^(-.5).* ...
  (2*rho)./(foci_distance*cos(gam)).^2;

theta = pi/2 - atan(dzdrho); % theta is relative to the z-axis
rho = sin(theta);
phi = atan(voxels(2,:)./voxels(1,:));

v = [rho.*cos(phi); rho.*sin(phi); cos(theta)];

v(1, xneg) = -v(1, xneg);
v(2, yneg) = -v(2, yneg);
v(3, zneg) = -v(3, zneg);

end
