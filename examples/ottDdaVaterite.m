% ottDdaVaterite.m    -- Simulate a birefringent particle using DDA
%
% This example shows how to simulate a birefringent spherical particle
% using the discrete dipole approximation.  For a simpler example showing
% how to simulate a homogeneous non-spherical particle,
% see `ottNonSpherical.m`
%
% The example calculates the T-matrix, the torque when the particle is
% placed in a Gaussian beam, and the runtime information.
%
% This script takes a couple of minutes on a i7 with 16 GB ram.
% If this is too slow on your system, try reducing the number
% of particles and the vaterite size.  If you run out of memory,
% try the slower low memory DDA (see comment bellow).
%
% This simulation uses a rather large spacing (lambda/10) and is likely
% very inacurate, a spacing of lambda/20 or smaller would be preferable.
%
% Copyright 2020 Isaac Lenton.
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../');

% On R2018a this script seems to crash matlab when variables/classes are not
% cleared before running, you may want to comment this line...
clearvars;

%% Define properties for particles/beams

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

%% Visualise the sheaf-of-wheat structure
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

radius = logspace(-7.0, -6.2147, 20);
times = zeros(size(radius));

% Clear/allocate memory for T-matrices
Ts = ott.tmatrix.Tmatrix();
Ts = repmat(Ts, 1, numel(radius));

particle = ott.particle.Fixed();
particle = repmat(particle, 1, numel(radius));

for ii = 1:length(times)
  
    disp(['Particle number : ' num2str(ii)]);

    % Pre-calculate voxel locations
    shape = ott.shape.Sphere(radius(ii));
    voxels = shape.voxels(spacing, 'even_range', true);

    if numel(voxels) == 0
      times(ii) = nan;
      continue;
    end

    % Calculate polarizability unit
    upol = ott.tmatrix.dda.polarizability.LDR(spacing ./ wavelength0, ...
        [index_o; index_o; index_e] ./ index_medium);

    % Calculate sheaf-of-wheta orientations
    dirs = sheafOfWheat(voxels, 0.5*radius(ii));

    % Rotate polarizability and index_relative units
    % We need index_relative to ensure correct scaling of the TE modes
    index_relative = ott.utils.rotate_3x3tensor(...
        diag([index_o, index_o, index_e] ./ index_medium), 'dir', dirs);
    polarizabilities = ott.utils.rotate_3x3tensor(...
        diag(upol), 'dir', dirs);

    tic();

    % Could also use 'ott.tmatrix.dda.Dda' which may run better
    % when memory is limited.
    dda = ott.tmatrix.dda.DdaHighMem(voxels./wavelength_medium, ...
        polarizabilities, 'xySymmetry', true, 'zRotSymmetry', 0);

    Ts(ii) = ott.tmatrix.Dda(dda);
    particle(ii) = ott.particle.Fixed('tmatrix', Ts(ii));

    times(ii) = toc();
end

%% Generate plot of runtime

figure();
loglog(radius, times, '*--');
xlabel('Radius [m]');
ylabel('Time [s]');
title('DDA Run-time');

%% Generate plot of the torque
%
% This is probably not very accurate for not very large or tiny particles.
%
% If we do larger particles and improve the accuracy we should be able
% to approximate figure 3 from
%
%   Highly birefringent vaterite microspheres: production,
%   characterization and applications for optical micromanipulation.
%   S. Parkin, et al., Optics Express Vol. 17, Issue 24, pp. 21944-21955 (2009)
%   https://doi.org/10.1364/OE.17.021944
%
% Tips to try to reproduce this figure: specify the ci-orders to use with
% TmatrixDda to match the ci-orders in the beam.  Increase resolution.

beam = ott.beam.Gaussian.FromNa(1.1, 'index_medium', index_medium, ...
  'power', 1.0, 'wavelength0', wavelength0, ...
  'polfield', [1, -1i], 'polbasis', 'cartesian');

tz = -beam.torque(particle, 'position', [0;0;0], 'rotation', eye(3));

figure();
plot(radius(end-length(Ts)+1:end), tz(3, :));
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
