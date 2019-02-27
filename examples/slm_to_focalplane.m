% Calculates how a SLM pattern would look at the focal plane of a microscope
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add the toolbox path
addpath('../');

% Close open figures
close all;

% Turn off change warnings for examples
ott.change_warnings('off');

% Numerical aperture of microscope
NA = 1.02;

% Refractive index of medium
index_medium = 1.33;

% Polarisation of incident beam
polarisation = [ 0 1 ];

% Range [in number of wavelengths in medium] we want to visualise fields
% near focal plane
range = 3.0;

% Describes how well matched the beam is to the SLM, to underfill
% changes this value to less than 1.
fillfactor = 1.0;

%% Setup the incident beam and SLM pattern

% Define a grid representing pixels on the SLM
x = linspace(-1, 1, 512);
y = linspace(-1, 1, 512);
x_idx = 1:length(x);
y_idx = 1:length(y);
[xx, yy] = meshgrid(x, y);
[xx_idx, yy_idx] = meshgrid(x_idx, y_idx);

% Calculate the incident beam amplitude (Gaussian beam, LG00)
incident_mode = ott.utils.lgmode(0,0, ...
    sqrt(xx.^2+yy.^2)/fillfactor, atan2(yy,xx));

% Generate the phase pattern
slm_phase = zeros(size(incident_mode));

% Set the entire pattern to checkerboard (zero amplitude)
slm_phase(logical(mod(xx_idx + yy_idx, 2))) = pi;

% Unmask the region we want our pattern on
xxyy = ott.utils.rotz(45) * [xx(:) yy(:) zeros(size(xx(:)))].';
slm_phase(sinc(sum([1 0 0] * xxyy, 1)*4) > 0) = pi/2;

%% Show the incident mode and the SLM pattern

figure(1);
imagesc(incident_mode);
title('Incident mode amplitude');
xlabel('X [pixels]'); ylabel('Y [pixels]');
axis('image');

figure(2);
imagesc(slm_phase);
title('SLM phase pattern');
xlabel('X [pixels]'); ylabel('Y [pixels]');
axis('image');

%% Calculate the beam shape coefficients for the paraxial beam

% For the visualisation we need to specify the range where we
% want to plot the fields, so we manually calculate Nmax
Nmax = ott.utils.ka2nmax(2*pi*index_medium*range/2.0);

% Combined the phase and incident beam and make a beam object
combined_beam = exp(1i*slm_phase) .* incident_mode;
beam = ott.BscPmParaxial(-NA, combined_beam, 'Nmax', Nmax, ...
    'polarisation', polarisation, 'index_medium', index_medium);

%% Visualise the near fields (transverse to beam axis)

% Create a grid of point we want to view
nx = 80;
ny = 80;
xrange = linspace(-0.5, 0.5, nx)*range/index_medium;
yrange = linspace(-0.5, 0.5, ny)*range/index_medium;
[xx, yy] = meshgrid(xrange, yrange);
xyz = [xx(:) yy(:) zeros(size(xx(:)))].';

% Calculate the electric field
E = beam.emFieldXyz(xyz);

% Calculate the intensity of the E field
Ei=reshape(sqrt(sum(real(E).^2,1)),[nx,ny]);

% Calculate the phase at the focal plane
Ep = reshape(angle(E(3, :)),[nx,ny]);

% Calculate the radiance
I=reshape(sum(abs(E).^2,1),[nx,ny]);

figure(3);
subplot(1, 3, 1);
imagesc(xrange, yrange, Ei);
axis('image'); xlabel('X [\lambda_m]'); ylabel('Y [\lambda_m]');
title('E field intensity');
subplot(1, 3, 2);
imagesc(xrange, yrange, Ep);
axis('image'); xlabel('X [\lambda_m]'); ylabel('Y [\lambda_m]');
title('E field phase');
subplot(1, 3, 3);
imagesc(xrange, yrange, I);
axis('image'); xlabel('X [\lambda_m]'); ylabel('Y [\lambda_m]');
title('radiance');
set(gcf, 'Name','At focal plane','NumberTitle','off');

%% Visualise the near fields (parallel to beam axis)

% Create a grid of point we want to view
nx = 80;
ny = 80;
xrange = linspace(-0.5, 0.5, nx)*range/index_medium;
yrange = linspace(-0.5, 0.5, ny)*range/index_medium;
[xx, yy] = meshgrid(xrange, yrange);
xyz = [xx(:) zeros(size(xx(:))) yy(:)].';

% Calculate the electric field
E = beam.emFieldXyz(xyz);

% Calculate the intensity of the E field
Ei=reshape(sqrt(sum(real(E).^2,1)),[nx,ny]);

% Calculate the phase at the focal plane
Ep = reshape(angle(E(3, :)),[nx,ny]);

% Calculate the radiance
I=reshape(sum(abs(E).^2,1),[nx,ny]);

figure(4);
subplot(1, 3, 1);
imagesc(xrange, yrange, Ei);
axis('image'); xlabel('X [\lambda_m]'); ylabel('Z [\lambda_m]');
title('E field intensity');
subplot(1, 3, 2);
imagesc(xrange, yrange, Ep);
axis('image'); xlabel('X [\lambda_m]'); ylabel('Z [\lambda_m]');
title('E field phase');
subplot(1, 3, 3);
imagesc(xrange, yrange, I);
axis('image'); xlabel('X [\lambda_m]'); ylabel('Z [\lambda_m]');
title('radiance');
set(gcf, 'Name','Along beam axis','NumberTitle','off');
