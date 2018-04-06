% Demonstrates generation of near and far field visualisations of the beams
%
% Shows the radiance at the focal plane, transverse to the focal
% plane and in the far field.
%
% Multiple beam configurations are included, use the beam_type
% variable to vary the beam that is being visualised.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

%% Setup the beam

% Select a beam type to visualise (gaussian, lg, hg, ig, addition, scattered)
beam_type = 'gaussian';

switch beam_type
  case 'gaussian'

    % Create a simple Gaussian beam with circular polarisation
    beam = ott.BscPmGauss('NA', 1.02, 'polarisation', [ 1 1i ]);

  case 'lg'

    % Create a simple LG03 beam with circular polarisation
    beam = ott.BscPmGauss('lg', [ 0 3 ], ...
        'NA', 1.02, 'polarisation', [ 1 1i ]);

  case 'hg'

    % Create a HG23 beam with circular polarisation
    %
    % The beam waist calculation doesn't yet handle HG beams so
    % we estimate the beam waist using a magic formula
    beam_angle = asin(NA/n_medium);
    w0=1.0/pi/tan(abs(theta/180*pi));
    beam = ott.BscPmGauss('hg', [ 2 3 ], ...
        'polarisation', polarisation, 'w0', w0);

  case 'ig'

    % Create a IG beam, not sure what would look good yet
    %
    % The beam waist calculation doesn't yet handle IG beams so
    % we estimate the beam waist using a magic formula
    beam_angle = asin(NA/n_medium);
    w0=1.0/pi/tan(abs(theta/180*pi));
    beam = ott.BscPmGauss('ig', [ 2 3 1 2 ], ...
        'polarisation', polarisation, 'w0', w0);

  case 'addition'

    % Separation between beams
    dx = 1.0;

    % Create a simple Gaussian beam with circular polarisation
    beam = ott.BscPmGauss('NA', 1.02, 'polarisation', [ 1 1i ]);
    beam.Nmax = beam.Nmax + ott.utils.ka2nmax(dx);

    % Displace the beams by +/- dx/2
    beam1 = beam.translateXyz([dx/2, 0, 0]);
    beam2 = beam.translateXyz([-dx/2, 0, 0]);

    % Join the beams to create our beam (adding a phase shift to the 2nd beam)
    beam = beam1 + beam2 * exp(2 * pi * 1i * 1/2);

  case 'scattered'

    % Create a simple Gaussian beam with circular polarisation
    ibeam = ott.BscPmGauss('NA', 1.02, 'polarisation', [ 1 1i ]);

    % Create a spherical particle to scatter the beam
    tmatrix = ott.Tmatrix.simple('sphere', 1.0, 'wavelength0', 1.0, ...
        'n_medium', 1.0, 'n_particle', 1.2);

    % Scatter the beam to create the final beam
    beam = tmatrix * ibeam;

  otherwise
    error('Unsupported beam type specified');
end

%% Generate plot of the intensity at the focal plane

% Create the grid of points we want to view
nx = 40;
ny = 40;
xrange = linspace(-8, 8, nx);
yrange = linspace(-8, 8, ny);
[xx, yy] = meshgrid(xrange, yrange);
xyz = [xx(:) yy(:) zeros(size(xx(:)))].';

% Calculate the E and H fields
[E, H] = beam.emFieldXyz(xyz);

% Calculate the intensity of the E field
Ei=reshape(sqrt(sum(real(E).^2,1)),[nx,ny]);

% Calculate the radiance
I=reshape(sum(abs(E).^2,1),[nx,ny]);

figure(1);
subplot(1, 2, 1);
imagesc(xrange, yrange, Ei);
title('E field intensity (focal plane)');
subplot(1, 2, 2);
imagesc(xrange, yrange, I);
title('radiance (focal plane)');

%% Generate plot of the intensity at the focal plane (same as above)

% Create the grid of points we want to view
nx = 40;
ny = 40;
xrange = linspace(-8, 8, nx);
yrange = linspace(-8, 8, ny);
[xx, yy] = meshgrid(xrange, yrange);
xyz = [xx(:) zeros(size(xx(:))) yy(:)].';

% Calculate the E and H fields
[E, H] = beam.emFieldXyz(xyz);

% Calculate the intensity of the E field
Ei=reshape(sqrt(sum(real(E).^2,1)),[nx,ny]);

% Calculate the radiance
I=reshape(sum(abs(E).^2,1),[nx,ny]);

figure(2);
subplot(1, 2, 1);
imagesc(xrange, yrange, Ei);
title('E field intensity (cross-section)');
subplot(1, 2, 2);
imagesc(xrange, yrange, I);
title('radiance (cross-section)');

%% Generate a figure showing the farfield

%build grid:
nt=80;
[x,y,z]=sphere(nt);

%generate angular points for farfield:
[~,theta,phi]=xyz2rtp(x,y,z);

%find far-field in theta, phi:
[E,H]=beam.farfield(theta(:),phi(:));

% Calculate the intensity of the E field
Ei=reshape(sqrt(sum(real(E).^2,1)),[nt+1,nt+1]);

% Calculate the average phase of the E field
% This should be similar to the pattern we put on a SLM
Ep=reshape(sum(arg(E),1)/2.0),[nt+1,nt+1]);

% Calculate the radiance
I=reshape(sum(abs(E).^2,1),[nt+1,nt+1]);

figure(3);
subplot(1, 2, 1);
surf(x,y,z,Ei,'facecolor','interp','edgecolor','none')
title('E field intensity (farfield)');
subplot(1, 2, 2);
surf(x,y,z,Ep,'facecolor','interp','edgecolor','none')
title('E field phase (farfield)');
subplot(1, 2, 3);
surf(x,y,z,I,'facecolor','interp','edgecolor','none')
title('radiance (farfield)');

