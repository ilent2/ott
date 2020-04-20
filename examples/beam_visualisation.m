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

% Add the toolbox to the path (assuming we are in ott/examples)
addpath('../');

% Close open figures
close all;

% Turn off change warnings for examples
ott.change_warnings('off');

%% Setup the beam

% Select a beam type to visualise (gaussian, lg, hg, ig, addition, scattered)
beam_type = 'lg';

% Specify Nmax explicity for our near-field visualisation
Nmax = ott.utils.ka2nmax(10*2*pi);

sca = 1.0;

switch beam_type
  case 'gaussian'

    % Create a simple Gaussian beam with circular polarisation
    beam = ott.BscPmGauss('NA', 0.9, 'polarisation', [ 1 0 ], ...
        'index_medium', 1, 'wavelength0', 1, 'Nmax', Nmax);

  case 'lg'

    % Create a simple LG03 beam with circular polarisation
    beam = ott.BscPmGauss('lg', [ 0 3 ], ...
        'index_medium', 1, 'wavelength0', 1, ...
        'NA', 0.8, 'polarisation', [ 1 1i ], 'Nmax', Nmax);

  case 'hg'

    % Create a HG23 beam with circular polarisation
    beam = ott.BscPmGauss('hg', [ 2 3 ], ...
        'index_medium', 1, 'wavelength0', 1, ...
        'NA', 0.8, 'polarisation', [ 0 1 ], 'Nmax', Nmax);

    sca = 2.0;

    % TODO: IG example

  case 'addition'

    % Separation between beams
    dx = 1.4;

    % Create a simple Gaussian beam with circular polarisation
    beam = ott.BscPmGauss('NA', 0.8, 'polarisation', [ 1 1i ], ...
        'index_medium', 1, 'wavelength0', 1, 'Nmax', Nmax);
    beam.Nmax = beam.Nmax + ott.utils.ka2nmax(dx);

    % Displace the beams by +/- dx/2
    beam1 = beam.translateXyz([dx/2; 0; 0]);
    beam2 = beam.translateXyz([-dx/2; 0; 0]);

    % Join the beams to create our beam (adding a phase shift to the 2nd beam)
    beam = beam1 + beam2 * exp(2 * pi * 1i * 1/2);

  case 'scattered'

    % Create a simple Gaussian beam with circular polarisation
    ibeam = ott.BscPmGauss('NA', 0.8, 'polarisation', [ 1 1i ], ...
        'index_medium', 1, 'wavelength0', 1, 'Nmax', Nmax);
    ibeam = ibeam.translateXyz([1; 0; 0]);

    % Create a spherical particle to scatter the beam
    tmatrix = ott.Tmatrix.simple('sphere', 0.5, 'wavelength0', 1.0, ...
        'index_medium', 1.0, 'index_particle', 1.5);

    % Scatter the beam to create the final beam
    beam = tmatrix * ibeam;
    
    % We can visualise the scattered or change the field to total
    beam = beam.totalField(ibeam);

  otherwise
    error('Unsupported beam type specified');
end

%% Generate plot of the intensity at the focal plane

% Create the grid of points we want to view
nx = 80;
ny = 80;
xrange = linspace(-2, 2, nx)*sca;
yrange = linspace(-2, 2, ny)*sca;
[xx, yy] = meshgrid(xrange, yrange);
xyz = [xx(:) yy(:) zeros(size(xx(:)))].';

% Calculate the E and H near-fields
% For this we need to use the regular VSWF basis (finite at origin)
% For a scattering particle, we could also choose to visulise the
% scattered field or the total field
beam.basis = 'regular';
[E, H] = beam.emFieldXyz(xyz);

% Calculate the intensity of the E field
Ei=reshape(sqrt(sum(real(E).^2,1)),[nx,ny]);

% Calculate the phase at the focal plane
Ep = reshape(angle(E(3, :)),[nx,ny]);

% Calculate the radiance
I=reshape(sum(abs(E).^2,1),[nx,ny]);

% The fields we calculate above are only valid outside the
% particle.  We could do a second scattering calculation to determine
% the fields inside the particle
if strcmp(beam_type, 'scattered')
  idx = ((xx - 1).^2 + yy.^2) <  0.5.^2;
  I(idx) = NaN;
  Ei(idx) = NaN;
  Ep(idx) = NaN;
end

figure(1);
subplot(1, 3, 1);
imagesc(xrange, yrange, Ei);
axis image;
title('E field intensity');
xlabel('X [\lambda_0]'); ylabel('Y [\lambda_0]');
subplot(1, 3, 2);
imagesc(xrange, yrange, Ep);
axis image;
title('E field phase');
xlabel('X [\lambda_0]'); ylabel('Y [\lambda_0]');
subplot(1, 3, 3);
imagesc(xrange, yrange, I);
axis image;
title('radiance');
xlabel('X [\lambda_0]'); ylabel('Y [\lambda_0]');
set(gcf, 'Name','Focal plane','NumberTitle','off');

%% Generate plot of the intensity along the beam (same as above)

% Create the grid of points we want to view
nx = 40;
ny = 40;
xrange = linspace(-2, 2, nx)*sca;
yrange = linspace(-2, 2, ny)*sca;
[xx, yy] = meshgrid(xrange, yrange);
xyz = [xx(:) zeros(size(xx(:))) yy(:)].';

% Calculate the E and H fields
% For this we need to use the regular VSWF basis (finite at origin)
% For a scattering particle, we could also choose to visulise the
% scattered field or the total field
beam.basis = 'regular';
[E, H] = beam.emFieldXyz(xyz);

% Calculate the intensity of the E field
Ei=reshape(sqrt(sum(real(E).^2,1)),[nx,ny]);

% Calculate the radiance
I=reshape(sum(abs(E).^2,1),[nx,ny]);

% The fields we calculate above are only valid outside the
% particle.  We could do a second scattering calculation to determine
% the fields inside the particle
if strcmp(beam_type, 'scattered')
  idx = ((xx - 1).^2 + yy.^2) <  0.5.^2;
  I(idx) = NaN;
  Ei(idx) = NaN;
  Ep(idx) = NaN;
end

figure(2);
subplot(1, 2, 1);
imagesc(xrange, yrange, Ei);
axis image;
title('E field intensity');
xlabel('X [\lambda_0]'); ylabel('Z [\lambda_0]');
subplot(1, 2, 2);
imagesc(xrange, yrange, I);
axis image;
title('radiance');
xlabel('X [\lambda_0]'); ylabel('Z [\lambda_0]');
set(gcf, 'Name','cross-section','NumberTitle','off');

%% Generate a figure showing the farfield

%build grid:
nt=160;
[x,y,z]=sphere(nt);

%generate angular points for farfield:
[~,theta,phi]=ott.utils.xyz2rtp(x,y,z);

% find far-field in theta, phi
% This requires an outgoing or incoming beam
beam.basis = 'outgoing';
[E,H]=beam.farfield(theta(:),phi(:));

% Calculate the phase of an E field component
Ep = reshape(angle(E(3, :)),[nt+1,nt+1]);

% Add a small translation, so it looks like a LG fork in farfield
% This is the pattern we would put on a SLM
% The translation requires a regular basis, the visualisation requires
% eitehr an outgoing or incoming basis.
beam.basis = 'regular';
slm_beam = beam.translateXyz([5; 0; 0]);
slm_beam.basis = 'outgoing';
[E_slm,~]=slm_beam.farfield(theta(:),phi(:));
Ep_slm = reshape(angle(E_slm(3, :)),[nt+1,nt+1]);

% Calculate the radiance
I=reshape(sum(abs(E).^2,1),[nt+1,nt+1]);

figure(3);
subplot(1, 3, 1);
surf(x,y,z,I,'facecolor','interp','edgecolor','none')
title('radiance');
xlabel('X [r]'); ylabel('Y [r]');
view(0, 90);
subplot(1, 3, 2);
surf(x,y,z,Ep,'facecolor','interp','edgecolor','none')
title('phase');
xlabel('X [r]'); ylabel('Y [r]');
view(0, 90);
subplot(1, 3, 3);
surf(x,y,z,Ep_slm,'facecolor','interp','edgecolor','none')
title('phase + offset');
xlabel('X [r]'); ylabel('Y [r]');
view(0, 90);
set(gcf, 'Name','farfield','NumberTitle','off');

% Make sure the user knows the beam is more interesting for scattered
if strcmp(beam_type, 'scattered')
  for ii = 1:3, subplot(1, 3, ii), view(3), end
end
