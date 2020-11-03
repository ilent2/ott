% ottLandscape.m -- Calculate trap stiffness for a range of particles
%
% This example shows how a trapping landscape plot showing trap stiffness
% for a range of different particle refractive indices and sizes can be
% generated.  This is based on Figure 3 from
%
%   Nieminen et al., 2007, Journal of Optics
%   https://doi.org/10.1088/1464-4258/9/8/S12
%
% Further examples can be found in
%
%   Stilgoe et al., 2008, Optics Express
%   https://doi.org/10.1364/OE.16.015039
%
% and discussion in
%
%   Bo Sun and David G. Grier, 2009, Optics Express
%   https://doi.org/10.1364/OE.17.002658
%
% As noted in the later, in a previous version of the toolbox and with
% certain versions of Matlab, small particles would be misrepresented.
% Regardless, care should be taken when interpreting these diagrams.
%
% This file is an example from the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../');

% Resolution for image
NR = 100;
NI = 50;

% Radius and RI values to match paper
radius = linspace(1e-2, 3.25, NR)*1e-6/2;   % Sphere radius [m]
index = linspace(1.34, 2.66, NI);           % Refractive index

% Setup the indecent beam
NA = 1.3;
beam = ott.beam.Gaussian.FromNa(NA, ...
    'index_medium', 1.33, 'wavelength0', 1064e-9, ...
    'polbasis', 'cartesian', 'polfield', [1, 1i], ...
    'truncation_angle', pi/2);

% Allocate memory for traps
traps = ott.tools.FindTraps1d();
traps(numel(index), numel(radius)) = traps;

%% Run calculation

tic

for ii = 1:numel(radius)

  % Display progress
  if mod(ii, 10) == 0
    disp(['Progress: ' num2str(ii) ' / ' num2str(numel(radius))]);
  end

  % Calculate beam translation distances
  if radius(ii) > 0.5*beam.wavelength
    zScale = (beam.index_medium./NA)^2*(0.5+radius(ii)./beam.wavelength ...
        * max(index ./ beam.index_medium));
  else
    zScale = (beam.index_medium./NA)^2*(1+3/4*radius(ii)./beam.wavelength);
  end
  z = linspace(-1, 1, 100)*zScale;

  % Translate beam
  % Use the maximum particle size for translation N-max (we don't care
  % about multipole modes above Nmax since they aren't involved in
  % scattering).
  tbeams = beam.getData([]).translateXyz([0;0;1].*z, ...
      'Nmax', ott.utils.ka2nmax(beam.wavenumber*max(radius)));

  for jj = 1:numel(index)

    particle = ott.tmatrix.Mie(radius(ii)./beam.wavelength, ...
        index(jj)./beam.index_medium);

    % Calculate force and find traps
    fz = tbeams.force(particle);
    our_traps = ott.tools.FindTraps1d.FromArray(z, fz(3, :));
    if ~isempty(our_traps)
      traps(jj, ii) = our_traps(1);
    end
  end
end

toc

%% Generate a plot

figure();
maxforce = [traps.max_force];
imagesc(radius, index, maxforce);
xlabel('Radius [m]');
ylabel('Refractive Index');
cb = colorbar();
yaxis(cb, 'Q_{max}');

