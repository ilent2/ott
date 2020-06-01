% Example demonstrating the different beams in the toolbox
%
% This file generates images demonstrating the beams included in the
% toolbox, for a full list of these beams and example output of this
% script, see the Beams section of the documentation.
%
% Most beams can be created using the +abstract beam classes.  These
% classes invoke the implementation classes for field visualisations.
% The implementation classes use various approximations and can produce
% different results.  Additional control over the generated beam can
% be achieved by either explicitly casting the abstract beams or
% directly invoking the beam implementation classes.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% Add the toolbox to the path
addpath('../');

%% Gaussian beams

figure();

% Create an abstract Gaussian beam
paraxial_waist = 0.5;
abstract_beam = ott.beam.abstract.Gaussian(paraxial_waist, ...
    'polarisation', [1, 0]);

% The default abstract visualisation method is the default vswf.Gaussian
% There are two kinds of VSWF representations of a Gaussian beam,
% they differ by the paraxial-to-farfield mapping used for point-matching.
% The default is 'sintheta' mapping, which produces results similar to
% many high-NA microscope objectives
subplot(1, 4, 1);
abstract_beam.visualise('range', [1, 1]);
title('Default: VSWF (mapping:sin)');

% The other mapping is 'tan', which should produce similar results
% for lower focussing (larger paraxial beam waist).
bsc_tantheta = ott.beam.vswf.Bsc(abstract_beam, 'mapping', 'tan');
subplot(1, 4, 2);
bsc_tantheta.visualise('range', [1, 1]);
title('VSWF (mapping:tan)');

% Paraxial beam approximation produces slightly different fields
subplot(1, 4, 3);
paraxial = ott.beam.paraxial.Paraxial(abstract_beam);
paraxial.visualise('range', [1, 1]);
title('Paraxial');

% We can also use a 5th order Davis approximation for the beam
davis = ott.beam.GaussianDavis5(abstract_beam);
subplot(1, 4, 4);
davis.visualise('range', [1, 1]);
title('5-th order Davis');

%% Laguerre-Gaussian beams

figure();

% LG beam with azimuthal mode (lmode = -5) and radial mode (rmode = 3)
paraxial_waist = 2;
abstract_beam = ott.beam.abstract.LaguerreGaussian(paraxial_waist, ...
    'lmode', -5, 'pmode', 3, 'polarisation', [1, 1i]);

% Default visualisation method is VSWF
subplot(1, 2, 1);
abstract_beam.visualise('range', 2*[1, 1]);
title('Default: VSWF');

% Can also construct a paraxial representation
subplot(1, 2, 2);
ott.beam.paraxial.Paraxial(abstract_beam).visualise('range', 2*[1, 1]);
title('Paraxial');

%% Hermite-Gaussian beams

figure();

% HG beam with azimuthal mode (mmode = -5) and radial mode (nmode = 3)
paraxial_waist = 3;
abstract_beam = ott.beam.abstract.HermiteGaussian(paraxial_waist, ...
    'mmode', 5, 'nmode', 3, 'polarisation', [1, 0]);
  
range = [1,1] * 2.5*paraxial_waist;

% Default visualisation method is VSWF
subplot(1, 2, 1);
abstract_beam.visualise('range', range);
title('Default: VSWF');

% Can also construct a paraxial representation
subplot(1, 2, 2);
ott.beam.paraxial.Paraxial(abstract_beam).visualise('range', range);
title('Paraxial');

%% Ince-Gaussian beam

figure();

% Construct an Ince-Gaussian beam
paraxial_waist = 2;
abstract_beam = ott.beam.abstract.InceGaussian(paraxial_waist, ...
    'lmode', 3, 'porder', 5, 'parity', 'even', 'polarisation', [1, 1i]);

% Default visualisation method is VSWF
abstract_beam.visualise('range', [3, 3]);

%% Plane wave beams

figure();

% A single plane wave can be constructed using the abstract class
% The visualisation method casts to the ott.beam.PlaneWave class
subplot(1, 4, 1);
abstract_beam = ott.beam.abstract.PlaneWave();
abstract_beam.visualise('range', [2, 2], 'axis', 'y', 'field', 'Re(Ex)');
title('Single');

% For arrays it is better to invoke ott.beam.PlaneWave directly
% The plane wave class has two constructors, one using directionSet
% and the other (used here) with a direction/polarisation vector.
origin = zeros(3, 5);
directions = randn(3, 5);
field = [ones(1, 5); zeros(1, 5)];
polarisation = cross(directions, randn(3, 5));
plane_waves = ott.beam.PlaneWave.FromDirection('origin', origin, ...
  'direction', directions, 'polarisation', polarisation, ...
  'field', field);

% Arrays of plane waves can be coherent or incoherent
subplot(1, 4, 2);
plane_waves.array_type = 'coherent';
plane_waves.visualise();
title('Coherent array');

% Incoherent plane waves doens't look very interesting
subplot(1, 4, 3);
plane_waves.array_type = 'incoherent';
plane_waves.visualise();
title('Incoherent array');

% Some beam approximations are only valid in a certain range, this is
% particularly true for PlaneWaves with VSWF representation, as demonstrated
% by the following visualisation
subplot(1, 4, 4);
beam = ott.beam.vswf.PlaneBasis(abstract_beam, 'Nmax', 12);
beam.position = [0;0;2.25];
beam = beam.applyTransformation();
beam.visualise('range', [2, 2], 'axis', 'y', 'field', 'Re(Ex)');
title('VSWF with Nmax=12');

%% Geometric ray beams
% Geometric rays are similar to Plane waves except their visualisation
% method draws rays instead of fields and their e/h-field methods default
% to a plane wave with a transverse intensity fall-off term.

figure();

% As with plane waves, the abstract class can be used to construct a
% single ray or for many rays it is better to call the Ray class directly
numrays = 5;
% direction = randn(3, numrays);
direction = [1;1;0].*randn(3, numrays);
origin = [1;1;0].*randn(3, numrays);
polarisation = [0;0;1].*ones(size(direction));
rays = ott.beam.Ray.FromDirection('direction', direction, ...
    'origin', origin, 'polarisation', polarisation);

% Default visualisation method draws rays
% Coherent/incoherent property has no effect for this visualisation
subplot(1, 2, 1);
rays.visualise();
title('Rays');

% A field can be specified, this effectively draws the rays as
% plane-waves with a transverse intensity fall-off
% Coherent/incoherent property changes the visualisation
subplot(1, 2, 2);
rays.visualise('field', 'Re(Ez)', 'range', 10*[1,1]);
caxis([-1, 1]);
title('Rays as Plane Waves');

%% Bessel beam

figure();

theta = pi/4;
abstract_beam = ott.beam.abstract.Bessel(theta);

abstract_beam.visualise();

%% Annular beam

figure();

theta = [pi/8, pi/4];
abstract_beam = ott.beam.abstract.Annular(theta);

abstract_beam.visualise();

%% Dipole beam

% Describes the radiation field of a single dipole
% beam.abstract.Dipole creates a single dipole, or many dipoles can be
% created with beam.Dipole
polarizations = randn(3*5, 1);
positions = randn(3, 5);
dipoles = ott.beam.Dipole('location', positions, ...
    'polarization', polarizations);

% The default visualisation sets the color range so that dipoles
% in the visualisation plane don't over-saturate the image.  This causes
% dipoles to appear as patches of constant intensity.
figure();
dipoles.visualise();

%% Top-Hat beams
% The toolbox includes two abstract top-hat beams: a collimated top-hat beam
% and a focussed top-hat beam.  There are no implementations for these
% beams, isntead the fields are calculated by casting to other beams
% such as vswf.FarfieldPm or Ray.  Be careful using these beams, the accuracy
% depends on the accuracy of the cast/approximation.

figure();

beam = ott.beam.abstract.TopHat('radius', 1.0);
subplot(1, 2, 1);
beam.visualise('range', 2*[1,1]);
title('Collimated Top-hat');

beam = ott.beam.abstract.FocussedTopHat('angle', pi/4);
subplot(1, 2, 2);
beam.visualise('axis', 'y');
title('Focussed Top-hat');

%% Masked beams
% The abstract Masked classes can be used to construct different combinations
% of beams with masked profiles.  The following example shows how to
% represent a Annular mask of a Gaussian beam.

figure();

% Paraxial beam to be masked
gaussian = ott.beam.abstract.Gaussian('waist', 1.0);
subplot(1, 3, 1);
gaussian.visualiseFarfield('direction', 'neg', 'basis', 'incoming');
title('Gaussian Far-field');
caxis([0, 2]);
colorbar();

% Annular mask (and visualisation using sin-theta mapping
% Mask function taks a [theta; phi] array as argument.
mask = @(tp) tp(1, :) > 3*pi/4 & tp(1, :) < 7*pi/8;
subplot(1, 3, 2);
[X, Y] = meshgrid(linspace(-1, 1), linspace(-1, 1));
T = pi - asin(sqrt(X.^2 + Y.^2));
P = atan2(Y, X);
h = pcolor(X, Y, reshape(double(mask([T(:), P(:)].')), size(X)));
axis image;
h.EdgeColor = 'none';
colorbar();
title('Mask');

% Combined beam + mask
beam = ott.beam.abstract.MaskedFarfield('mask', mask, ...
    'masked_beam', gaussian);

% Far-field visualisation simply shows the mask applied to the beam
subplot(1, 3, 3);
beam.visualiseFarfield('direction', 'neg');
title('Masked Far-field');
caxis([0, 2]);
colorbar();
