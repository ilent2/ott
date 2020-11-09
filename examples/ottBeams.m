% ottBeam.m -- Shows how to create and visualise beams
%
% This example demonstrates how to create and use different kinds of
% optical tweezers toolbox (OTT) beams.  For more information, see the
% corresponding example in the OTT documentation and the `beam` package
% reference section.  An interactive live script of this example is also
% available in the ``examples/liveScripts/`` directory.
%
% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is an example from the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../');

%% Creating a beam
% OTT provides two main methods for creating beams: using the classes in
% the `ott.beam` package and using the classes in the `ott.bsc` package.
% The `ott.beam` package is intended to provide a easy to use, high-level
% interface for interacting with beams.  Internally, the `ott.beam` classes
% use the functions declared in `ott.bsc`, and in some cases, using the
% `ott.bsc` classes directly can give better run-times.  However, for most
% users, the `ott.beam` classes should be adequate.
%
% In this example, we will focus on the `ott.beam` package.

% The `beam` package provides several methods classes representing common
% beams, for example, to create a Gaussian beam, call
beam = ott.beam.Gaussian();

% This create the default Gaussian beam (see the documentation for the
% exact specification).  We can change the beam properties such as the
% medium speed and power using the beam's properties or methods, for example
beam.power = 0.1;           % Set power [Watts]
beam.speed = 3e8/1.33;      % Set speed [meters/second]

% Some properties, such as wavelength, need to be set using beam methods.
% Beams are not handle classes, as such, it is important to store the
% resulting beam object.  The following sets the wavelength by keeping the
% speed fixed:
beam = beam.setWavelength(532e-9, 'fixedSpeed'); % Set wavelength [meters]

% Alternatively, we can create the beam with the desired properties at
% the start by passing a list of named arguments to the constructor:
beam = ott.beam.Gaussian('power', 0.1, ...
    'index_medium', 1.33, 'omega', 2*pi*3e8/532e-9);

% To view our beam, we can use one of the beam visualisation methods,
% for example, the following creates a XY slice through the beam focus:
figure();
beam.visNearfield();

%% Translating and rotating beams
% Beams have a ``position`` and ``rotation`` property.  These properties
% are applied to the beam whenever the beam is used (for example, when a
% visualisation method is called).

% The position property is a 3x1 numeric vector.  To shift the beam by
% 1 wavelength in the x direction, we can directly set the position property:
beam.position = [1;0;0]*beam.wavelength;
figure();
beam.visNearfield();

% Alternatively, we can use the `translateXyz` method.  The translation
% method applies the translation on top of any existing displacement and
% returns a new copy of the beam, for example, to translate our previous beam
% along the Y-direction, we could use:
tbeam = beam.translateXyz([0;1;0]*beam.wavelength);
figure();
tbeam.visNearfield();

% Rotations are stored as 3x3 rotation matrices.  As with the ``position``
% property, we can also directly set the ``rotation`` property, however
% it is often easier to use the ``rotate*`` methods.
%
% The following rotates the beam pi/2 radians about the Y axis.
% When the beam is used, the rotation is applied before the translation.
rbeam = beam.rotateY(pi/2);
figure();
rbeam.visNearfield();

%% Creating arrays of beams
% The toolbox supports three kinds of arrays: Coherent arrays, Incoherent
% arrays, and generic arrays.  Incoherent/Coherent arrays represent beams
% which can be represented by a finite set of sub-beams.  Generic arrays
% are simply collections of multiple beams.

% To create a generic array, simply use Matlab's array syntax, for example
beams = [ott.beam.Gaussian(), ...
    ott.beam.LaguerreGaussian('lmode', 10)];

% Most operations can be applied to generic arrays.  The result is the
% same as applying the operation to each element of the array.  For example,
% to translate the array of beams
tbeams = beams.translateXyz([1;0;0]*beam.wavelength);

% Or to set the position of each beam
[tbeams.position] = deal([1;0;0]*beam.wavelength);

% Coherent and Incoherent arrays can be created using the Coherent
% and Incoherent classes, for example
cbeams = ott.beam.Coherent(beams);
figure();
cbeams.visNearfield();

%% Calculating the change in momentum between two beams
% The ``ott.beam.Beam`` class provides methods for calculating the change
% in momentum between two beams.  Although it is more common to calculate
% the force acting on a particle (see the ottForce.m example), the follow
% shows how to calculate the change in momentum between two beams.

beam1 = ott.beam.Gaussian();
beam2 = beam1.rotateY(pi/2);
force = beam1.force(beam2)

%% Creating a custom beam
% Although the toolbox has several different beams commonly used in
% optical trapping (for a complete list, see the `beam` package
% reference section), it is sometimes necessary to create a custom beam.
% The most common scenario is when modelling an SLM or the experimentally
% measured field at the back aperture of the focussing objective.  For this
% task we can use the `PmParaxial` class (for more control over the fields
% we could also use the `ott.bsc` classes).  The following example shows
% how we could model a phase-only SLM illuminated by a Gaussian-like beam:

% Generate coordinates for pattern
x = linspace(-1, 1, 20);
y = linspace(-1, 1, 20);
[X, Y] = ndgrid(x, y);

% Calculate incident field
E0 = exp(-(X.^2 + Y.^2)./4);

% Calculate SLM-like pattern
kx = 2;
phi = 2*pi*kx*x;

% Calculate field at back aperture
E = E0 .* exp(1i*phi);

% Calculate beam
beam = ott.beam.PmParaxial.InterpProfile(X, Y, E);
figure();
beam.visNearfield();

