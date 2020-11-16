% ottNeuralNetwork.m -- Train and use a neural network for simulation
%
% This example shows how a neural network can be trained on simulated
% force (or torque) data.  The trained neural network can then be used
% for fast dynamics simulations or for sharing/distribution.
% This example shows a overly simplified case with a small neural network,
% however the same approach could be extended to a much more complicated
% problem.  The example requires the Matlab neural network toolbox.
%
% This approach is similar to interpolation, except, unlike interpolation,
% the neural network can be easily scaled up to include multiple inputs
% and outputs, and the resulting network often has a much smaller memory
% footprint (making it easy to share).  This idea is described in
%
%   Lenton et al., Machine learning reveals complex behaviours in
%   optically trapped particles.  MLST, 2020
%   https://doi.org/10.1088/2632-2153/abae76
%
% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% Add toolbox to path (uncomment this line if OTT is not already on the path)
%addpath('../');

%% Setup the incident beam and particle

index_medium = 1.33;    % Refractive index of medium (water)
index_particle = 1.59;  % Refractive index of particle (polystyrene)

wavelength0 = 1064.0e-9;   % Wavelength in vacuum [m]

radius = 1.0e-6;      % Radius of particle [m]

NA = 1.2;       % Numerical aperture of beam
lmode = 10;     % Azimuthal mode of LG beam

beam = ott.beam.LaguerreGaussian.FromNa(NA, 'lmode', lmode, ...
    'wavelength0', wavelength0, 'index_medium', index_medium, ...
    'power', 1.0);

sphere = ott.shape.Sphere(radius);
particle = ott.particle.Fixed.FromShape(sphere, ...
    'index_relative', index_particle ./ index_medium, ...
    'wavelength0', beam.wavelength0);

%% Calculate the optical force
% Takes about 1 minute with 100x50 grid.
%
% This example assumes the beam-particle system is axially symmetric
% However a similar approach could be done for any arbitrary particle/beam.

z = linspace(-5, 5, 100)*1e-6;
r = linspace(0, 5, 50)*1e-6;
[rr, zz] = meshgrid(r, z);

tic

fxyz = -beam.force(particle, 'position', {rr, 0*rr, zz});

disp(['Force calculation took : ' num2str(toc()) ' seconds']);

%% Visualise the forces we calculated

figure();
subplot(1, 2, 1);
imagesc(r, z, squeeze(fxyz(1, :, :)));
hold on;
contour(r, z, squeeze(fxyz(3, :, :)), [0, 0], 'w--');
contour(r, z, squeeze(fxyz(1, :, :)), [0, 0], 'w-');
hold off;
title('Radial force');
axis image;
xlabel('radial osition [m]');
ylabel('axial osition [m]');

subplot(1, 2, 2);
imagesc(r, z, squeeze(fxyz(3, :, :)));
hold on;
contour(r, z, squeeze(fxyz(3, :, :)), [0, 0], 'w--');
contour(r, z, squeeze(fxyz(1, :, :)), [0, 0], 'w-');
hold off;
title('Axial force');
axis image;
xlabel('radial osition [m]');
ylabel('axial osition [m]');

%% Train a network with this information
% The following code trains a shallow neural network.  This takes about
% 1-2 minutes for a [16,16,16] network.  It is useful to scale the force
% and position to be around unity.  Adding an offset to the force can
% sometimes be useful too.
%
% This gives decent accuracy, high enough for most dynamics simulations,
% however in some cases this won't be accurate enough.  Consider using
% different network architectures or different structuring for the grid.
%
% For this 2-D problem, interpolation would also work pretty well.

% Scale position/force values so they are approximately unity
posScale = 5e-6;
forceScale = 20e-11;
forceOffset = 0.5;
X = [rr(:), zz(:)].'./posScale;
Y = fxyz(1:3, :)./forceScale - forceOffset;

net = feedforwardnet([16, 16, 16]);
net = train(net,X,Y);
% net = train(net,X,Y,'useGPU','yes');   % Faster with a GPU
predY = net(X);
perf = perform(net,Y,predY)

% Save the data (optional)
% save('network.mat', 'net', 'beam', 'tmatrix', '-v7.3');

% Save the network (for compatibility without the NN-toolbox)
% If you wish to share the generated NN, use this option.
% genFunction(net,'LgBeamNetwork','MatrixOnly','yes');

%% Generate another dataset for validation

v_xyz = (rand(3, 1000).*[5; 0; 10] - [0; 0; 5]).*1e-6;

tic
v_fxyz = -beam.force(particle, 'position', v_xyz);
disp(['Force calculation (OTT) took : ' num2str(toc()) ' seconds']);

tic
predY = (net(v_xyz([1, 3], :)./posScale)+forceOffset).*forceScale;
% (LgBeamNetwork(v_xyz([1, 3], :)./posScale)+forceOffset).*forceScale;
disp(['Force calculation (NN) took : ' num2str(toc()) ' seconds']);
perf = perform(net,v_fxyz(1:3, :),net(v_xyz([1, 3], :)./posScale))

% Calculate a error (perhaps there are better ones we could use)
mae = mean(abs((v_fxyz([1,3],:) - predY([1,3],:))), 2)

