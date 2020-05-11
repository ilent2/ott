% Demonstrate different scattering methods for force calculation
%
% This example shows some of the different scattering methods included
% in the toolbox for simulating optical scattering.

% Add toolbox path
addpath('../');

%% Describe system

% Describe shape
radius = 1.0;
shape = ott.shapes.Sphere(radius);
index_relative = 1.2;

% Describe beam
% By using an abstract beam we let the method choose an appropriate
% approximation for the fields.
waist = 1.0;
beam = ott.beam.abstract.Gaussian(waist);

% Range of axial positions for graph
z = [0;0;1].*linspace(-10, 10, 80);

%% Calculate force using shape-induced force method

tic

particle = ott.scat.shapeforce.Shape(shape, index_relative);
f_shapeforce = particle.force(beam, 'position', z);

disp(['Shape-induced force took ' num2str(toc()) ' seconds']);

%% Calculate force using geometric optics

tic

particle = ott.scat.geometric.Shape(shape, index_relative);
f_geometric = particle.force(baem, 'position', z);

disp(['Geometric optics took ' num2str(toc()) ' seconds']);

%% T-matrix method

tic

particle = ott.scat.vswf.Mie(shape.radius, index_relative);
f_tmatrix = particle.force(beam, 'position', z);

disp(['T-matrix took ' num2str(toc()) ' seconds']);

%% Dipole approximation

tic

particle = ott.scat.dipole.Dielectric(shape.radius, index_relative);
f_dipole = particle.force(beam, 'position', z);

disp(['Dipole took ' num2str(toc()) ' seconds']);

%% Local interpolation

tic

% Only use a few samples for interpolation (to demonstrate capabilities)
% Don't down sample in actual simulations
sampling = 5;

particle = ott.scat.interp.Force(z(:, 1:sampling:end), ...
    f_tmatrix(:, 1:sampling:end));
f_interp = particle.force('position', z);

disp(['Local interpolation took ' num2str(toc()) ' seconds']);

%% Global interpolation (harmonic model)

tic

particle = ott.scat.interp.harmonic.Force(z, f_tmatrix);
f_harmonic = particle.force('position', z);

disp(['Global interpolation took ' num2str(toc()) ' seconds']);

%% Generate figure comparing different approximations

figure();
plot(z(3, :), [f_shapeforce(3, :); f_geometric(3, :); ...
    f_tmatrix(3, :); f_dipole(3, :), f_interp(3, :), f_harmonic(3, :)]);
title('Scattering Methods');
legend({'shapeforce', 'geometric', 'tmatrix', 'dipole', 'interp', 'harmonic'});
xlabel('Axial position');
xlabel('Force position');

