function tests = tmatrixmie
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../../');
end

function testSimple(testCase)

  index_relative = 1.2;
  radius = 1.0;

  T = ott.optics.vswf.tmatrix.Mie(radius, 'index_relative', index_relative);
  
  % Check all getters
  testCase.verifyEqual(T.radius, 1.0, 'r');
  testCase.verifyEqual(T.permittivity_relative, 1.44, 'eps_r');
  testCase.verifyEqual(T.permeability_relative, 1.0, 'mu_r');
  testCase.verifyEqual(T.permittivity_medium, 1.0, 'eps_m');
  testCase.verifyEqual(T.permeability_medium, 1.0, 'mu_m');
  testCase.verifyEqual(T.wavelength_medium, 1.0, 'lambda');
  testCase.verifyEqual(T.index_relative, 1.2, 'n_r');
  testCase.verifyEqual(T.index_medium, 1.0, 'n_m');
  testCase.verifyEqual(T.index_particle, 1.2, 'n_p');
  testCase.verifyEqual(T.wavenumber_medium, 2*pi, 'k_m');
  testCase.verifyEqual(T.wavenumber_particle, 2*pi.*1.2, 'k_p');
  testCase.verifyEqual(T.permittivity_particle, 1.44, 'eps_particle');
  testCase.verifyEqual(T.permeability_particle, 1, 'mu_particle');
  testCase.verifyEqual(T.wavelength_particle, 1.0 ./ 1.2, 'wavelength');
  testCase.verifyEqual(T.Nmax, [12, 12], 'Nmax');

end

function testLayered(testCase)

  radii = [0.5, 1.0];
  indexes = [1.4, 1.2];

  T = ott.optics.vswf.tmatrix.Mie(radii, 'index_relative', indexes);
  
  testCase.verifyEqual(T.index_relative, indexes, ...
    'Index_relative not set correctly');
  testCase.verifyEqual(T.radius, radii, ...
    'Radius not set correctly');
end

function testShrink(testCase)

  % A test case to check if shrinking the Nmax after
  % constructing the layered sphere produces a similar force/torque

  % Generate the layered spheres
  import ott.optics.vswf.*;
  T0 = tmatrix.Mie([0.5, 1.0], 'index_relative', [1.4, 1.2], ...
      'shrink', false);
  T1 = tmatrix.Mie([0.5, 1.0], 'index_relative', [1.4, 1.2], ...
      'shrink', true);
    
  testCase.verifyNotEqual(T0.Nmax, T1.Nmax, 'Nmaxes are equal');

  % Generate a beam for testing
  beam = bsc.PmGauss('power', 1.0);
  
  % Calculate force
  f0 = beam.force(T0);
  f1 = beam.force(T1);
  
  testCase.verifyEqual(f1, f0, 'AbsTol', 1.0e-6, 'Shrinking failed');
end

function testFields(testCase)

  % Test internal fields for non-contrast particle
  T = ott.optics.vswf.tmatrix.Mie(0.5, 'index_relative', 1.0, 'internal', true);
  beam = ott.optics.vswf.bsc.PmGauss('power', 1.0);
  
  sbeam = T * beam;
  
  xyz = [0.2; 0.0; 0.0];
  
  E1 = beam.emFieldXyz(xyz);
  E2 = sbeam.emFieldXyz(xyz);

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  testCase.verifyThat(E1, IsEqualTo(E2, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Internal field has incorrect value');
    
  % Test external fields for non-contrast particle
  T = ott.optics.vswf.tmatrix.Mie(0.5, 'index_relative', 1.0, 'internal', false);
  beam = ott.optics.vswf.bsc.PmGauss('power', 1.0);
  
  sbeam = T * beam;
  
  xyz = [0.7; 0.0; 0.0];
  
  E1 = beam.emFieldXyz(xyz);
  E2 = E1 + 2*sbeam.emFieldXyz(xyz);

  testCase.verifyThat(E1, IsEqualTo(E2, ...
      'Within', AbsoluteTolerance(tol)), ...
      'External field has incorrect value');
    
  % Test internal and external fields satisfy boundary conditions
  % corss(n, Ei - Ee) = 0 and dot(n, Ei.*n_particle^2 - Ee) = 0
  R = 0.5;
  nrel = 1.2;
  Te = ott.optics.vswf.tmatrix.Mie(R, 'index_relative', nrel, 'internal', false);
  Ti = ott.optics.vswf.tmatrix.Mie(R, 'index_relative', nrel, 'internal', true);
  beam = ott.optics.vswf.bsc.PmGauss('power', 1.0);
  
  % Scatter the beams
  sbeam = Te * beam;
  ibeam = Ti * beam;
  
  npts = 10;
  xyz = rand(3, npts);
  xyz = xyz ./ sqrt(dot(xyz, xyz));
  
  % Calculate the fields
  Einc = beam.emFieldXyz(R*xyz);
  Eint = ibeam.emFieldXyz(R*xyz);
  Esca = sbeam.emFieldXyz(R*xyz);
  Eext1 = Esca + Einc;
  
  % Test fields are continous
  testCase.verifyThat(dot(xyz, Eext1 - Eint.*nrel^2), ...
      IsEqualTo(zeros(1, size(xyz, 2)), ...
      'Within', AbsoluteTolerance(tol)), ...
      'D field not continuous across boundary (inc+sca)');
  testCase.verifyThat(vecnorm(cross(xyz, Eext1 - Eint)), ...
      IsEqualTo(zeros(1, size(xyz, 2)), ...
      'Within', AbsoluteTolerance(tol)), ...
      'E field not continuous across boundary (inc+sca)');
    
% This seems to be limited by precision (big fields)
%   % Do same calculation again but in the in-out basis
%   beam.basis = 'incoming';
%   tbeam = sbeam.totalField(beam);
%   Ein = beam.emFieldXyz(R*xyz);
%   Eout = tbeam.emFieldXyz(R*xyz);
%   Eext2 = Ein + Eout;
%   
%   testCase.verifyThat(dot(xyz, Eext2 - Eint.*nrel^2), ...
%       IsEqualTo(zeros(1, size(xyz, 2)), ...
%       'Within', AbsoluteTolerance(tol)), ...
%       'D field not continuous across boundary (total)');
%   testCase.verifyThat(vecnorm(cross(xyz, Eext2 - Eint)), ...
%       IsEqualTo(zeros(1, size(xyz, 2)), ...
%       'Within', AbsoluteTolerance(tol)), ...
%       'E field not continuous across boundary (total)');
end
