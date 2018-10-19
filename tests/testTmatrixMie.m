function tests = tmatrixmie
  tests = functiontests(localfunctions);
end

function testConstruct(testCase)

  addpath('../');

  % Simple sphere
  Tsimple = ott.TmatrixMie(1.0, 'index_relative', 1.2);

  % Layered sphere
  Tsimple = ott.TmatrixMie([0.5, 1.0], 'index_relative', [1.4, 1.2]);

end

function testShrink(testCase)

  % A test case to check if shrinking the Nmax after
  % constructing the layered sphere produces a similar force/torque

  addpath('../');

  % Generate the layered spheres
  T0 = ott.TmatrixMie([0.5, 1.0], 'index_relative', [1.4, 1.2], ...
      'shrink', false);
  T1 = ott.TmatrixMie([0.5, 1.0], 'index_relative', [1.4, 1.2], ...
      'shrink', true);

  % Generate a beam for testing
  beam = ott.BscPmGauss('power', 1.0);

  f0 = ott.forcetorque(beam, T0);
  f1 = ott.forcetorque(beam, T1);

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  testCase.verifyThat(f1, IsEqualTo(f0, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Shrinking layered T-matrix does not work');

end

function testFields(testCase)

  % Test internal fields for non-contrast particle
  T = ott.TmatrixMie(0.5, 'index_relative', 1.0, 'internal', true);
  beam = ott.BscPmGauss('power', 1.0);
  
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
  T = ott.TmatrixMie(0.5, 'index_relative', 1.0, 'internal', false);
  beam = ott.BscPmGauss('power', 1.0);
  
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
  Te = ott.TmatrixMie(R, 'index_relative', nrel, 'internal', false);
  Ti = ott.TmatrixMie(R, 'index_relative', nrel, 'internal', true);
  beam = ott.BscPmGauss('power', 1.0);
  
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
