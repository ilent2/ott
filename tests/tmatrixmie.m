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
