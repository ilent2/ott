function tests = bsc
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)

  % Ensure the ott package is in our path
  addpath('../');

  % Generate a gaussian beam for testing
  testCase.TestData.beam = ott.BscPmGauss();

end

function testTranslation(testCase)
  % Check that the translation functions run

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  beam = testCase.TestData.beam;
  beam.power = 1.0;
  dz = pi/2;
  tol = 1.0e-6;

  tbeam1 = beam.translateZ(dz);

  [~, AB] = beam.translateZ(dz);
  tbeam2 = AB * beam;

  [~, A, B] = beam.translateZ(dz);
  tbeam3 = [ A B ; B A ] * beam;
  tbeam4 = beam.translate(A, B);

  % Check all beams are equal
  target = tbeam1.getCoefficients;
  testCase.verifyThat(tbeam2.getCoefficients, IsEqualTo(target, ...
    'Within', AbsoluteTolerance(tol)), ...
    'AB * beam does not match translateZ');
  testCase.verifyThat(tbeam3.getCoefficients, IsEqualTo(target, ...
    'Within', AbsoluteTolerance(tol)), ...
    '[A, B; B, A] * beam does not match translateZ');
  testCase.verifyThat(tbeam4.getCoefficients, IsEqualTo(target, ...
    'Within', AbsoluteTolerance(tol)), ...
    'beam.translate(A, B) does not match translateZ');

  xbeam1 = beam.translateXyz([dz, dz, 0]);

  [~, AB] = beam.translateXyz([dz, dz, 0]);
  xbeam2 = AB * beam;

  [~, A, B] = beam.translateXyz([dz, dz, 0]);
  xbeam3 = beam.translate(A, B);

  [~, A, B, D] = beam.translateXyz([dz, dz, 0]);
  xbeam4 = beam.translateXyz(A, B, D);

  % Check all beams are equal
  target = xbeam1.getCoefficients;
  testCase.verifyThat(xbeam2.getCoefficients, IsEqualTo(target, ...
    'Within', AbsoluteTolerance(tol)), ...
    'AB * beam does not match translateXyz');
  testCase.verifyThat(xbeam3.getCoefficients, IsEqualTo(target, ...
    'Within', AbsoluteTolerance(tol)), ...
    'beam.translate(A, B) does not match translateXyz');
  testCase.verifyThat(xbeam4.getCoefficients, IsEqualTo(target, ...
    'Within', AbsoluteTolerance(tol)), ...
    'beam.translateXyz(A, B, D) does not match translateXyz');
end

function testUnevenTranslation(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  beam = testCase.TestData.beam;
  beam.power = 1.0;
  dz = pi/2;
  tol = 1.0e-6;

  tbeam = beam.translateZ(dz, 'Nmax', beam.Nmax - 5);
  testCase.verifyThat(tbeam.Nmax, IsEqualTo(beam.Nmax - 5), ...
    'Translated beam does not have correct Nmax');

end

