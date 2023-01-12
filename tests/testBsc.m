function tests = bsc
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)

  % Ensure the ott package is in our path
  addpath('../');

  % Generate a gaussian beam for testing
  testCase.TestData.beam = ott.BscPmGauss();

  % Generate a empty beam for testing
  testCase.TestData.emtpyBeam = ott.Bsc([], [], ...
      'incoming', 'incident');

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

  xbeam1 = beam.translateXyz([dz; dz; 0]);

  [~, AB] = beam.translateXyz([dz; dz; 0]);
  xbeam2 = AB * beam;

  [~, A, B] = beam.translateXyz([dz; dz; 0]);
  xbeam3 = beam.translate(A, B);

  [~, A, B, D] = beam.translateXyz([dz; dz; 0]);
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

function testEmptyTranslation(testCase)

  beam = testCase.TestData.emtpyBeam;
  tbeam = beam.translateZ(1.0);
  testCase.verifyEqual(tbeam.getCoefficients, beam.getCoefficients, ...
    'Empty beam should still be an empty beam');

end

function testMakeBeamVectorEmpty(testCase)
% Check that make_beam_vector works for empty beams

  import matlab.unittest.constraints.IsEqualTo;

  a = [];
  b = [];
  nn = [];
  mm = [];
  
  [a1, b1] = ott.Bsc.make_beam_vector(a, b, nn, mm);
  
  testCase.verifyThat([size(a1), size(b1)], IsEqualTo([0, 0, 0, 0]), ...
    'Size of beam vectors incorrect with implicit Nmax');
  
  [a2, b2] = ott.Bsc.make_beam_vector(a, b, nn, mm, 1);
  
  testCase.verifyThat([size(a2), size(b2)], IsEqualTo([3, 0, 3, 0]), ...
    'Size of beam vectors incorrect with explicit Nmax');
end

function testMakeBeamVectorMulti(testCase)
% Check to make sure make_beam_vector functions with multiple beams
% with the same nn and mm indices.

  a = [1; 2; 3];
  b = [4; 5; 6];
  
  nn = [1; 2; 3];
  mm = [0; 0; 0];
  
  [a1, b1] = ott.Bsc.make_beam_vector(a, b, nn, mm);
  [a2, b2] = ott.Bsc.make_beam_vector(a+6, b+6, nn, mm);
  
  [ac, bc] = ott.Bsc.make_beam_vector([a, a+6], [b, b+6], nn, mm);
  
  import matlab.unittest.constraints.IsEqualTo;
  
  testCase.verifyThat([ac(:, 1); bc(:, 1)], IsEqualTo([a1; b1]), ...
    'First beam doesn''t match');
  
  testCase.verifyThat([ac(:, 2); bc(:, 2)], IsEqualTo([a2; b2]), ...
    'Second beam doesn''t match');

end

function testSum(testCase)

  beam1 = ott.BscPmGauss('polarisation', [0, 1i]);
  beam2 = ott.BscPmGauss('polarisation', [1, 0]);
  beam3 = beam1 + beam2;
  
  beamU = beam1.append(beam2);
  testCase.verifyEqual(beamU.sum(), beam3, ...
    'beamU.sum() incorrect');
  testCase.verifyEqual(sum(beamU), beam3, ...
    'sum(beamU) incorrect');
  
  % Test array sum
  beamarr = [beam1, beam2];
  testCase.verifyEqual(sum(beamarr), beam3, ...
    'sum([beam1, beam2]) incorrect');

  % Test empty
  ebeam = testCase.TestData.emtpyBeam;
  testCase.verifyEqual(ebeam + ebeam, ebeam, ...
    'sum of empty beams should be empty beam');
  
end

function testShrinkNmax(testCase)

  ebeam = testCase.TestData.emtpyBeam;
  testCase.verifyEqual(ebeam.shrink_Nmax(), ebeam, ...
    'Empty beam should not shrink/change');

end

