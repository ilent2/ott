function tests = testPlaneWave
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testFromDirection(testCase)

  dir = [0;0;1];
  pol = [1;0;0];
  Nmax = 20;
  beam = ott.bsc.PlaneWave.FromDirection(Nmax, dir, pol);

  testCase.verifyClass(beam, 'ott.bsc.PlaneWave');
  testCase.verifyEqual(beam.Nmax, Nmax, 'Nmax');
  testCase.verifyEqual(beam.direction, dir, 'direction');

end

function testTranslateZ(testCase)
  % Test translation in comparision to normal bsc translation

  Nmax = 40;
  dir = [0;0;1];
  pol = [1;0;0];
  beam = ott.bsc.PlaneWave.FromDirection(Nmax, dir, pol);

  % Calculate translation with Annular
  tbeam = beam.translateZ(1);

  % Calculate translation with Bsc
  bsc = ott.bsc.Bsc(beam);
  tbsc = bsc.translateZ(1);

  % Compare coefficients
  ab1 = tbeam.getCoefficients(1:ott.utils.combined_index(10, 10));
  ab2 = tbsc.getCoefficients(1:ott.utils.combined_index(10, 10));
  testCase.verifyEqual(ab1, ab2, 'RelTol', 1.0e-6);

end

function testRotate(testCase)

  Nmax = 20;
  dir = [0;0;1];
  pol = [1;0;0];
  beam = ott.bsc.PlaneWave.FromDirection(Nmax, dir, pol);

  beam = beam.rotateY(pi/2);
  testCase.verifyEqual(beam.direction, [1;0;0], ...
      'AbsTol', 1.0e-12, 'new direction');

end
