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
  
  % Test fields are correct in near-field
  Exyz = beam.efieldRtp([0;0;0]).vxyz;
  testCase.verifyEqual(abs(Exyz), abs(pol), 'AbsTol', 1e-15, 'field direction');
  
  % Test fields are correct in far-field
  rtp = ott.utils.xyz2rtp([dir, -dir]);
  Exyz = beam.efarfield(rtp, 'basis', 'outgoing').vxyz;
  testCase.verifyEqual(abs(Exyz)./abs(Exyz(1)), [pol, 0*pol], ...
      'AbsTol', 1e-15, 'far-field direction');

end

function testZeroNmax(testCase)

  dir = randn(3, 2);
  pol = randn(3, 2);
  Nmax = 0;
  beam = ott.bsc.PlaneWave.FromDirection(Nmax, dir, pol);
  
  testCase.verifySize(beam, [1, size(dir, 2)], 'size');
  testCase.verifyEqual([beam.Nmax], [0, 0], 'nmax');
end

function testSetPower(testCase)

  Nmax = 40;
  dir = [0;0;1];
  pol = [1;0;0];
  beam = ott.bsc.PlaneWave.FromDirection(Nmax, dir, pol);
  
  beam.power = 2;
  testCase.verifyEqual(beam.power, 2, 'RelTol', 1e-15, 'set power');

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
    
  target = ott.bsc.PlaneWave.FromDirection(Nmax, [1;0;0], [0;0;-1]);
  testCase.verifyEqual(beam.getCoefficients(), target.getCoefficients(), ...
      'RelTol', 1.0e-14, 'AbsTol', 1.0e-13, 'beam data');

end
