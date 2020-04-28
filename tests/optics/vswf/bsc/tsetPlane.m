function tests = testBscPlane
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../../');
end

function testConstruct(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import ott.optics.vswf.*;
  beam = bsc.Plane([0.0, pi/4], [0.0, 0.0], 'radius', 1.0);

  testCase.verifyThat(beam.Nbeams, IsEqualTo(2), ...
    'Incorrect number of beams stored');

  nmax = 12;
  beam = bsc.Plane([0.0, pi/4], [0.0, 0.0], 'Nmax', nmax);
  testCase.verifyThat(beam.Nmax, IsEqualTo(nmax), ...
    'Incorrect Nmax for beam');

end

function testTranslate(testCase)

  import ott.optics.vswf.*;
  beam = bsc.Plane([0.0, pi/4], [0.0, 0.0], 'radius', 1.0);
  dz = 0.5;

  tbeam1 = beam.translateZ(dz);

  import matlab.unittest.constraints.IsEqualTo;
  testCase.verifyThat(tbeam1.Nbeams, IsEqualTo(2), ...
    'Incorrect number of beams after translation');

  [~, A, B] = beam.translateZ(dz);

  tbeam2 = beam.translate(A, B);

  tbeam3 = beam.translateZ(dz, 'Nmax', beam.Nmax - 5);

end

function testMultiple(testCase)
% Combined beams should produce the same output as individual beams

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  import ott.optics.vswf.*;
  tol = 1.0e-14;
  Nmax = 10;
  
  theta1 = 0.2*pi/2;
  theta2 = 0.6*pi/2;
  phi1 = 0.1*pi/2;
  phi2 = 0.3*pi/2;

  beam1 = bsc.Plane(theta1, phi1, 'Nmax', Nmax);
  beam2 = bsc.Plane(theta2, phi2, 'Nmax', Nmax);
  beamc = bsc.Plane([theta1, theta2], [phi1, phi2], 'Nmax', Nmax);
  
  
  testCase.verifyThat(beamc.beam(1).getCoefficients, IsEqualTo(beam1.getCoefficients, ...
    'Within', AbsoluteTolerance(tol)), ...
    'First beam doesnt match');
  
  testCase.verifyThat(beamc.beam(2).getCoefficients, IsEqualTo(beam2.getCoefficients, ...
    'Within', AbsoluteTolerance(tol)), ...
    'Second beam doesnt match');

end

