function tests = bscplane
  tests = functiontests(localfunctions);
end

function testConstruct(testCase)

  addpath('../');
  import matlab.unittest.constraints.IsEqualTo;
  beam = ott.BscPlane([0.0, pi/4], [0.0, 0.0], 'radius', 1.0);

  testCase.verifyThat(beam.Nbeams, IsEqualTo(2), ...
    'Incorrect number of beams stored');

  nmax = 12;
  beam = ott.BscPlane([0.0, pi/4], [0.0, 0.0], 'Nmax', nmax);
  testCase.verifyThat(beam.Nmax, IsEqualTo(nmax), ...
    'Incorrect Nmax for beam');

end

function testTranslate(testCase)

  addpath('../');
  beam = ott.BscPlane([0.0, pi/4], [0.0, 0.0], 'radius', 1.0);
  dz = 0.5;

  tbeam1 = beam.translateZ(dz);

  import matlab.unittest.constraints.IsEqualTo;
  testCase.verifyThat(tbeam1.Nbeams, IsEqualTo(2), ...
    'Incorrect number of beams after translation');

  [~, A, B] = beam.translateZ(dz);

  tbeam2 = beam.translate(A, B);

  tbeam3 = beam.translateZ(dz, 'Nmax', beam.Nmax - 5);

end

