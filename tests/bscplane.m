function tests = bscplane
  tests = functiontests(localfunctions);
end

function testConstruct(testCase)

  addpath('../');
  beam = ott.BscPlane([0.0, pi/4], [0.0, 0.0], 'radius', 1.0);
  
  import matlab.unittest.constraints.IsEqualTo;
  testCase.verifyThat(beam.beams, IsEqualTo(2), ...
    'Incorrect number of beams stored');

end

function testTranslate(testCase)

  addpath('../');
  beam = ott.BscPlane([0.0, pi/4], [0.0, 0.0], 'radius', 1.0);
  dz = 0.5;

  tbeam1 = beam.translateZ(dz);
  
  import matlab.unittest.constraints.IsEqualTo;
  testCase.verifyThat(tbeam1.beams, IsEqualTo(2), ...
    'Incorrect number of beams after translation');

  [~, A, B] = beam.translateZ(dz);

  tbeam2 = beam.translateZ(A, B);

end

