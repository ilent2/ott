function tests = testPlaneWave
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructOptional(testCase)

  direction = [0;0;1];
  polarisation = [0;1;0];
  beam = ott.beam.vswf.PlaneWave('direction', direction, ...
      'polarisation', polarisation);

  testCase.verifyEqual(beam.direction, direction, 'dir');
  testCase.verifyEqual(beam.polarisation, polarisation, 'pol');

end

function testConstructDefault(testCase)

  beam = ott.beam.vswf.PlaneWave();

  testCase.verifyEqual(beam.direction, [0;0;1], 'dir');
  testCase.verifyEqual(beam.polarisation, [1;0;0], 'pol');

end

