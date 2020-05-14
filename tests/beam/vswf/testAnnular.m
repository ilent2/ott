function tests = testAnnular
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructDefault(testCase)

  beam = ott.beam.vswf.Annular();

  testCase.verifyEqual(beam.angle_start, angle_start, 'start');
  testCase.verifyEqual(beam.angle_end, angle_end, 'end');

end

function testConstructOptional(testCase)

  angle_start = 0;
  angle_end = 0.1;

  beam = ott.beam.vswf.Annular(angle_start, angle_end);

  testCase.verifyEqual(beam.angle_start, angle_start, 'start');
  testCase.verifyEqual(beam.angle_end, angle_end, 'end');

end

