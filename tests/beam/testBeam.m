function tests = testBeam
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testProperties(testCase)

  beam = ott.beam.Gaussian();

  testCase.verifyEqual(beam.speed0, 3e8, 'speed0');

end

