function tests = testGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  waist = 1.0;
  beam = ott.beam.abstract.Gaussian(waist);
  testCase.verifyEqual(beam.waist, waist);
end
