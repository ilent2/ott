function tests = testGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  waist = 1.0;
  index_medium = 1.5;
  beam = ott.optics.beam.Gaussian(waist, ...
    'index_medium', index_medium);
  
  testCase.verifyEqual(beam.waist, waist);
  testCase.verifyEqual(beam.index_medium, index_medium);
  testCase.verifyEqual(beam.permittivity, index_medium.^2);
end
