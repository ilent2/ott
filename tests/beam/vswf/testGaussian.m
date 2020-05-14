function tests = testGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructDefault(testCase)

  beam = ott.beam.vswf.Gaussian();

  testCase.verifyEqual(beam.waist, 1.0, 'waist');
end

function testConstructOptional(testCase)

  waist = 0.1;
  beam = ott.beam.vswf.Gaussian(waist);

  testCase.verifyEqual(beam.waist, waist, 'waist');

end

function testConstructArray(testCase)

  waist = [0.1, 0.2];
  beam = ott.beam.vswf.Gaussian(waist);

  testCase.verifyEqual(beam.waist, waist, 'waist');

end

