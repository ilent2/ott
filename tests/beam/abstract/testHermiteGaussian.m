function tests = testHermiteGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  waist = 1.0;
  mmode = 2;
  nmode = 3;
  beam = ott.beam.abstract.HermiteGaussian(waist, mmode, nmode);
  testCase.verifyEqual(beam.waist, waist, 'w');
  testCase.verifyEqual(beam.mmode, mmode, 'm');
  testCase.verifyEqual(beam.nmode, nmode, 'n');
end
