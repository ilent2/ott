function tests = testLaguerreGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  waist = 1.0;
  lmode = 1;
  pmode = 2;
  beam = ott.optics.beam.LaguerreGaussian(waist, lmode, pmode);
  testCase.verifyEqual(beam.waist, waist);
  testCase.verifyEqual(beam.lmode, lmode);
  testCase.verifyEqual(beam.pmode, pmode);
end
