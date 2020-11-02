function tests = testHermiteGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  waist = 1.0e-6;
  mmode = 1;
  nmode = 2;
  beam = ott.beam.HermiteGaussian(waist, mmode, nmode);
  testCase.verifyEqual(beam.mmode, mmode, 'mmode');
  testCase.verifyEqual(beam.nmode, nmode, 'nmode');
  testCase.verifyEqual(beam.power, 1.0, 'power');
  testCase.verifyEqual(beam.mapping, 'sin', 'mapping');
  testCase.verifyEqual(beam.polfield, [1, 1i], 'polfield');
  testCase.verifyEqual(beam.polbasis, 'cartesian', 'polbasis');

end

