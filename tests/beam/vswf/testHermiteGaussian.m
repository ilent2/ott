function tests = testHermiteGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructDefault(testCase)

  beam = ott.beam.vswf.HermiteGaussian();

  testCase.verifyEqual(beam.waist, 1.0, 'waist');
  testCase.verifyEqual(beam.mmode, 0, 'mmode');
  testCase.verifyEqual(beam.nmode, 0, 'mmode');
end

function testConstructOptional(testCase)

  waist = 0.1;
  mmode = 1;
  nmode = 2;
  beam = ott.beam.vswf.HermiteGaussian(waist, mmode, nmode);

  testCase.verifyEqual(beam.waist, waist, 'waist');
  testCase.verifyEqual(beam.mmode, mmode, 'mmode');
  testCase.verifyEqual(beam.nmode, nmode, 'mmode');

end

