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
  testCase.verifyEqual(beam.isGaussian, false, 'isGauss');
  
  % Check beam direction
  moment = beam.intensityMoment();
  testCase.verifyGreaterThan(moment(3), 0, 'beam should go in +z');

end

function testFromNa(testCase)

  NA = 0.9;
  beam = ott.beam.HermiteGaussian.FromNa(NA);
  testCase.verifyLessThan(beam.waist, beam.wavelength*10, 'waist may be wrong');

end

