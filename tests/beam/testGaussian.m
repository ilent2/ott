function tests = testGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  beam = ott.beam.Gaussian();
  testCase.verifyEqual(beam.power, 1.0, 'power');
  testCase.verifyEqual(beam.mapping, 'sin', 'mapping');
  testCase.verifyEqual(beam.polfield, [1, 1i], 'polfield');
  testCase.verifyEqual(beam.polbasis, 'cartesian', 'polbasis');
  testCase.verifyEqual(beam.isGaussian, true, 'isGauss');
  
  % Check beam direction
  moment = beam.intensityMoment();
  testCase.verifyGreaterThan(moment(3), 0, 'beam should go in +z');

end

function testFromNa(testCase)

  NA = 0.9;
  beam = ott.beam.Gaussian.FromNa(NA);
  testCase.assertInstanceOf(beam, 'ott.beam.Gaussian');
  testCase.verifyLessThan(beam.waist, beam.wavelength*10, 'waist may be wrong');
  
end

