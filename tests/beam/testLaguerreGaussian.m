function tests = testLaguerreGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  waist = 1.0e-6;
  lmode = 2;
  pmode = 3;
  beam = ott.beam.LaguerreGaussian(waist, lmode, pmode);

  testCase.verifyEqual(beam.waist, waist, 'w');
  testCase.verifyEqual(beam.lmode, lmode, 'l');
  testCase.verifyEqual(beam.pmode, pmode, 'p');
  testCase.verifyEqual(beam.isGaussian, false, 'isGauss');

  testCase.verifyEqual(beam.power, 1.0, 'power');
  testCase.verifyEqual(beam.mapping, 'sin', 'mapping');
  testCase.verifyEqual(beam.polfield, [1, 1i], 'polfield');
  testCase.verifyEqual(beam.polbasis, 'cartesian', 'polbasis');
  
  % Check beam direction
  moment = beam.intensityMoment();
  testCase.verifyGreaterThan(moment(3), 0, 'beam should go in +z');

end

function testFromNa(testCase)

  NA = 0.9;
  beam = ott.beam.LaguerreGaussian.FromNa(NA);
  testCase.verifyLessThan(beam.waist, beam.wavelength*10, 'waist may be wrong');

end
