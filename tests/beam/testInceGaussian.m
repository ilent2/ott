function tests = testInceGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  waist = 1.0e-6;
  lmode = 1;
  porder = 3;
  ellip = 3;
  parity = 'odd';
  beam = ott.beam.InceGaussian(waist, lmode, porder, ellip, parity);

  testCase.verifyEqual(beam.waist, waist, 'w');
  testCase.verifyEqual(beam.lmode, lmode, 'l');
  testCase.verifyEqual(beam.porder, porder, 'p');
  testCase.verifyEqual(beam.ellipticity, ellip, 'e');
  testCase.verifyEqual(beam.parity, parity, 'parity');

  testCase.verifyEqual(beam.power, 1.0, 'power');
  testCase.verifyEqual(beam.mapping, 'sin', 'mapping');
  testCase.verifyEqual(beam.polfield, [1, 1i], 'polfield');
  testCase.verifyEqual(beam.polbasis, 'cartesian', 'polbasis');

end

function testFromNa(testCase)

  NA = 0.9;
  beam = ott.beam.InceGaussian.FromNa(NA);

end

