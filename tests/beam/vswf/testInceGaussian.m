function tests = testInceGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructDefault(testCase)

  beam = ott.beam.vswf.InceGaussian();

  testCase.verifyEqual(beam.waist, 1.0, 'waist');
  testCase.verifyEqual(beam.porder, 1, 'parax');
  testCase.verifyEqual(beam.lmode, 0, 'azum');
  testCase.verifyEqual(beam.parity, 'even', 'parity');
  testCase.verifyEqual(beam.ellipticity, 1, 'ellip');

end

function testConstructOptional(testCase)

  waist = 0.1;
  parax = 5;
  azim = 3;
  parity = 'odd';
  ellip = 1.0;
  beam = ott.beam.vswf.InceGaussian(waist, azim, parax, parity, ellip);

  testCase.verifyEqual(beam.waist, waist, 'waist');
  testCase.verifyEqual(beam.porder, parax, 'parax');
  testCase.verifyEqual(beam.lmode, azim, 'azum');
  testCase.verifyEqual(beam.parity, parity, 'parity');
  testCase.verifyEqual(beam.ellipticity, ellip, 'ellip');

end
