function tests = testInceGaussian
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructDefault(testCase)

  beam = ott.beam.vswf.InceGaussian();

  testCase.verifyEqual(beam.waist, 1.0, 'waist');
  testCase.verifyEqual(beam.parax, 0, 'parax');
  testCase.verifyEqual(beam.azim, 0, 'azum');
  testCase.verifyEqual(beam.parity, 1, 'parity');
  testCase.verifyEqual(beam.ellip, 1, 'ellip');

end

function testConstructOptional(testCase)

  waist = 0.1;
  parax = 1;
  azim = 2;
  partity = 2;
  ellip = 2;
  beam = ott.beam.vswf.InceGaussian(waist, parax, azim, parity, ellip);

  testCase.verifyEqual(beam.waist, waist, 'waist');
  testCase.verifyEqual(beam.parax, parax, 'parax');
  testCase.verifyEqual(beam.azim, azum, 'azum');
  testCase.verifyEqual(beam.parity, parity, 'parity');
  testCase.verifyEqual(beam.ellip, ellip, 'ellip');

end

function testConstructArray(testCase)

  waist = [0.1, 0.2];
  parax = [1, 3];
  azim = [2, -2];
  partity = [2, -2];
  ellip = [2, -2];
  beam = ott.beam.vswf.InceGaussian(waist, parax, azim, parity, ellip);

  testCase.verifyEqual(beam.waist, waist, 'waist');
  testCase.verifyEqual(beam.parax, parax, 'parax');
  testCase.verifyEqual(beam.azim, azum, 'azum');
  testCase.verifyEqual(beam.parity, parity, 'parity');
  testCase.verifyEqual(beam.ellip, ellip, 'ellip');

end

