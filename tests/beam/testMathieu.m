function tests = testMathieu
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  theta = pi/4;
  morder = 1;
  ellip = 1;
  parity = 'even';
  Nmax = 2;
  beam = ott.beam.Mathieu(theta, morder, ellip, parity, ...
      'Nmax', Nmax);

  testCase.verifyEqual(beam.parity, parity, 'parity');
  testCase.verifyEqual(beam.ellipticity, ellip, 'ellipticity');
  testCase.verifyEqual(beam.Nmax, Nmax, 'Nmax');
  testCase.verifyEqual(beam.morder, morder, 'mmode');
  testCase.verifyEqual(beam.theta, theta, 'theta');

end

