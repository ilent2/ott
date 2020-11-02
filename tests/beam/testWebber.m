function tests = testWebber
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  theta = pi/4;
  alpha = 1.0;
  parity = 'even';
  Nmax = 2;
  beam = ott.beam.Webber(theta, alpha, parity, ...
      'initial_Nmax', Nmax);

  testCase.verifyEqual(beam.alpha, alpha, 'alpha');
  testCase.verifyEqual(beam.data.Nmax, Nmax, 'Nmax');
  testCase.verifyEqual(beam.parity, parity, 'parity');
  testCase.verifyEqual(beam.theta, theta, 'theta');

end

